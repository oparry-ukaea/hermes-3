#pragma once
#ifndef RATE_HELPER_H
#define RATE_HELPER_H

#include <functional>
#include <variant>

#include <bout/bout_types.hxx>
#include <bout/region.hxx>

#include "component.hxx"
#include "hermes_build_config.hxx"
#include "hermes_utils.hxx"
#include "integrate.hxx"

BOUT_ENUM_CLASS(RateParamsTypes, T, ET, nT)

/// Struct to hold pre-averaged data for each cell
struct CellData {
  CellData() : centre(0), left(0), right(0) {}
  CellData(BoutReal c, BoutReal l, BoutReal r) : centre(c), left(l), right(r) {}
  BoutReal centre, left, right;
};

/// Signatures for different rate calculations.
/// N.B. one extra arg required for the mass action factor.
using OneDRateFunc = std::function<BoutReal(BoutReal, BoutReal)>;
using TwoDRateFunc = std::function<BoutReal(BoutReal, BoutReal, BoutReal)>;
using RateFuncVariant = std::variant<OneDRateFunc, TwoDRateFunc>;
using RatesMap = std::map<std::string, Field3D>;

// This is a workaround before CWG2518/P2593R1, taken from cppreference.com
template <RateParamsTypes>
constexpr bool dependent_false = false;

// Name used to store effective temperature (RateParamsTypes::T)
static const std::string Teff_name = "Teff";

/**
 * @brief Struct to encapsulate reaction rate and collision frequency calculations for a
 * number of different parameterisations.
 *
 * @tparam RateParamsType type identifying the reaction rate function parameters.
 */
template <RateParamsTypes RateParamsType>
struct RateHelper {
  /**
   * @brief Construct a new RateHelper, extracting and storing some fields from the state
   * for use later in the rate calculation.
   *
   * @param state
   * @param reactant_names vector of reactant names
   * @param rate_calc_func function with which to compute the rate from the mass action
   * factor, n_e and T_e
   * @param region the region in which to calculate the rate
   */
  RateHelper(const GuardedOptions state, const Options& units,
             const std::vector<std::string>& reactant_names, const Region<Ind3D> region)
      : region(region) {

    // Compute / extract fields that are required as parameters for the rate calculations
    if constexpr (RateParamsType == RateParamsTypes::ET) {
      static_assert(dependent_false<RateParamsType>,
                    "RateParamsTypes::ET not implemented");
      // add_rate_param("e:density", get<Field3D>(state["species"]["e:density"]));
      // Field3D energy = ;
      // add_rate_param("e:energy", energy);
    } else if constexpr (RateParamsType == RateParamsTypes::nT) {
      for (auto field_lbl : {"e:density", "e:temperature"}) {
        add_rate_param(field_lbl,
                       state["species"][field_lbl].template get_ref<Field3D>());
      }
    } else if constexpr (RateParamsType == RateParamsTypes::T) {
      calc_Teff(state, units, reactant_names, Teff_storage);
      add_rate_param("Teff", Teff_storage);
    } else {
      // Compile-time error if any other RateParamsType enum exists
      static_assert(dependent_false<RateParamsType>, "Unhandled RateParamsType");
    }

    // Extract and store reactant densities
    this->reactant_names = reactant_names;
    for (const auto& reactant_name : reactant_names) {
      this->reactant_densities[reactant_name] =
          &state["species"][reactant_name]["density"].get_ref<Field3D>();
      this->freq_labels.push_back(freq_lbl(reactant_name));
    }
  }

  /**
   * @brief Compute the cell-averaged reaction rate and collision frequencies, accounting
   * for reactant densities.
   *
   * @param rate_calc_func_variant a function that calculates the rate. Typed as
   * std::variant to easily switch between different rate parameterisations.
   * @param result a std::map containing the calculated rates and collision frequencies.
   */
  void calc_rates(const RateFuncVariant& rate_calc_func_variant, RatesMap& result,
                  bool do_averaging = true) {

    // Set up map to store results: one collision frequency per reactant + reaction rate
    std::string first_key = str_keys(this->rate_params)[0];
    result["rate"] = emptyFrom(*this->rate_params[first_key]);
    for (const std::string& reactant_name : this->reactant_names) {
      result[freq_lbl(reactant_name)] = emptyFrom(*this->rate_params[first_key]);
    }

    // Temporary storage for central,left,right vals of each property at each cell index
    std::map<std::string, CellData> cell_data;
    std::transform(result.begin(), result.end(),
                   std::inserter(cell_data, cell_data.end()),
                   [&](const auto& entry) { return make_pair(entry.first, CellData()); });

    // Populate cell_data differently according to rate function type
    std::visit(
        [this, &cell_data, &do_averaging, &result](auto&& rate_calc_func) {
          auto J = result["rate"].getCoordinates()->J;
          std::vector<std::string> prop_keys = str_keys(result);
          BOUT_FOR(i, region) {

            auto yp = i.yp();
            auto ym = i.ym();
            auto Ji = J[i];

            using RateFuncType = std::decay_t<decltype(rate_calc_func)>;
            if constexpr (std::is_same_v<RateFuncType, OneDRateFunc>) {
              if constexpr (RateParamsType == RateParamsTypes::T) {
                BoutReal T = get_rate_param(Teff_name, i);
                BoutReal TL = get_rate_param_left(Teff_name, i, ym, yp);
                BoutReal TR = get_rate_param_right(Teff_name, i, ym, yp);

                CellData mass_actions = compute_mass_actions(i, ym, yp, do_averaging);
                cell_data["rate"].centre = rate_calc_func(mass_actions.centre, T);
                if (do_averaging) {
                  cell_data["rate"].left = rate_calc_func(mass_actions.left, TL);
                  cell_data["rate"].right = rate_calc_func(mass_actions.right, TR);
                }

                // collision freqs. at centre, left, right
                for (size_t j = 0; j < this->reactant_names.size(); ++j) {
                  const std::string& reactant_name = this->reactant_names[j];
                  const std::string& lbl = freq_labels[j];
                  CellData dens_prods =
                      compute_density_products(i, ym, yp, do_averaging, reactant_name);
                  cell_data[lbl].centre = rate_calc_func(dens_prods.centre, T);
                  if (do_averaging) {
                    cell_data[lbl].left = rate_calc_func(dens_prods.left, TL);
                    cell_data[lbl].right = rate_calc_func(dens_prods.right, TR);
                  }
                }
              } else {
                throw BoutException(
                    "Unhandled RateParamsType (1D rate function being passed)");
              }
            } else if constexpr (std::is_same_v<RateFuncType, TwoDRateFunc>) {
              if constexpr (RateParamsType == RateParamsTypes::nT) {
                BoutReal ne = get_rate_param("e:density", i);
                BoutReal ne_left = get_rate_param_left("e:density", i, ym, yp);
                BoutReal ne_right = get_rate_param_right("e:density", i, ym, yp);
                BoutReal Te = get_rate_param("e:temperature", i);
                BoutReal Te_left = get_rate_param_left("e:temperature", i, ym, yp);
                BoutReal Te_right = get_rate_param_right("e:temperature", i, ym, yp);

                // reaction rates at centre, left, right
                CellData mass_actions = compute_mass_actions(i, ym, yp, do_averaging);
                cell_data["rate"].centre = rate_calc_func(mass_actions.centre, ne, Te);
                if (do_averaging) {
                  cell_data["rate"].left =
                      rate_calc_func(mass_actions.left, ne_left, Te_left);
                  cell_data["rate"].right =
                      rate_calc_func(mass_actions.right, ne_right, Te_right);
                }

                // collision freqs. at centre, left, right
                for (size_t j = 0; j < this->reactant_names.size(); ++j) {
                  const std::string& reactant_name = this->reactant_names[j];
                  const std::string& lbl = freq_labels[j];
                  CellData dens_prods =
                      compute_density_products(i, ym, yp, do_averaging, reactant_name);
                  cell_data[lbl].centre = rate_calc_func(dens_prods.centre, ne, Te);
                  if (do_averaging) {
                    cell_data[lbl].left =
                        rate_calc_func(dens_prods.left, ne_left, Te_left);
                    cell_data[lbl].right =
                        rate_calc_func(dens_prods.right, ne_right, Te_right);
                  }
                }

              } else if constexpr (RateParamsType == RateParamsTypes::ET) {
                static_assert(
                    dependent_false<RateParamsType>,
                    "RateHelper::calc_rate not set up for RateParamsTypes::ET yet");
              } else {
                throw BoutException(
                    "Unhandled RateParamsType (2D rate function being passed)");
              }
            }

            // Compute averages for each property and store in result map
            for (const std::string& prop : prop_keys) {
              if (do_averaging) {
                result[prop][i] = 4. / 6 * cell_data[prop].centre
                                  + (Ji + J[ym]) / (12. * Ji) * cell_data[prop].left
                                  + (Ji + J[yp]) / (12. * Ji) * cell_data[prop].right;
              } else {
                result[prop][i] = cell_data[prop].centre;
              }
            }
          }
        },
        rate_calc_func_variant);
  }

private:
  /// region in which to calculate the rate
  const Region<Ind3D> region;
  /// Function to calculate reaction rate as a function of n_e, T_e
  RateFuncVariant rate_calc_func;

  /// Reactant densities, keyed by species name (stored as pointers to avoid copying)
  std::map<std::string, const Field3D*> reactant_densities;

  /// Reactant names in the order provided to the constructor
  std::vector<std::string> reactant_names;

  /// Pre-computed collision frequency labels for each reactant
  std::vector<std::string> freq_labels;

  // Rate parameter fields (stored as pointers to avoid copying)
  std::map<std::string, const Field3D*> rate_params;

  // Storage for Teff when RateParamsType == T (needed to keep the object alive)
  Field3D Teff_storage;

  /**
   * @brief Store a field associated with a rate parameter
   *
   * @param field_lbl Label/tag associated with the field
   * @param fld reference to the field object
   */
  void add_rate_param(const std::string& field_lbl, const Field3D& fld) {
    this->rate_params.insert(std::make_pair(field_lbl, &fld));
  }

  /**
   * @brief Compute the effective temperature (in eV) of heavy reactants.
   *
   * @details Used to scale different isotope masses and finite neutral particle
   temperatures by using the effective temperature (Amjuel p43) T_eff = (M/M_1)T_1 +
   (M/M_2)T_2
   *
   * @param[in] state
   * @param[in] reactant_names names of all reactant species
   * @param[inout] Teff Field3D object in which to store the result
   *
   * @todo read clamp values from json?
   */
  void calc_Teff(const GuardedOptions state, const Options& units,
                 const std::vector<std::string>& reactant_names, Field3D& Teff) {

    std::vector<std::string> heavy_reactant_names;
    std::copy_if(reactant_names.begin(), reactant_names.end(),
                 std::back_inserter(heavy_reactant_names),
                 [](const std::string s) { return s.compare("e") != 0; });
    Teff = 0.0;
    BoutReal Tnorm = get<BoutReal>(units["eV"]);
    for (auto& sp : heavy_reactant_names) {
      const Field3D& temperature =
          state["species"][sp]["temperature"].template get_ref<Field3D>();
      BoutReal AA = get<BoutReal>(state["species"][sp]["AA"]);
      Teff += (temperature / AA) * Tnorm;
    }

    // Clamp values
    constexpr BoutReal Teff_min = 0.01;
    constexpr BoutReal Teff_max = 10000;
    for (const auto& i : Teff.getRegion("RGN_NOBNDRY")) {
      Teff[i] = std::clamp(Teff[i], Teff_min, Teff_max);
    }
  }

  /**
   * @brief Compute the product of all reactant densities at centre, left, and right
   * positions in a single pass, optionally excluding one species.
   *
   * @param i central index
   * @param ym neighbour 1 index
   * @param yp neighbour 2 index
   * @param do_averaging whether to compute left and right values
   * @param exclude_sp optional species name to exclude from the product
   * @return CellData struct containing centre, left, and right products
   */
  CellData compute_density_products(Ind3D i, Ind3D ym, Ind3D yp, bool do_averaging,
                                    std::string exclude_sp = "") {
    CellData result{1, 1, 1};
    for (const auto& [sp_name, dens] : this->reactant_densities) {
      if (sp_name == exclude_sp) {
        continue;
      }
      result.centre *= (*dens)[i];
      if (do_averaging) {
        result.left *= cellLeft<hermes::Limiter>((*dens)[i], (*dens)[ym], (*dens)[yp]);
        result.right *= cellRight<hermes::Limiter>((*dens)[i], (*dens)[ym], (*dens)[yp]);
      }
    }
    return result;
  }

  /**
   * @brief Compute the mass action factor (product of all reactant densities) at
   * centre, left, and right positions.
   *
   * @param i central index
   * @param ym neighbour 1 index
   * @param yp neighbour 2 index
   * @param do_averaging whether to compute left and right values
   * @return CellData struct containing centre, left, and right mass action factors
   */
  CellData compute_mass_actions(Ind3D i, Ind3D ym, Ind3D yp, bool do_averaging) {
    return compute_density_products(i, ym, yp, do_averaging);
  }

  /**
   * @brief Function to get consistent labels for collision freqency variables.
   *
   * @param reactant_name Name of the reactant associated with this collision frequency.
   * @return std::string The label
   */
  std::string freq_lbl(const std::string& reactant_name) {
    return fmt::format("{:s}:collision_frequency", reactant_name);
  }

  /**
   * @brief Extract the (cell-centre) value of a rate parameter.
   *
   * @param name name of the parameter (label in state["species"])
   * @param i central index
   * @return BoutReal
   */
  BoutReal get_rate_param(const std::string& name, Ind3D i) {
    return (*this->rate_params[name])[i];
  }

  /**
   * @brief Extract the value of a rate parameter at the left of a cell.
   *
   * @param name name of the parameter (label in state["species"])
   * @param i central index
   * @param ym neighbour 1 index
   * @param yp neighbour 2 index
   * @return BoutReal
   */
  BoutReal get_rate_param_left(const std::string& name, Ind3D i, Ind3D ym, Ind3D yp) {
    return cellLeft<hermes::Limiter>((*this->rate_params[name])[i],
                                     (*this->rate_params[name])[ym],
                                     (*this->rate_params[name])[yp]);
  }

  /**
   * @brief Extract the value of a rate parameter at the right of a cell.
   * @param name name of the parameter (label in state["species"])
   * @param i central index
   * @param ym neighbour 1 index
   * @param yp neighbour 2 index
   * @return BoutReal
   */
  BoutReal get_rate_param_right(const std::string& name, Ind3D i, Ind3D ym, Ind3D yp) {
    return cellRight<hermes::Limiter>((*this->rate_params[name])[i],
                                      (*this->rate_params[name])[ym],
                                      (*this->rate_params[name])[yp]);
  }
};

#endif
