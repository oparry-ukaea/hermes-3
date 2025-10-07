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

/// Signatures for different rate calculations.
/// N.B. one extra arg required for the mass action factor.
using OneDRateFunc = std::function<BoutReal(BoutReal, BoutReal)>;
using TwoDRateFunc = std::function<BoutReal(BoutReal, BoutReal, BoutReal)>;
using RateFuncVariant = std::variant<OneDRateFunc, TwoDRateFunc>;

// This is a workaround before CWG2518/P2593R1, taken from cppreference.com
template <RateParamsTypes>
constexpr bool dependent_false = false;

// Name used to store effective temperature (RateParamsTypes::T)
static const std::string Teff_name = "Teff";

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
  RateHelper(const Options& state, const std::vector<std::string>& reactant_names,
             const Region<Ind3D> region)
      : region(region) {

    // Compute / extract fields that are required as parameters for the rate calculations
    if constexpr (RateParamsType == RateParamsTypes::ET) {
      static_assert(dependent_false<RateParamsType>,
                    "RateParamsTypes::ET not implemented");
      // add_rate_param("e:density", get<Field3D>(state["species"]["e:density"]));
      // Field3D energy = ;
      // add_rate_param("e:energy", energy);
    } else if constexpr (RateParamsType == RateParamsTypes::nT) {
      for (auto field_id : {"e:density", "e:temperature"}) {
        add_rate_param(field_id, get<Field3D>(state["species"][field_id]));
      }
    } else if constexpr (RateParamsType == RateParamsTypes::T) {
      Field3D Teff;
      calc_Teff(state, reactant_names, Teff);
      add_rate_param("Teff", Teff);
    } else {
      // Compile-time error if any other RateParamsType enum exists
      static_assert(dependent_false<RateParamsType>, "Unhandled RateParamsType");
    }

    // Extract and store reactant densities
    std::transform(
        reactant_names.begin(), reactant_names.end(),
        std::inserter(this->reactant_densities, this->reactant_densities.end()),
        [&](const auto& reactant_name) {
          return make_pair(reactant_name,
                           get<Field3D>(state["species"][reactant_name]["density"]));
        });
  }

  void add_rate_param(const std::string& field_id, const Field3D& fld) {
    this->rate_params.insert(std::make_pair(field_id, fld));
  }

  /**
   * @brief Compute the cell-averaged reaction rate, accounting for the mass action
   * factor (product of reactant densities)
   *
   * @param rate_calc_func_variant a function that calculates the rate. Typed as
   * std::variant to easily switch between different rate parameterisations.
   * @param result a std::map containing the calculated rates and collision frequencies.
   */
  void calc_rates(const RateFuncVariant& rate_calc_func_variant,
                  std::map<std::string, Field3D>& result) {
    // Set up one collision rate per reactant, plus rate in map and return it
    std::string first_key = str_keys(this->rate_params)[0];
    result["rate"] = emptyFrom(this->rate_params[first_key]);

    std::visit(
        [this, &result](auto&& rate_calc_func) {
          auto J = result["rate"].getCoordinates()->J;
          BOUT_FOR(i, region) {

            auto yp = i.yp();
            auto ym = i.ym();
            auto Ji = J[i];

            // Calc rate_central, rate_left, rate_right at this index
            BoutReal rate_central, rate_left, rate_right;

            using RateFuncType = std::decay_t<decltype(rate_calc_func)>;
            if constexpr (std::is_same_v<RateFuncType, OneDRateFunc>) {
              if constexpr (RateParamsType == RateParamsTypes::T) {
                BoutReal T = get_rate_param(Teff_name, i);
                BoutReal TL = get_rate_param_left(Teff_name, i, ym, yp);
                BoutReal TR = get_rate_param_right(Teff_name, i, ym, yp);
                rate_central = rate_calc_func(mass_action(i), T);
                rate_left = rate_calc_func(mass_action_left(i, ym, yp), TL);
                rate_right = rate_calc_func(mass_action_right(i, ym, yp), TR);
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
                rate_central = rate_calc_func(mass_action(i), ne, Te);
                rate_left = rate_calc_func(mass_action_left(i, ym, yp), ne_left, Te_left);
                rate_right =
                    rate_calc_func(mass_action_right(i, ym, yp), ne_right, Te_right);
              } else if constexpr (RateParamsType == RateParamsTypes::ET) {
                static_assert(
                    dependent_false<RateParamsType>,
                    "RateHelper::calc_rate not set up for RateParamsTypes::ET yet");
              } else {
                throw BoutException(
                    "Unhandled RateParamsType (2D rate function being passed)");
              }
            }

            // Overall rate at this index
            result["rate"][i] = 4. / 6 * rate_central
                                + (Ji + J[ym]) / (12. * Ji) * rate_left
                                + (Ji + J[yp]) / (12. * Ji) * rate_right;
          }
        },
        rate_calc_func_variant);
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
  void calc_Teff(const Options& state, const std::vector<std::string>& reactant_names,
                 Field3D& Teff) {

    std::vector<std::string> heavy_reactant_names;
    std::copy_if(reactant_names.begin(), reactant_names.end(),
                 std::back_inserter(heavy_reactant_names),
                 [](const std::string s) { return s.compare("e") != 0; });
    Teff = 0.0;
    BoutReal Tnorm = get<BoutReal>(state["units"]["eV"]);
    for (auto& sp : heavy_reactant_names) {
      Teff += (get<Field3D>(state["species"][sp]["temperature"])
               / get<Field3D>(state["species"][sp]["AA"]))
              * Tnorm;
    }

    // Clamp values
    constexpr BoutReal Teff_min = 0.01;
    constexpr BoutReal Teff_max = 10000;
    for (const auto& i : Teff.getRegion("RGN_NOBNDRY")) {
      Teff[i] = std::clamp(Teff[i], Teff_min, Teff_max);
    }
  }

  /**
   * @brief Extract the (cell-centre) value of a rate parameter.
   *
   * @param name name of the parameter (label in state["species"])
   * @param i central index
   * @return BoutReal
   */
  BoutReal get_rate_param(const std::string& name, Ind3D i) {
    return this->rate_params[name][i];
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
    return cellLeft<hermes::Limiter>(this->rate_params[name][i],
                                     this->rate_params[name][ym],
                                     this->rate_params[name][yp]);
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
    return cellRight<hermes::Limiter>(this->rate_params[name][i],
                                      this->rate_params[name][ym],
                                      this->rate_params[name][yp]);
  }

private:
  /// region in which to calculate the rate
  const Region<Ind3D> region;
  /// Function to calculate reaction rate as a function of n_e, T_e
  RateFuncVariant rate_calc_func;

  /// Reactant densities, keyed by species name
  std::map<std::string, Field3D> reactant_densities;

  // Rate parameter fields
  std::map<std::string, Field3D> rate_params;

  /**
   * @brief Compute the product of all reactant densities at a cell centre, optionally
   * excluding \p exclude_sp.
   *
   * @param i cell index
   * @param exclude_sp optional species name to exclude from the density product.
   * @return BoutReal the density product
   */
  BoutReal density_product(Ind3D i, std::string exclude_sp = "") {
    BoutReal result = 1;
    for (const auto& [sp_name, dens] : this->reactant_densities) {
      if (sp_name.compare(exclude_sp) == 0) {
        continue;
      }
      result *= dens[i];
    }
    return result;
  }

  /**
   * @brief Compute the product of all reactant densities at the left edge of a cell,
   * optionally excluding \p exclude_sp.
   *
   * @param i central index
   * @param ym neighbour 1 index
   * @param yp neighbour 2 index
   * @param exclude_sp optional species name to exclude from the density product.
   * @return BoutReal the density product
   */
  BoutReal density_product_left(Ind3D i, Ind3D ym, Ind3D yp,
                                std::string exclude_sp = "") {
    BoutReal result = 1;
    for (const auto& [sp_name, dens] : this->reactant_densities) {
      if (sp_name.compare(exclude_sp) == 0) {
        continue;
      }
      result *= cellLeft<hermes::Limiter>(dens[i], dens[ym], dens[yp]);
    }
    return result;
  }

  /**
   * @brief Compute the product of all reactant densities at the right edge of a cell,
   * optionally excluding \p exclude_sp.
   *
   * @param i central index
   * @param ym neighbour 1 index
   * @param yp neighbour 2 index
   * @param exclude_sp optional species name to exclude from the density product.
   * @return BoutReal the density product
   */
  BoutReal density_product_right(Ind3D i, Ind3D ym, Ind3D yp,
                                 std::string exclude_sp = "") {
    BoutReal result = 1;
    for (const auto& [sp_name, dens] : this->reactant_densities) {
      if (sp_name.compare(exclude_sp) == 0) {
        continue;
      }
      result *= cellRight<hermes::Limiter>(dens[i], dens[ym], dens[yp]);
    }
    return result;
  }

  /**
   * @brief Compute the mass action factor (product of all reactant densities) at a cell.
   * centre.
   *
   * @param i cell index
   * @return BoutReal the mass action factor
   */
  BoutReal mass_action(Ind3D i) { return density_product(i); }

  /**
   * @brief Compute the mass action factor (product of all reactant densities) at the left
   * edge of a cell.
   *
   * @param i central index
   * @param ym neighbour 1 index
   * @param yp neighbour 2 index
   * @return BoutReal the mass action factor
   */
  BoutReal mass_action_left(Ind3D i, Ind3D ym, Ind3D yp) {
    return density_product_left(i, ym, yp);
  }

  /**
   * @brief Compute the mass action factor (product of all reactant densities) at the
   * right edge of a cell.
   *
   * @param i central index
   * @param ym neighbour 1 index
   * @param yp neighbour 2 index
   * @return BoutReal the mass action factor
   */
  BoutReal mass_action_right(Ind3D i, Ind3D ym, Ind3D yp) {
    return density_product_right(i, ym, yp);
  }
};

#endif