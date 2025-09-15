#pragma once
#ifndef RATE_HELPER_H
#define RATE_HELPER_H

#include <functional>

#include <bout/bout_types.hxx>
#include <bout/region.hxx>

#include "component.hxx"
#include "hermes_build_config.hxx"
#include "integrate.hxx"

enum class RateParamsTypes { T, ET, nT };

// This is a workaround before CWG2518/P2593R1, taken from cppreference.com
template <RateParamsTypes>
constexpr bool dependent_false = false;

typedef std::function<BoutReal(BoutReal, BoutReal, BoutReal)> RateFunctionType;
template <RateParamsTypes RateParamsType = RateParamsTypes::nT,
          typename LimiterType = hermes::Limiter, typename IdxType = Ind3D>

struct RateHelper {
  /**
   * @brief Construct a new RateHelper, extracting and storing some fields from the state
   * for use later in the rate calculation.
   *
   * @tparam LimiterType
   * @tparam RegionType
   * @param state
   * @param reactant_names vector of reactant names
   * @param rate_calc_func function with which to compute the rate from the mass action
   * factor, n_e and T_e
   * @param region the region in which to calculate the rate
   */
  RateHelper(const Options& state, const std::vector<std::string>& reactant_names,
             RateFunctionType rate_calc_func, const Region<IdxType> region)
      : rate_calc_func(rate_calc_func), region(region) {

    // Electron temperature is always a rate parameter
    std::vector<std::string> rate_param_names = {"e:temperature"};

    if constexpr (RateParamsType == RateParamsTypes::ET) {
      static_assert(dependent_false<RateParamsType>,
                    "RateParamsTypes::ET not implemented");
    } else if constexpr (RateParamsType == RateParamsTypes::nT) {
      rate_param_names.push_back("e:density");
    } else if constexpr (RateParamsType == RateParamsTypes::T) {
      // pass
    } else {
      // Compile-time error if any other RateParamsType enum exists
      static_assert(dependent_false<RateParamsType>, "Unhandled RateParamsType");
    }

    // Extract and store rate param fields
    std::transform(rate_param_names.begin(), rate_param_names.end(),
                   std::inserter(this->rate_params, this->rate_params.end()),
                   [&](const std::string& rate_param_name) {
                     return std::make_pair(
                         rate_param_name,
                         get<Field3D>(state["species"][rate_param_name]));
                   });

    // Extract and store reactant densities
    std::transform(reactant_names.begin(), reactant_names.end(),
                   std::back_inserter(this->n_reactants),
                   [&](const std::string& reactant_name) {
                     return get<Field3D>(state["species"][reactant_name]["density"]);
                   });
  }

  /**
   * @brief Compute the cell-averaged reaction rate, accounting for the mass action
   * factor (product of reactant densities)
   *
   * @tparam LimiterType
   * @tparam IdxType
   * @return Field3D the cell-averaged reaction rate
   */
  Field3D calc_rate() {
    Field3D reaction_rate{emptyFrom(this->rate_params["e:temperature"])};
    auto J = reaction_rate.getCoordinates()->J;
    BOUT_FOR(i, region) {

      auto yp = i.yp();
      auto ym = i.ym();
      auto Ji = J[i];

      // Calc rate_central, rate_left, rate_right at this index
      BoutReal rate_central, rate_left, rate_right;
      if constexpr (RateParamsType == RateParamsTypes::T) {
        rate_central = rate_calc_func(mass_action(i), get_rate_param("e:temperature", i));
        rate_left = rate_calc_func(mass_action_left(i, ym, yp),
                                   get_rate_param_left("e:temperature", i, ym, yp));
        rate_right = rate_calc_func(mass_action_right(i, ym, yp),
                                    get_rate_param_right("e:temperature", i, ym, yp));
      } else if constexpr (RateParamsType == RateParamsTypes::ET) {
        static_assert(dependent_false<RateParamsType>,
                      "RateParamsTypes::ET rate not implemented");
      } else if constexpr (RateParamsType == RateParamsTypes::nT) {
        rate_central = rate_calc_func(mass_action(i), get_rate_param("e:density", i),
                                      get_rate_param("e:temperature", i));
        rate_left = rate_calc_func(mass_action_left(i, ym, yp),
                                   get_rate_param_left("e:density", i, ym, yp),
                                   get_rate_param_left("e:temperature", i, ym, yp));
        rate_right = rate_calc_func(mass_action_right(i, ym, yp),
                                    get_rate_param_right("e:density", i, ym, yp),
                                    get_rate_param_right("e:temperature", i, ym, yp));
      } else {
        // Compile-time error if any other RateParamsType enum exists
        static_assert(dependent_false<RateParamsType>, "Unhandled RateParamsType");
      }

      // Overall rate at this index
      reaction_rate[i] = 4. / 6 * rate_central + (Ji + J[ym]) / (12. * Ji) * rate_left
                         + (Ji + J[yp]) / (12. * Ji) * rate_right;
    }
    return reaction_rate;
  }

  /**
   * @brief Extract the (cell-centre) value of a rate parameter.
   *
   * @param name name of the parameter (label in state["species"])
   * @param i central index
   * @return BoutReal
   */
  BoutReal get_rate_param(const std::string& name, IdxType i) {
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
  BoutReal get_rate_param_left(const std::string& name, IdxType i, IdxType ym,
                               IdxType yp) {
    return cellLeft<LimiterType>(this->rate_params[name][i], this->rate_params[name][ym],
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
  BoutReal get_rate_param_right(const std::string& name, IdxType i, IdxType ym,
                                IdxType yp) {
    return cellRight<LimiterType>(this->rate_params[name][i], this->rate_params[name][ym],
                                  this->rate_params[name][yp]);
  }

private:
  /// region in which to calculate the rate
  const Region<IdxType> region;
  /// Function to calculate reaction rate as a function of n_e, T_e
  RateFunctionType rate_calc_func;

  /// Reactant densities
  std::vector<Field3D> n_reactants;

  // Rate parameter fields
  std::map<std::string, Field3D> rate_params;

  BoutReal mass_action(IdxType i) {
    BoutReal result = 1;
    for (const auto& n : n_reactants) {
      result *= n[i];
    }
    return result;
  }

  BoutReal mass_action_left(IdxType i, IdxType ym, IdxType yp) {
    BoutReal result = 1;
    for (const auto& n : n_reactants) {
      result *= cellLeft<LimiterType>(n[i], n[ym], n[yp]);
    }
    return result;
  }

  BoutReal mass_action_right(IdxType i, IdxType ym, IdxType yp) {
    BoutReal result = 1;
    for (const auto& n : n_reactants) {
      result *= cellRight<LimiterType>(n[i], n[ym], n[yp]);
    }
    return result;
  }
};

#endif