#ifndef TEST_REACTIONS_H__
#define TEST_REACTIONS_H__

#include "gtest/gtest.h"
#include <sstream>

#include "component.hxx"
#include "test_extras.hxx" // FakeMesh
#include "bout/options_io.hxx"
#include <bout/constants.hxx>
#include <bout/field_factory.hxx> // For generating functions

#include "../../include/amjuel_hyd_recombination.hxx"
#include "bout/options_io.hxx"
#include <filesystem>

/// Global mesh
namespace bout::globals {
extern Mesh* mesh;
} // namespace bout::globals

// The unit tests use the global mesh
using namespace bout::globals;

static Options default_opts =
    Options({{"units", {{"eV", 1.0}, {"inv_meters_cubed", 1.0}, {"seconds", 1.0}}}});

template <typename RTYPE>
class ReactionTest : public FakeMeshFixture {

  static_assert(std::is_base_of<Component, RTYPE>(),
                "Template arg to ReactionTest must derive from Component");

protected:
  ReactionTest(std::string lbl, std::string reaction_str)
      : lbl(lbl), reaction_str(reaction_str), component("test", default_opts, nullptr) {
    std::cout << "Creating reg. test for " << lbl << "(" << reaction_str << ")"
              << std::endl;
  };
  std::string lbl;
  RTYPE component;

  virtual Options generate_state() = 0;

  /**
   * @brief Util to generate a 3D field with values that increase linearly along axis \p
   * axis_str (which has \p axis_ngrid elements) from \p v_min to \p v_max.
   *
   * @param v_min
   * @param v_max
   * @param axis_str
   * @param axis_ngrid
   * @return Field3D
   */
  Field3D gen_lin_field(BoutReal v_min, BoutReal v_max, std::string axis_str,
                        int axis_ngrid) {

    BoutReal axis_min = 0.0, axis_max(axis_ngrid);
    BoutReal axis_range = axis_max - axis_min;
    BoutReal v_range = v_max - v_min;
    std::stringstream expression;
    expression << v_min << " + (" << axis_str << "-" << axis_min << ")/" << axis_range
               << "*" << v_range;
    return FieldFactory::get()->create3D(expression.str(), &default_opts, mesh);
  }

  std::filesystem::path ref_data_path() {
    return std::filesystem::path(__FILE__).parent_path() / "reactions" / (lbl + ".nc");
  }

  void sources_regression_test() {
    std::filesystem::path ref_path = ref_data_path();
    Options ref_state = bout::OptionsIO::create(ref_path)->read();

    // TODO - check ref state mesh matches current mesh nx,ny,mz

    // TODO extract test input from ref_state
    Options test_state = ref_state.copy();

    // Run reaction
    component.transform(test_state);

    // // TODO
    // // Check that test_state matches ref state for all species, sources
    // ASSERT_TRUE(IsFieldEqual(get<Field3D>(test_state["species"][sp][source_type]),
    //                          ref_state["species"][sp][source_type], "RGN_NOBNDRY"));
  }

  void generate_data() {
    // Generate input state
    Options state = generate_state();

    // Run reaction
    component.transform(state);

    // Write output state
    std::filesystem::path outpath = ref_data_path();
    std::filesystem::remove(outpath);
    bout::OptionsIO::create(std::string(outpath))->write(state);

    // // Copy fields of interest
    // Options fields;
    // fields["n_e"] = state["species"]["e"]["density"].as<Field3D>(mesh);
    // fields["T_e"] = state["species"]["e"]["temperature"].as<Field3D>(mesh);
    // std::filesystem::remove(outpath);
    // bout::OptionsIO::create(std::string(outpath))->write(fields);
  }

private:
  std::string reaction_str;
};

template <typename RTYPE>
class AmjuelReactionTest : public ReactionTest<RTYPE> {
public:
  AmjuelReactionTest(std::string lbl, std::string reaction_str, std::string isotope)
      : ReactionTest<RTYPE>(lbl, reaction_str), isotope(isotope) {}

protected:
  std::string isotope;

  Options generate_state() override {
    std::string ion = isotope + "+";

    // Density and Temperature ranges (log vals)
    constexpr BoutReal logn_min = 14, logn_max = 22;
    constexpr BoutReal logT_min = -1, logT_max = 4;
    Field3D ne = this->gen_lin_field(logn_min, logn_max, "y", this->ny);
    Field3D Te = this->gen_lin_field(logT_min, logT_max, "z", this->nz);

    // Construct state
    Options state{
        {"species",
         {{"e", {{"density", ne}, {"temperature", Te}}},
          {"h", {{"AA", 1.0}, {"density", 1.0}, {"temperature", 1.0}, {"velocity", 1.0}}},
          {"h+",
           {{"AA", 1.0},
            {"charge", 1.0},
            {"density", 1.0},
            {"temperature", 1.0},
            {"velocity", 1.0}}}}}};

    return state;
  }
};
class AmjuelDRecTest : public AmjuelReactionTest<AmjuelHydRecombinationIsotope<'d'>> {
public:
  AmjuelDRecTest()
      : AmjuelReactionTest<AmjuelHydRecombinationIsotope<'d'>>("Drec", "d+ + e -> d",
                                                               "d") {}
};

#endif