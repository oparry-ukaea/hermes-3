
#include "gtest/gtest.h"

#include "fake_mesh_fixture.hxx"
#include "test_extras.hxx" // FakeMesh

#include "../../include/recycling.hxx"

/// Global mesh
namespace bout {
namespace globals {
extern Mesh* mesh;
} // namespace globals
} // namespace bout

// The unit tests use the global mesh
using namespace bout::globals;

#include <bout/field_factory.hxx> // For generating functions

// Reuse the "standard" fixture for FakeMesh
// using RecyclingTest = FakeMeshFixture;

class RecyclingTest : public FakeMeshFixture {
public:
  RecyclingTest()
      : FakeMeshFixture(),
        state({{"units",
                {{"eV", 1.0},
                 {"meters", 1.0},
                 {"seconds", 1.0},
                 {"inv_meters_cubed", 1e19}}},
               {"test", {{"species", "d+"}}},
               {"d+",
                {{"recycle_as", "d"},
                 {"target_recycle", true},
                 {"density", 1},
                 {"temperature", 1},
                 {"velocity", 1},
                 {"AA", 2}}},
               {"d", {{"density", 1}, {"temperature", 1}, {"velocity", 1}, {"AA", 2}}}}),
        component("test", state, nullptr) {}
  Options state;
  Recycling component;
};

TEST_F(RecyclingTest, CreateComponent) {
  Options state;
  state["units"]["eV"] = 5;           // Normalisation temperature
  state["recycling"]["species"] = ""; // No species to recycle

  Recycling component("recycling", state, nullptr);
}

// Make sure that increasing the recycling fraction
// will increase the recycling source
TEST_F(RecyclingTest, RecycleFractionChange) {
  Options state1;
  Options state2;

  state1 = state.copy();
  state2 = state.copy();

  state1["species"]["d+"]["target_recycle_multiplier"] = 0.5;
  state2["species"]["d+"]["target_recycle_multiplier"] = 1.0;

  component.declareAllSpecies({"d+", "d"});
  component.transform(state1);
  component.transform(state2);

  Field3D source1 = state1["species"]["d"]["density_source"];
  Field3D source2 = state2["species"]["d"]["density_source"];

  BOUT_FOR_SERIAL(i, source1.getRegion("RGN_NOBNDRY")) {
    ASSERT_GT(source2[i], source1[i]);
  }
}