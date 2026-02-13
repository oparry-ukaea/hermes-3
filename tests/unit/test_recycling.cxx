
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
using RecyclingTest = FakeMeshFixture;

// class RecyclingTest : public FakeMeshFixture {
// public:
//   RecyclingTest()
//       : FakeMeshFixture(),
//         state({
//             {"units",
//              {{"eV", 1.0},
//               {"meters", 1.0},
//               {"seconds", 1.0},
//               {"inv_meters_cubed", 1e19}}},

//              {"test", {{"species", "d+"}}},

//             //  {"species",

//             //   {"d+",
//             //    {{"recycle_as", "d"},
//             //     {"target_recycle", true},
//             //     {"density", 1},
//             //     {"temperature", 1},
//             //     {"velocity", 1},
//             //     {"AA", 2}}},

//             //   {"d", {{"density", 1}, {"temperature", 1}, {"velocity", 1}, {"AA",
//             2}}}}}

//         }),
//         component("test", state, nullptr) {}
//   Options state;
//   Recycling component;
// };

TEST_F(RecyclingTest, CreateComponent) {
  Options options;
  options["units"]["eV"] = 5;           // Normalisation temperature
  options["recycling"]["species"] = ""; // No species to recycle

  Recycling component("recycling", options, nullptr);
}

// Increase velocity, pressure and temperature
// and check that recycled density and energy sources increase
TEST_F(RecyclingTest, RecycleSourceChange) {
  Options options;
  options["units"]["eV"] = 1.0;
  options["test"]["species"] = "d+";
  options["d+"]["recycle_as"] = "d";
  options["d+"]["target_recycle"] = true;

  Options state1;
  Options state2;

  state1["species"]["d+"]["AA"] = 2.0;
  state1["species"]["d+"]["density"] = 1.0;

  state1["species"]["d"]["AA"] = 2.0;
  state1["species"]["d"]["density"] = 1.0;
  state1["species"]["d"]["velocity"] = 1.0;
  state1["species"]["d"]["temperature"] = 1.0;
  state1["species"]["d"]["pressure"] = 1.0;

  state2 = state1.copy();
  state1["species"]["d+"]["velocity"] = 1.0;
  state1["species"]["d+"]["temperature"] = 1.0;
  state1["species"]["d+"]["pressure"] = 1.0;

  state2["species"]["d+"]["velocity"] = 2.0;
  state2["species"]["d+"]["pressure"] = 2.0;
  state2["species"]["d+"]["temperature"] = 2.0;

  Recycling component("test", options, nullptr);
  component.declareAllSpecies({"d+", "d"});
  component.transform(state1);
  component.transform(state2);

  Field3D density_source1 = state1["species"]["d"]["density_source"];
  Field3D density_source2 = state2["species"]["d"]["density_source"];

  Field3D energy_source1 = state1["species"]["d"]["energy_source"];
  Field3D energy_source2 = state2["species"]["d"]["energy_source"];

  for (RangeIterator r = mesh->iterateBndryUpperY(); !r.isDone(); r++) {
    for (int jz = 0; jz < mesh->LocalNz; jz++) {

      auto i = indexAt(density_source1, r.ind, mesh->yend, jz);

      ASSERT_GT(density_source2[i], density_source1[i]);
      ASSERT_GT(energy_source2[i], energy_source1[i]);
    }
  }
}