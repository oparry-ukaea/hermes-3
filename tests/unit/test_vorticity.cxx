#include "gtest/gtest.h"

#include "test_extras.hxx" // FakeMesh
#include "fake_solver.hxx"
#include "fake_mesh_fixture.hxx"

#include "../../include/vorticity.hxx"

/// Global mesh
namespace bout{
namespace globals{
extern Mesh *mesh;
} // namespace globals
} // namespace bout

// The unit tests use the global mesh
using namespace bout::globals;

// Reuse the "standard" fixture for FakeMesh
using VorticityTest = FakeMeshFixture;

TEST_F(VorticityTest, CreateComponent) {
  FakeSolver solver;

  Options::root()["mesh"]["paralleltransform"]["type"] = "shifted";
  Options options {{"units",
                    {{"seconds", 1.0},
                     {"Tesla", 1.0},
                     {"meters", 1.0}}}};

  Vorticity component("test", options, &solver);

  Options state = solver.getState();

  ASSERT_TRUE(state.isSet("Vort"));
}

TEST_F(VorticityTest, Transform) {
  FakeSolver solver;

  Options::root()["mesh"]["paralleltransform"]["type"] = "shifted";
  Options options {{"units",
                    {{"seconds", 1.0},
                     {"Tesla", 1.0},
                     {"meters", 1.0}}}};

  Vorticity component("test", options, &solver);

  Options state;
  component.transform(state);

  ASSERT_TRUE(state["fields"].isSet("vorticity"));
  ASSERT_TRUE(state["fields"].isSet("phi"));
}
