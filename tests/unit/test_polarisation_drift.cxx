#include "gtest/gtest.h"

#include "test_extras.hxx" // FakeMesh
#include "fake_mesh_fixture.hxx"

#include "../../include/polarisation_drift.hxx"
#include <bout/invert_laplace.hxx>

/// Global mesh
namespace bout{
namespace globals{
extern Mesh *mesh;
} // namespace globals
} // namespace bout

// The unit tests use the global mesh
using namespace bout::globals;

// Reuse the "standard" fixture for FakeMesh
using PolarisationDriftTest = FakeMeshFixture;

TEST_F(PolarisationDriftTest, CreateComponent) {
  Options options;

  PolarisationDrift component("test", options, nullptr);
}

