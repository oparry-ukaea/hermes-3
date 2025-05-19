#ifndef TEST_REACTIONS_H__
#define TEST_REACTIONS_H__

#include "gtest/gtest.h"

#include "component.hxx"
#include "test_extras.hxx" // FakeMesh
#include <bout/constants.hxx>
#include <bout/field_factory.hxx> // For generating functions

#include "../../include/amjuel_hyd_recombination.hxx"

/// Global mesh
namespace bout::globals {
extern Mesh* mesh;
} // namespace bout::globals

// The unit tests use the global mesh
using namespace bout::globals;

class ReactionTestData {
public:
  ReactionTestData(std::string reaction_str) {
    std::cout << "Reading test data for " << reaction_str << std::endl;
  }
};

static Options default_opts =
    Options({{"units", {{"eV", 1.0}, {"inv_meters_cubed", 1.0}, {"seconds", 1.0}}}});

template <typename RTYPE>
class ReactionTest : public FakeMeshFixture {

  static_assert(std::is_base_of<Component, RTYPE>(),
                "Template arg to ReactionTest must derive from Component");

protected:
  ReactionTest(std::string reaction_str)
      : reaction_str(reaction_str), data(reaction_str),
        component("test", default_opts, nullptr) {
    std::cout << "Creating reg. test for reaction " << reaction_str << std::endl;
  };
  RTYPE component;
  const ReactionTestData data;

private:
  std::string reaction_str;
};

class AmjuelHRecTest : public ReactionTest<AmjuelHydRecombinationIsotope<'d'>> {
public:
  AmjuelHRecTest() : ReactionTest("d+ + e -> d") {}
};

#endif