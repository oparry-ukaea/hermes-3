#include "test_reactions.hxx"

// ReactionTest::ReactionTest(std::string reaction_str)
//     : reaction_str(reaction_str), data(reaction_str) {
//   std::cout << "Creating reg. test for reaction " << reaction_str << std::endl;
// }

TEST_F(AmjuelHRecTest, CreateComponent) {
  Options options{{"units", {{"eV", 1.0}, {"inv_meters_cubed", 1.0}, {"seconds", 1.0}}}};

  AmjuelHydRecombinationIsotope<'h'> component("test", options, nullptr);
}