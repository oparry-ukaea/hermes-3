#include "test_reactions.hxx"

TEST_F(AmjuelDRecTest, SourcesRegression) {
  generate_data(); /// TMP!
  sources_regression_test();
}