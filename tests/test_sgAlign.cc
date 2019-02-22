#define CATCH_CONFIG_MAIN
#include <catch.hpp>

#include "../tssv/sgAlign.h"


TEST_CASE("Test", "[matrix]") {
  alignment a = align(
      (char *)"GCCAACTGTTACCAAGGTCCCTCCCATGCATGCTGCTCCCTACAGAGGCATGTGCACAGT",
      (char *)"CTGTTTCCAAGG",
      1);

  REQUIRE(a.distance == 1);
}
