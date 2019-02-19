#define CATCH_CONFIG_MAIN
#include <catch.hpp>

#include "../tssv/sg_align.h"


TEST_CASE("Test", "[matrix]") {
  alignment a;
  int i;

  for (i = 0; i < 10000000; i++) {
    a = align(
      (char *)"GCCAACTGTTACCAAGGTCCCTCCCATGCATGCTGCTCCCTACAGAGGCATGTGCACAGT",
      (char *)"CTGTTTCCAAGG",
      1);
  }

  REQUIRE(a.distance == 3);
}
