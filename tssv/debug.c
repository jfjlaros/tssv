#include <stdio.h>

#include "sg_align.h"

int main(void) {
  char read[] = "GCCAACTGTTTCCAAGGTCCCTCCCATGCATGCTGCTCTCTACAGAGGCATGTGCACAGT",
       anchor[] = "AGGTCGCTCC";
  int i;

  for (i = 0; i < 10000000; i++) {
    align(read, anchor, 1);
  }

  return 0;
}//main
