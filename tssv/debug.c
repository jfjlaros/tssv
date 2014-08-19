#include <stdio.h>

#include "sg_align.h"

int main(void) {
  char read[] = "GCCAACTGTTTCCAAGGTCCCTCCCATGCATGCTGCTCTCTACAGAGGCATGTGCACAGT",
       anchor[] = "AGGTCGCTCC";

  align(read, anchor);

  return 0;
}//main
