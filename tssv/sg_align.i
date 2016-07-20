%module sg_align

%{
#include "sg_align.h"
%}

typedef struct {
  unsigned int distance,
               position;
} alignment;

extern alignment align(char *, char *, unsigned int);
