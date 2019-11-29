#ifndef __SG_ALIGN_H__
#define __SG_ALIGN_H__

typedef struct {
  int distance,
      position;
} alignment;

alignment align(char *, char *, int);

#endif
