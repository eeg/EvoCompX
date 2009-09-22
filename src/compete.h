#ifndef __COMPETE_H__
#define __COMPETE_H__

#include "landscape.h"
#include "input.h"

void competition_happens(Cell space[][2], int old, Params *params);
void comp_sel(double nz_new[2], int sp, int i, double opt, Cell space[][2],
              int old, Params *p);

#endif

