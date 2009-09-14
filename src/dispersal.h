#ifndef __DISPERSAL_H__
#define __DISPERSAL_H__

#include "landscape.h"
#include "input.h"

void nearest_neighbor(Cell space[][2], int old, int i, int j, 
                      double delta, int sp);
void dispersal_happens(Cell space[][2], int old, Params *params);

#endif
