#ifndef __LANDSCAPE_H__
#define __LANDSCAPE_H__

#include "input.h"

#define MAX_NUM_SP 20
#define MAX_SPACE_SIZE 1000

/*** spatial units on the landscape ***/
typedef struct cell
{
	double num[MAX_NUM_SP];
	double zbar[MAX_NUM_SP];
	double ztotal[MAX_NUM_SP];
} Cell;

void record_landscape(FILE **fp_num, FILE **fp_zbar, Cell space[][2],
                      Params *params, int old_new);

void initialize_landscape(Cell space[][2], Params *params);

#endif
