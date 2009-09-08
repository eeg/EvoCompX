#ifndef __LANDSCAPE_H__
#define __LANDSCAPE_H__

#include "input.h"

/*** spatial units on the landscape ***/
typedef struct cell
{
	double *kernel[2];
	double num[2];
	double zbar[2];
	double ztotal[2];
} Cell;

void empty_border(Cell space[][2], int REAL_START, int REAL_STOP, 
		int TOTAL_COLS, int UNDEF);

void record_landscape(
		FILE *fp_num1, FILE *fp_num2, FILE *fp_zbar1, FILE *fp_zbar2, 
		Cell space[][2], int old_new, int REAL_START, int REAL_STOP);

void initialize_landscape(Cell space[][2], Params *params, int TOTAL_COLS, 
		int REAL_START, int REAL_STOP, int UNDEF);

#endif
