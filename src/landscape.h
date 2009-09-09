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

void record_landscape(FILE *fp_num1, FILE *fp_num2, FILE *fp_zbar1, 
                      FILE *fp_zbar2, Cell space[][2], int old_new, 
                      int space_size);

void initialize_landscape(Cell space[][2], Params *params, int space_size, 
                          int UNDEF);

#endif
