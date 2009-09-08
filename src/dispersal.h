#ifndef __DISPERSAL_H__
#define __DISPERSAL_H__

#include "landscape.h"

double *allocate_kernel(int kernel_width);
void normalize_kernel(double new_kernel[], int kernel_halfwidth);
void print_kernel(double kernel[], int kernel_width);
int read_kernel(char *kernel_file, double **base_kernel, int *kernel_width,
		int MAX_KERN_WIDTH);

void apply_border_kernel(Cell space[][2], int kernel_width, int MAX_KERN_WIDTH,
		int REAL_START, int REAL_STOP, int TOTAL_COLS, int m);
void assign_kernels(Cell space[][2], double base_kernel[], int kernel_width,
		int MAX_KERN_WIDTH, int TOTAL_COLS, int REAL_START, int REAL_STOP,
		int m);
void dispersal_happens(Cell space[][2], int kernel_halfwidth[], int old,
		int REAL_START, int REAL_STOP, int UNDEF);

#endif
