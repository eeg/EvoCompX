#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "dispersal.h"


/***
 * These functions allow dispersal beyond the neighboring cells, according to
 * the specified dispersal kernel.  This makes the code a lot messier.
 ***/


/*** functions for manipulating the kernel itself ***/


/*** allocate memory for the dispersal kernel ***/
// FIXME: should free all the kernels allocated before exiting
double *allocate_kernel(int kernel_width)
{
	double *new_kernel;

	new_kernel = (double *) malloc(sizeof(double)*kernel_width);

	if (new_kernel == NULL)
	{
		printf("new_kernel malloc failed\n");
		exit(1);
	}

	return new_kernel;
}

/*** assign the central element so that the kernel elements sum to 0 ***/
void normalize_kernel(double new_kernel[], int kernel_halfwidth)
{
	int k;
	double sum;

	new_kernel[kernel_halfwidth] = 0;
	sum = 0;
	for (k=0; k<kernel_halfwidth*2+1; k++)
		sum += new_kernel[k];
	new_kernel[kernel_halfwidth] = sum * (-1);
}

/*** report the kernel ***/
void print_kernel(double kernel[], int kernel_width)
{
	int i;

	for (i=0; i<kernel_width; i++)
		printf("%f  ", kernel[i]);
	printf("\n");
}

/*** read in the kernel file ***/
int read_kernel(char *kernel_file, double **base_kernel, int *kernel_width, 
		int MAX_KERN_WIDTH)
{
	FILE *fp;
	int countme;
	int i;
	double readin[MAX_KERN_WIDTH];

	if ( (fp = fopen(kernel_file, "r")) == NULL)
	{
		printf("ERROR: kernel file not found\n");
		return -1;
	}

	for (countme=0;countme<MAX_KERN_WIDTH;countme++)
		if (fscanf(fp, "%lf", &readin[countme]) < 1) break;

	fclose(fp);

	*kernel_width = countme;

	if (*kernel_width%2 == 0)
	{
		printf("ERROR: the kernel width is not an odd integer\n");
		return -1;
	}

	*base_kernel = allocate_kernel(*kernel_width);
	for (i=0; i<*kernel_width; i++)
		(*base_kernel)[i] = readin[i];

	normalize_kernel(*base_kernel, (*kernel_width-1)/2);

	// get here if everything is fine
	return 0;
}


/*** functions for using the kernel ***/


/*** apply the no-dispersal kernel to the landscape border ***/
void apply_border_kernel(Cell space[][2], int kernel_width, int MAX_KERN_WIDTH, 
		int REAL_START, int REAL_STOP, int TOTAL_COLS, int m)
{
	int i, j;
	double *null_kernel;

	null_kernel = allocate_kernel(kernel_width);

	// make the null kernel
	for (i=0; i<MAX_KERN_WIDTH; i++)
		null_kernel[i] = 0;

	// left side
	for (i=0; i<REAL_START; i++)
		for (j=0; j<2; j++)
			space[i][j].kernel[m] = null_kernel;

	// right side
	for (i=REAL_STOP+1; i<TOTAL_COLS; i++)
		for (j=0; j<2; j++)
			space[i][j].kernel[m] = null_kernel;
}


/*** apply dispersal kernels to all cells ***/
void assign_kernels(Cell space[][2], double base_kernel[], int kernel_width, 
		int MAX_KERN_WIDTH, int TOTAL_COLS, int REAL_START, int REAL_STOP, 
		int m)
{
	int i, j, k, b;
	double *old_kernel;
	int kernel_halfwidth = (kernel_width-1)/2;

	old_kernel = allocate_kernel(kernel_width);

	/* apply the base kernel to all real cells */
	for (i=REAL_START; i<=REAL_STOP; i++)
		for (j=0; j<2; j++)
			space[i][j].kernel[m] = base_kernel;

	/* apply a null kernel to all border cells */
	apply_border_kernel(space, kernel_width, MAX_KERN_WIDTH, REAL_START, 
			REAL_STOP, TOTAL_COLS, m);


	/*** modify the kernels near the edges of space ***/

	// allocate new kernels for all affected cells
	for (i=REAL_START; i<REAL_START+kernel_halfwidth; i++)
	{
		if (space[i][0].kernel[m] == base_kernel)
		{
			space[i][0].kernel[m] = allocate_kernel(kernel_width);
			for (k=0; k<kernel_width; k++)
				(space[i][0].kernel[m])[k] = base_kernel[k];
		}
	}
	for (i=REAL_STOP; i>REAL_STOP-kernel_halfwidth; i--)
	{
		if (space[i][0].kernel[m] == base_kernel)
		{
			space[i][0].kernel[m] = allocate_kernel(kernel_width);
			for (k=0; k<kernel_width; k++)
				(space[i][0].kernel[m])[k] = base_kernel[k];
		}
	}

	// left edge
	for (i=REAL_START; i<REAL_START+kernel_halfwidth; i++)
	{
		b = i - REAL_START + 1;
		for (j=1; j<=kernel_halfwidth-b+1; j++)
		{
			(space[i][0].kernel[m])[kernel_halfwidth-b+j] += 
						(space[i][0].kernel[m])[kernel_halfwidth-b+1-j];
			(space[i][0].kernel[m])[kernel_halfwidth-b+1-j] = 0;
		}
		space[i][1].kernel[m] = space[i][0].kernel[m];
	}

	// right edge
	for (i=REAL_STOP; i>REAL_STOP-kernel_halfwidth; i--)
	{
		b = REAL_STOP - i + 1;
		for (j=1; j<=kernel_halfwidth+1-b; j++)
		{
			(space[i][0].kernel[m])[kernel_halfwidth+b-j] += 
						(space[i][0].kernel[m])[kernel_halfwidth+b+j-1];
			(space[i][0].kernel[m])[kernel_halfwidth+b+j-1] = 0;
		}
		space[i][1].kernel[m] = space[i][0].kernel[m];
	}
}


void dispersal_happens(Cell space[][2], int kernel_halfwidth[], int old, 
		int REAL_START, int REAL_STOP, int UNDEF)
{
	int new = (old+1)%2;
	int i, j, k, sp;

	/* put a copy of the old landscape into new */
	for (sp=0; sp<2; sp++)
	{
		for (i=REAL_START; i<=REAL_STOP; i++)
		{
			space[i][new].num[sp] = space[i][old].num[sp];
			space[i][new].zbar[sp] = space[i][old].zbar[sp];
			space[i][new].ztotal[sp] = space[i][old].ztotal[sp];
		}
	}

	for (sp=0; sp<2; sp++)
	{
		/* dispersal for all real cells (from i to j) */
		for (i=REAL_START; i<=REAL_STOP; i++)
		{
			k = 0;

			for (j=i-kernel_halfwidth[sp]; j<=i+kernel_halfwidth[sp]; j++)
			{
				space[j][new].num[sp] += (space[i][old].num[sp]) * 
							(space[i][old].kernel[sp][k]);
				space[j][new].ztotal[sp] += (space[i][old].num[sp]) * 
							(space[i][old].kernel[sp][k]) * 
							space[i][old].zbar[sp];
				k++;
			}
		}

		/* compute mean phenotype */
		for (i=REAL_START; i<=REAL_STOP; i++)
		{
			if (space[i][new].num[sp] == 0)
				space[i][new].zbar[sp] = UNDEF;
			else
				space[i][new].zbar[sp] = space[i][new].ztotal[sp] / 
									space[i][new].num[sp];
		}
	}
}
