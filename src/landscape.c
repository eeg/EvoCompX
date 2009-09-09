#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "input.h"
#include "landscape.h"
#include "optimum.h"


/***
 * These functions initialize and report the state of the landscape.  
 ***/

/*** write the abundances and mean phenotypes to output files ***/
void record_landscape(FILE *fp_num1, FILE *fp_num2, FILE *fp_zbar1, 
                      FILE *fp_zbar2, Cell space[][2], int old_new, 
                      int space_size)
{
	int i;

	for (i=0; i<space_size; i++)
	{
		fprintf(fp_num1, "%3.3e\t", space[i][old_new].num[0]);
		fprintf(fp_num2, "%3.3e\t", space[i][old_new].num[1]);
		fprintf(fp_zbar1, "%3.3e\t", space[i][old_new].zbar[0]);
		fprintf(fp_zbar2, "%3.3e\t", space[i][old_new].zbar[1]);
	}

	fprintf(fp_num1, "\n");
	fprintf(fp_num2, "\n");
	fprintf(fp_zbar1, "\n");
	fprintf(fp_zbar2, "\n");
}

/*** apply the initial condition ***/
void initialize_landscape(Cell space[][2], Params *params, int space_size, 
                          int UNDEF)
{
	int i, j, sp;
	FILE *num_fp, *zbar_fp;

	/*** clear the landscape ***/
	for (i=0; i<space_size; i++)
		for (j=0; j<2; j++)
			for (sp=0; sp<2; sp++)
			{
				space[i][j].num[sp] = 0;
				space[i][j].zbar[sp] = UNDEF;
				space[i][j].ztotal[sp] = 0;
			}

	/*** put in the initial individuals ***/

	num_fp = fopen(params->initial_num, "r");
	zbar_fp = fopen(params->initial_zbar, "r");

	for (i=0; i<space_size; i++)
	{
		for (sp=0; sp<2; sp++)
		{
			fscanf(num_fp, "%lf", &space[i][0].num[sp]);
			fscanf(zbar_fp, "%lf", &space[i][0].zbar[sp]);
			space[i][0].ztotal[sp] = space[i][0].num[sp] *
								space[i][0].zbar[sp];
		}
	}

	fclose(num_fp);
	fclose(zbar_fp);
}
