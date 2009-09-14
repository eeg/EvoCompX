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
void record_landscape(FILE **fp_num, FILE **fp_zbar, Cell space[][2],
                      Params *params, int old_new)
{
	int i, sp;

	for (i=0; i<params->space_size; i++)
	{
		for (sp=0; sp<params->num_sp; sp++)
		{
			fprintf(fp_num[sp], "%3.3e\t", space[i][old_new].num[sp]);
			fprintf(fp_zbar[sp], "%3.3e\t", space[i][old_new].zbar[sp]);
		}
	}
	for (sp=0; sp<params->num_sp; sp++)
	{
		fprintf(fp_num[sp], "\n");
		fprintf(fp_zbar[sp], "\n");
	}
}

/*** apply the initial condition ***/
void initialize_landscape(Cell space[][2], Params *params)
{
	int i, j, sp;
	FILE *num_fp, *zbar_fp;

	/*** clear the landscape ***/

	for (i=0; i<params->space_size; i++)
	{
		for (j=0; j<2; j++)
		{
			for (sp=0; sp<params->num_sp; sp++)
			{
				space[i][j].num[sp] = 0;
				space[i][j].zbar[sp] = UNDEF;
				space[i][j].ztotal[sp] = 0;
			}
		}
	}

	/*** put in the initial individuals ***/

	num_fp = fopen(params->initial_num, "r");
	zbar_fp = fopen(params->initial_zbar, "r");

	for (i=0; i<params->space_size; i++)
	{
		for (sp=0; sp<params->num_sp; sp++)
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
