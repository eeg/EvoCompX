#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "input.h"
#include "landscape.h"
#include "optimum.h"


/***
 * These functions initialize and report the state of the landscape.  
 * There is a "border" around real space into which edge dispersers disappear.
 ***/


/*** remove individuals and their phenotypes from the landscape border ***/
void empty_border(Cell space[][2], int REAL_START, int REAL_STOP, int TOTAL_COLS, int UNDEF)
{
	int i, j, sp;

	// left side
	for (i=0; i<REAL_START; i++)
		for (j=0; j<2; j++)
			for (sp=0; sp<2; sp++)
			{
				space[i][j].num[sp] = 0;
				space[i][j].zbar[sp] = UNDEF;
				space[i][j].ztotal[sp] = 0;
			}

	// right side
	for (i=REAL_STOP+1; i<TOTAL_COLS; i++)
		for (j=0; j<2; j++)
			for (sp=0; sp<2; sp++)
			{
				space[i][j].num[sp] = 0;
				space[i][j].zbar[sp] = UNDEF;
				space[i][j].ztotal[sp] = 0;
			}
}

/*** write the abundances and mean phenotypes to output files ***/
void record_landscape(
		FILE *fp_num1, FILE *fp_num2, FILE *fp_zbar1, FILE *fp_zbar2, 
		Cell space[][2], int old_new, int REAL_START, int REAL_STOP)
{
	int i;

	for (i=REAL_START; i<=REAL_STOP; i++)
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
void initialize_landscape(Cell space[][2], Params *params, int TOTAL_COLS, 
		int REAL_START, int REAL_STOP, int UNDEF)
{
	int i, j, sp;
	FILE *num_fp, *zbar_fp;
	//FILE *fp; // if using useme.dat below

	/*** clear the landscape ***/
	for (i=0; i<TOTAL_COLS; i++)
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

	for (i=REAL_START; i<=REAL_STOP; i++)
	{
		for (sp=0; sp<2; sp++)
		{
			fscanf(num_fp, "%lf", &space[i][0].num[sp]);
			fscanf(zbar_fp, "%lf", &space[i][0].zbar[sp]);
			space[i][0].ztotal[sp] = space[i][0].num[sp] * space[i][0].zbar[sp];
		}
	}

	fclose(num_fp);
	fclose(zbar_fp);
/*
	fp = fopen("useme.dat", "r");
	for (i=REAL_START; i<=REAL_STOP; i++)
	{
		fscanf(fp, "%lf\t%lf\t%lf\t%lf\n", &space[i][0].num[0], &space[i][0].num[1], &space[i][0].zbar[0], &space[i][0].zbar[1]);
		space[i][0].ztotal[0] = space[i][0].num[0] * space[i][0].zbar[0];
		space[i][0].ztotal[1] = space[i][0].num[1] * space[i][0].zbar[1];
	}
	fclose(fp);
*/
}

