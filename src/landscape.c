/**********************************************************************
 * Copyright 2009 Emma Goldberg
 * 
 * This file is part of EvoCompX.
 * 
 * EvoCompX is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * EvoCompX is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with EvoCompX.  If not, see <http://www.gnu.org/licenses/>.
 *********************************************************************/


#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "input.h"
#include "landscape.h"
#include "optimum.h"


/***
 * initialize and report the state of the landscape
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
				space[i][j].zbar[sp] = UNDEF_PHEN;
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
			if (fscanf(num_fp, "%lf", &space[i][0].num[sp]) != 1)
			{
				fprintf(stderr, "Error: invalid input in initial_num\n");
				exit(1);
			}
			if (fscanf(zbar_fp, "%lf", &space[i][0].zbar[sp]) != 1)
			{
				fprintf(stderr, "Error: invalid input in initial_zbar\n");
				exit(1);
			}
			space[i][0].ztotal[sp] = space[i][0].num[sp] *
								space[i][0].zbar[sp];
		}
	}

	fclose(num_fp);
	fclose(zbar_fp);
}
