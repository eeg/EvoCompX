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


/***
  * keep track of numbers and mean phenotypes of each of several species
  * local adaptation
  * density dependent growth (equal intra- and inter-specific competition)
  * one-dimensional space
  * discrete time
  * nearest-neighbor dispersal (though could be repeated within a generation)
  * no dispersal barriers
***/


#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "compete.h"
#include "dispersal.h"
#include "landscape.h"
#include "optimum.h"
#include "keyvalue.h"
#include "input.h"


int main(int argc, char *argv[])
{
	Params *params;

	/* keep two copies of the space so events can happen simultaneously */
	Cell space[MAX_SPACE_SIZE][2];
	int old, new;

	/* output files */
	FILE *fp_time, *fp_num[MAX_NUM_SP], *fp_zbar[MAX_NUM_SP];
	char str[1000];

	int i, t, sp;

	int t_steps;
	int recorded;


	/*** read in the parameter file ***/

	params = GetParams(argc, argv);
	/* exits GetParams() if parameters are insufficient */

	/* construct the matrix of competition coefficients */
	make_alpha(params);

	/* the number of generations to run */
	t_steps = params->stop_t - params->start_t + 1;


	/*** clear the landscape and place initial individuals ***/

	initialize_landscape(space, params);

	/* indices for the two copies of the landscape */
	old = 0;
	new = 1;

	/*** write the initial condition to the output files ***/

	fp_time = fopen("time.dat", "w");
	for (sp=0; sp<params->num_sp; sp++)
	{
		sprintf(str, "num%d.dat", sp+1);
		fp_num[sp] = fopen(str, "w");

		sprintf(str, "zbar%d.dat", sp+1);
		fp_zbar[sp] = fopen(str, "w");
	}

	recorded = 0;
	fprintf(fp_time, "%d\n", params->start_t);
	record_landscape(fp_num, fp_zbar, space, params, old);
	recorded++;

	printf("t = %d\n", params->start_t);


	/*** actually run the model ***/

	for (t=1; t<t_steps; t++)
	{
		/*** dispersal ***/
		dispersal_happens(space, old, params);
		/* the updated state is now in new */

		/*** competition, stabilizing selection, hybridization ***/
		competition_happens(space, new, params);
		/* the updated state is now in old, ready for the next round */

		/*** record the new state, if it's time to ***/
		if (t == recorded * params->record_interval)
		{
			fprintf(fp_time, "%d\n", t + params->start_t);
			record_landscape(fp_num, fp_zbar, space, params, old);
			recorded++;
			printf("t = %d\n", t + params->start_t);
		}
	}

	/* got most of the results now */
	fclose(fp_time);
	for (sp=0; sp<params->num_sp; sp++)
	{
		fclose(fp_num[sp]);
		fclose(fp_zbar[sp]);
	}

	/*** record the final state ***/

	fp_num[0] = fopen("num_final.dat", "w");
	fp_zbar[0] = fopen("zbar_final.dat", "w");

	for (i=0; i<params->space_size; i++)
	{
		for (sp=0; sp<params->num_sp; sp++)
			fprintf(fp_num[0], "%e\t", space[i][old].num[sp]); 
		fprintf(fp_num[0], "\n");

		for (sp=0; sp<params->num_sp; sp++)
			fprintf(fp_zbar[0], "%e\t", space[i][old].zbar[sp]); 
		fprintf(fp_zbar[0], "\n");
	}

	fclose(fp_num[0]);
	fclose(fp_zbar[0]);

	FreeParams(params);

	return 0;
}
