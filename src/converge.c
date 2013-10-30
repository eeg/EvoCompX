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

#include "converge.h"
#include "input.h"
#include "landscape.h"


/***
 * assess whether convergence to equilibrium has been reached
 ***/


/*** replicate one copy of space, to be used to assess convergence with a time lag ***/
void copy_converge_landscape(Cell space[][2], int old_new, Cell converge_space[], Params *params)
{
	int i, sp;

	for (i=0; i<params->space_size; i++)
	{
		for (sp=0; sp<params->num_sp; sp++)
		{
			converge_space[i].num[sp] = space[i][old_new].num[sp];
			converge_space[i].abar[sp] = space[i][old_new].abar[sp];
			/* not using these to assess convergence:
			converge_space[i].zbar[sp] = space[i][old_new].zbar[sp];
			converge_space[i].ztotal[sp] = space[i][old_new].ztotal[sp];
			converge_space[i].atotal[sp] = space[i][old_new].atotal[sp];
			*/
		}
	}
}


/*** see how much num and abar have changed ***/
double assess_convergence(Cell space[][2], int old_new, Cell converge_space[], Params *params)
{
	double max_change = 0;
	double change;
	int sp, i;
	double num[2], abar[2];

	for (sp=0; sp<params->num_sp; sp++)
	{
		for (i=0; i<params->space_size; i++)
		{
			/* 0 = old state, for comparison; 1 = current state */
			num[0]  = converge_space[i].num[sp];
			abar[0] = converge_space[i].abar[sp];
			num[1]  = space[i][old_new].num[sp];
			abar[1] = space[i][old_new].abar[sp];

			if (num[0] == 0)
			{
				if (num[1] == 0)
					change = 0;
				else
					change = 1;
				if (change > max_change)
					max_change = change;
			}
			else
			{
				change = fabs((num[0] - num[1]) / num[0]);
				if (change > max_change)
					max_change = change;
				change = fabs((abar[0] - abar[1]) / abar[0]);
				if (change > max_change)
					max_change = change;
			}
		}
	}

	return max_change;
}
