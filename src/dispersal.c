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
#include "dispersal.h"


/*** 
 * movement of individuals on the landscape
 * (note: will need more care if delta varies with cell or direction)
 ***/


/* exchange between cells i and j, from the perspective of i 
 * adjusts num and ztotal, but not zbar */
void nearest_neighbor(Cell space[][2], int old, int i, int j, 
                      double delta, int sp)
{
	int new = (old+1)%2;

	double n_gain, n_loss, z_gain, z_loss, a_gain, a_loss;

	n_gain = space[j][old].num[sp] * delta;
	z_gain = n_gain * space[j][old].zbar[sp];
	a_gain = n_gain * space[j][old].abar[sp];

	n_loss = space[i][old].num[sp] * delta;
	z_loss = n_loss * space[i][old].zbar[sp];
	a_loss = n_loss * space[i][old].abar[sp];

	space[i][new].num[sp] += n_gain - n_loss;
	space[i][new].ztotal[sp] += z_gain - z_loss;
	space[i][new].atotal[sp] += a_gain - a_loss;
}


/* run the exchange over the landscape */
void dispersal_happens(Cell space[][2], int old, Params *params)
{
	int new = (old+1)%2;
	int i, sp;

	/* put a copy of the original landscape into new */
	for (sp=0; sp<params->num_sp; sp++)
	{
		for (i=0; i<params->space_size; i++)
		{
			space[i][new].num[sp] = space[i][old].num[sp];
			space[i][new].zbar[sp] = space[i][old].zbar[sp];
			space[i][new].ztotal[sp] = space[i][old].ztotal[sp];
			space[i][new].abar[sp] = space[i][old].abar[sp];
			space[i][new].atotal[sp] = space[i][old].atotal[sp];
		}
	}

	for (sp=0; sp<params->num_sp; sp++)
	{
		if (params->delta[sp] > 0)
		{
			/* nearest-neighbor dispersal for all cells except edges */
			for (i=1; i<params->space_size-1; i++)
			{
				nearest_neighbor(space, old, i, i-1, params->delta[sp], sp);
				nearest_neighbor(space, old, i, i+1, params->delta[sp], sp);
			}

			/* now the two edge cells */
			nearest_neighbor(space, old, 0, 1, params->delta[sp], sp);
			nearest_neighbor(space, old, params->space_size-1, 
			                 params->space_size-2, params->delta[sp], sp);


			/* update mean phenotype, i.e. get zbar from ztotal; 
			 * same for mean breeding value */
			for (i=0; i<params->space_size; i++)
			{
				if (space[i][new].num[sp] == 0)
				{
					space[i][new].zbar[sp] = UNDEF_PHEN;
					space[i][new].abar[sp] = UNDEF_PHEN;
				}
				else
				{
					space[i][new].zbar[sp] = space[i][new].ztotal[sp] / 
										space[i][new].num[sp];
					space[i][new].abar[sp] = space[i][new].atotal[sp] / 
										space[i][new].num[sp];
				}
			}
		}
	}
}
