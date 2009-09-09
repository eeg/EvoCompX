#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "dispersal.h"

/* note: will need more care if delta varies with cell or direction */

/* exchange between cells i and j, from the perspective of i 
 * adjusts num and ztotal, but not zbar */
void nearest_neighbor(Cell space[][2], int old, int i, int j, double delta, 
                      int sp)
{
	int new = (old+1)%2;

	double n_gain, n_loss, z_gain, z_loss;

	n_gain = space[j][old].num[sp] * delta;
	z_gain = n_gain * space[j][old].zbar[sp];

	n_loss = space[i][old].num[sp] * delta;
	z_loss = n_loss * space[i][old].zbar[sp];

	space[i][new].num[sp] += n_gain - n_loss;
	space[i][new].ztotal[sp] += z_gain - z_loss;
}


/* run the exchange over the landscape */
void dispersal_happens(Cell space[][2], int old, int space_size, double delta,
                       int UNDEF)
{
	int new = (old+1)%2;
	int i, sp;

	/* put a copy of the old landscape into new */
	for (sp=0; sp<2; sp++)
	{
		for (i=0; i<space_size; i++)
		{
			space[i][new].num[sp] = space[i][old].num[sp];
			space[i][new].zbar[sp] = space[i][old].zbar[sp];
			space[i][new].ztotal[sp] = space[i][old].ztotal[sp];
		}
	}

	for (sp=0; sp<2; sp++)
	{
		/* nearest-neighbor dispersal for all cells except at either edge */
		for (i=1; i<space_size-1; i++)
		{
			nearest_neighbor(space, old, i, i-1, delta, sp);
			nearest_neighbor(space, old, i, i+1, delta, sp);
		}

		/* now the two edge cells */
		nearest_neighbor(space, old, 0, 1, delta, sp);
		nearest_neighbor(space, old, space_size-1, space_size-2, delta, sp);


		/* update mean phenotype, i.e. get zbar from ztotal */
		for (i=0; i<space_size; i++)
		{
			if (space[i][new].num[sp] == 0)
				space[i][new].zbar[sp] = UNDEF;
			else
				space[i][new].zbar[sp] = space[i][new].ztotal[sp] / 
									space[i][new].num[sp];
		}
	}
}
