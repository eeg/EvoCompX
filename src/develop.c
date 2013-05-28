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

#include "develop.h"
#include "optimum.h"

/***
 * development, and the action of plasticity
 ***/

/* fill in zbar and ztotal, based on abar and num */
void development_happens(Cell space[][2], int old, Params *params)
{
	int i, sp;

	for (i=0; i<params->space_size; i++)
	{
		for (sp=0; sp<params->num_sp; sp++)
		{
			space[i][old].zbar[sp] = space[i][old].abar[sp] + 
			     params->bbar[sp] * get_environment(i, params->env_slope);
			space[i][old].ztotal[sp] = space[i][old].num[sp] *
			                           space[i][old].zbar[sp];
		}
	}
}
