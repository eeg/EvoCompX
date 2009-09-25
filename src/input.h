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


#ifndef __INPUT_H__
#define __INPUT_H__ 

/* for undefined mean phenotypes */
#define UNDEF_PHEN -9999

#include "keyvalue.h"

typedef struct
{
	int num_sp;

	/* biology */
	double r;      /* growth rate                                            */
	double K;      /* carrying capacity, per cell                            */
	double h2;     /* heritability, h^2                                      */
	double V_s;    /* variance of stabilizing selection function, \sigma_s^2 */
	double V_p;    /* variance of phenotypic distribution, \sigma_p^2        */
	double V_u;    /* variance of competition function, \sigma_u^2           */
	double beta;   /* hybridization consideration                            */
	double delta;  /* probability of dispersal into neighboring cell         */

	/* landscape */
	int space_size;         /* number of cells, x (1-D) */
	double opt_slope;       /* slope of optimum phenotype function, theta(x) */
	const char *initial_num;        /* files with initial conditions... */
	const char *initial_zbar;       /* ...for n and zbar                */

	/* record-keeping */
	int start_t;
	int stop_t;
	int record_interval;

} Params;

Params *NewParams();
void FreeParams(Params *params);
int AcquireParams(struct KeyValue *kv, Params *parameters);
Params *GetParams(int argc, char *argv[]);

#endif
