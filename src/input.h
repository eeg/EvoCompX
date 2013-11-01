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

/* for near-zero abundances */
#define TINY 1e-16

#define MAX_NUM_SP 20
#define MAX_SPACE_SIZE 1000

#define ERROR -1

#include "vector-sm.h"
#include "keyvalue.h"

typedef struct
{
	/* number of species */
	int num_sp;

	/* biology */
	Vector r;       /* growth rate                                           */
	Vector K;       /* carrying capacity, per cell                           */
	Vector h2;      /* heritability, h^2                                     */
	Vector V_s;     /* variance of stabilizing selection func, \sigma_s^2    */
	Vector V_p;     /* variance of phenotypic distribution, \sigma_p^2       */
	Vector V_u;     /* variance of competition function, \sigma_u^2          */
	Vector beta;    /* hybridization consideration                           */
	Vector delta;   /* probability of dispersal into neighboring cell        */
	const char *alpha_file;  /* file with matrix of competition coefficients */
	double alpha[MAX_NUM_SP][MAX_NUM_SP];  /* the actual competition matrix  */
	Vector bbar;    /* degree of plasticity                                  */

	/* landscape */
	int space_size;             /* number of cells, x (1-D)                  */
	double opt_slope;           /* slope of optimum phenotype func, theta(x) */
	double env_slope;           /* slope of environment func, epsilon(x)     */
	const char *initial_num;    /* files with initial conditions...          */
	const char *initial_abar;   /* ...for n and abar                         */

	/* record-keeping */
	long long int start_t;
	long long int stop_t;
	/* [long long allows larger than 2^32-1 (~10 digits)] */
	int record_interval;
	int converge_interval;
	double converge_tolerance;

} Params;

Params *NewParams();
void FreeParams(Params *params);
int AcquireParams(struct KeyValue *kv, Params *parameters);
Params *GetParams(int argc, char *argv[]);
int CheckParam(Vector v, int nsp, const char *msg);
void VectorizeParams(Params *p);
Vector VecPar(Vector v, int n);
void PrintParams(Params *p);

#endif
