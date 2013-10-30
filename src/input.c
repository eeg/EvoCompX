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
#include <stdlib.h>
#include <string.h>
#include <glib.h>

#include "input.h"
#include "vector-sm.h"
#include "keyvalue.h"
#include "landscape.h"


/***
 * user input of parameters
 ***/


/*** create a Params structure to hold the input parameter values ***/
Params *NewParams()
{
	Params *p;
	p = (Params *) malloc(sizeof(Params));

	if (p == NULL)
	{
		fprintf(stderr, "NewParams malloc failed\n");
		exit(1);
	}

	return p;
}

void FreeParams(Params *p)
{
	if (p->r != NULL)
		deleteVector(p->r);
	if (p->K != NULL)
		deleteVector(p->K);
	if (p->h2 != NULL)
		deleteVector(p->h2);
	if (p->V_s != NULL)
		deleteVector(p->V_s);
	if (p->V_p != NULL)
		deleteVector(p->V_p);
	if (p->V_u != NULL)
		deleteVector(p->V_u);
	if (p->beta != NULL)
		deleteVector(p->beta);
	if (p->delta != NULL)
		deleteVector(p->delta);
	if (p->bbar != NULL)
		deleteVector(p->bbar);
	free(p);
}


/*** read in the user's input parameter values ***/
int AcquireParams(struct KeyValue *kv, Params *p)
{
	int nsp;

	p->num_sp = getKeyValueint(kv, "num_species");
     if (p->num_sp == KV_INTERR)
	{
		fprintf(stderr, "need to specify number of species, num_species\n");
		return ERROR;
	}
     if (p->num_sp < 2 || p->num_sp > MAX_NUM_SP)
	{
		fprintf(stderr, "need 2 <= num_species <= %d\n", MAX_NUM_SP);
		return ERROR;
	}
	nsp = p->num_sp;

	/*** biology ***/

	p->r = getKeyValueVector(kv, "r");
	if (CheckParam(p->r, nsp, "growth rate, r") == ERROR)
		return ERROR;

	p->K = getKeyValueVector(kv, "K");
	if (CheckParam(p->K, nsp, "carrying capacity, K") == ERROR)
		return ERROR;

	p->h2 = getKeyValueVector(kv, "h2");
	if (CheckParam(p->h2, nsp, "heritability, h2") == ERROR)
		return ERROR;

	p->V_s = getKeyValueVector(kv, "V_s");
	if (CheckParam(p->V_s, nsp, "variance of stabilizing selection, V_s") == ERROR)
		return ERROR;

	p->V_p = getKeyValueVector(kv, "V_p");
	if (CheckParam(p->V_p, nsp, "phenotypic variance, V_p") == ERROR)
		return ERROR;

	p->V_u = getKeyValueVector(kv, "V_u");
	if (CheckParam(p->V_u, nsp, "variance of competition function, V_u") == ERROR)
		return ERROR;

	p->alpha_file = getKeyValuestring(kv, "alpha_file");
     if (p->alpha_file == 0)
	{
		fprintf(stderr, "competition matrix not specified, using 1\n");
	}

	p->beta = getKeyValueVector(kv, "beta");
     if (p->beta == 0)
	{
		/* default is no hybridization */
		p->beta = newVector(1);
		p->beta[0] = 0;
	}
	else
	{
		if (CheckParam(p->beta, nsp, "hybridization consideration, beta") == ERROR)
			return ERROR;
	}

	p->delta = getKeyValueVector(kv, "delta");
	if (CheckParam(p->delta, nsp, "dispersal probability, delta") == ERROR)
		return ERROR;

	p->bbar = getKeyValueVector(kv, "bbar");
	if (CheckParam(p->bbar, nsp, "plasticity, bbar") == ERROR)
		return ERROR;

	/*** landscape ***/

	p->space_size = getKeyValueint(kv, "space_size");
     if (p->space_size == KV_INTERR || p->space_size < 1)
	{
		fprintf(stderr, "valid space_size not specified, using 100\n");
		p->space_size = 100;
	}

	p->opt_slope = getKeyValuedouble(kv, "opt_slope");
     if (p->opt_slope == KV_FLOATERR)
	{
		fprintf(stderr, "need to specify slope of optimum phenotype, "
			    "opt_slope\n");
		return -1;
	}

	/* TODO allow different slopes for theta and epsilon, i.e., B and C */
	p->env_slope = 1;

	p->initial_num = getKeyValuestring(kv, "initial_num");
     if (p->initial_num == 0)
	{
		fprintf(stderr, "need to specify initial abundances, initial_num\n");
		return -1;
	}

	p->initial_abar = getKeyValuestring(kv, "initial_abar");
     if (p->initial_abar == 0)
	{
		fprintf(stderr, "need to specify initial breeding values, initial_abar\n");
		return -1;
	}

	/*** record-keeping ***/

	p->start_t = getKeyValueint(kv, "start_t");
     if (p->start_t == KV_INTERR)
	{
		fprintf(stderr, "valid start_t not specified, using 0\n");
		p->start_t = 0;
	}

	p->stop_t = getKeyValueint(kv, "stop_t");
     if (p->stop_t == KV_INTERR || p->stop_t < p->start_t)
	{
		fprintf(stderr, "valid stop_t not specified, using start_t + 1000\n");
		p->stop_t = p->start_t + 1000;
	}

	p->record_interval = getKeyValueint(kv, "record_interval");
     if (p->record_interval == KV_INTERR || p->record_interval <= 0)
	{
		fprintf(stderr, "valid record_interval not specified, using elapsed "
		                "time/10\n");
		p->record_interval = (p->stop_t - p->start_t)/10;
	}

	return 0;		/* returns -1 elsewhere if there's an error */
}


Params *GetParams(int argc, char *argv[])
{
	struct KeyValue *kv;
	char k[100], v[100];
	Params *p;
	int i, j;

	/* open the user's file */

	if (argc < 2)
	{
		fprintf(stderr, "\n           *** This is EvoCompX ver 0.3 ***\n"
		     "  Input parameters are required.  Please see the README and\n"
		     "  example params.in file for usage examples and options.\n\n");
		exit(1);
	}

	p = NewParams();
	kv = loadKeyValue(argv[1]);
	if (kv == 0)
	{
		fprintf(stderr, "unable to load parameters file\n");
		exit(ERROR);
	}

	/* overwrite parameter values with those specified on the command line */

     if (argc > 2) for (i = 2; i < argc; i++)
	{
		if (argv[i][0] == '=')
		{
			fprintf(stderr, "Warning -- option begins with = "
					"(be sure to use option=value, without spaces)\n");
		}
		for (j = 0; argv[i][j] != 0; j++)
			if (argv[i][j] == '=')
				argv[i][j] = ' ';
		if (sscanf(argv[i], "%s %s", k, v) != 2)
			continue;

		j = KeyValuekeyindex(kv, k);
		if (j < 0)
			KeyValueaddparm(kv, k, v);
		else
		{
			free(kv->value[j]);
			kv->value[j] = strdup(v);
		}
	}

	if (AcquireParams(kv, p) == -1)
	{
		fprintf(stderr, "unable to proceed -- "
				"not all required parameters specified\n");
		exit(ERROR);
	}

	/* convert any scalar parameters into a vector of length num_species */
	VectorizeParams(p);

	return p;
}

int CheckParam(Vector v, int nsp, const char *msg)
{
     if (v == 0)
	{
		fprintf(stderr, "need to specify %s\n", msg);
		return ERROR;
	}
	if (VectorSize(v) != 1 && VectorSize(v) != nsp)
	{
		fprintf(stderr, "invalid number of values for %s\n", msg);
		return ERROR;
	}
	return 0;
}

void VectorizeParams(Params *p)
{
	int nsp = p->num_sp;

	if (VectorSize(p->r) == 1)
		p->r = VecPar(p->r, nsp);
	if (VectorSize(p->K) == 1)
		p->K = VecPar(p->K, nsp);
	if (VectorSize(p->h2) == 1)
		p->h2 = VecPar(p->h2, nsp);
	if (VectorSize(p->V_s) == 1)
		p->V_s = VecPar(p->V_s, nsp);
	if (VectorSize(p->V_p) == 1)
		p->V_p = VecPar(p->V_p, nsp);
	if (VectorSize(p->V_u) == 1)
		p->V_u = VecPar(p->V_u, nsp);
	if (VectorSize(p->beta) == 1)
		p->beta = VecPar(p->beta, nsp);
	if (VectorSize(p->delta) == 1)
		p->delta = VecPar(p->delta, nsp);
	if (VectorSize(p->bbar) == 1)
		p->bbar = VecPar(p->bbar, nsp);
}

/* v has length 1, vn has length n */
Vector VecPar(Vector v, int n)
{
	double x = v[0];
	Vector vn;
	int i;

	deleteVector(v);
	vn = newVector(n);

	for (i=0; i<n; i++)
	{
		vn[i] = x;
	}

	return (vn);
}

void PrintParams(Params *p)
{
	printf("num_sp = %d\n", p->num_sp);
	printf("r = ");     printVector(p->r);
	printf("K = ");     printVector(p->K);
	printf("h2 = ");    printVector(p->h2);
	printf("V_s = ");   printVector(p->V_s);
	printf("V_p = ");   printVector(p->V_p);
	printf("V_u = ");   printVector(p->V_u);
	printf("beta = ");  printVector(p->beta);
	printf("delta = "); printVector(p->delta);
	printf("bbar = ");  printVector(p->bbar);
}
