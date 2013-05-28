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

void FreeParams(Params *params)
{
	if (params->r != NULL)
		deleteVector(params->r);
	free(params);
}


/*** read in the user's input parameter values ***/
int AcquireParams(struct KeyValue *kv, Params *parameters)
{
	int nsp;

	parameters->num_sp = getKeyValueint(kv, "num_species");
     if (parameters->num_sp == KV_INTERR)
	{
		fprintf(stderr, "need to specify number of species, num_species\n");
		return ERROR;
	}
     if (parameters->num_sp < 2 || parameters->num_sp > MAX_NUM_SP)
	{
		fprintf(stderr, "need 2 <= num_species <= %d\n", MAX_NUM_SP);
		return ERROR;
	}
	nsp = parameters->num_sp;

	/*** biology ***/

	parameters->r = getKeyValueVector(kv, "r");
	if (CheckParam(parameters->r, nsp, "growth rate, r") == ERROR)
		return ERROR;

	parameters->K = getKeyValueVector(kv, "K");
	if (CheckParam(parameters->K, nsp, "carrying capacity, K") == ERROR)
		return ERROR;

	parameters->h2 = getKeyValueVector(kv, "h2");
	if (CheckParam(parameters->h2, nsp, "heritability, h2") == ERROR)
		return ERROR;

	parameters->V_p = getKeyValuedouble(kv, "V_p");
	parameters->V_u = getKeyValuedouble(kv, "V_u");

	parameters->V_s = getKeyValueVector(kv, "V_s");
	if (CheckParam(parameters->V_s, nsp, "variance of stabilizing selection, V_s") == ERROR)
		return ERROR;

/*
	parameters->V_p = getKeyValueVector(kv, "V_p");
	if (CheckParam(parameters->V_p, nsp, "phenotypic variance, V_p") == ERROR)
		return ERROR;

	parameters->V_u = getKeyValueVector(kv, "V_u");
	if (CheckParam(parameters->V_u, nsp, "variance of competition function, V_u") == ERROR)
		return ERROR;
*/

	parameters->alpha_file = getKeyValuestring(kv, "alpha_file");
     if (parameters->alpha_file == 0)
	{
		fprintf(stderr, "competition matrix not specified, using 1\n");
	}

	parameters->beta = getKeyValueVector(kv, "beta");
     if (parameters->beta == 0)
	{
		/* default is no hybridization */
		parameters->beta = newVector(1);
		// default is no hybridization
		parameters->beta[0] = 0;
	}
	else
	{
		if (CheckParam(parameters->beta, nsp, "hybridization consideration, beta") == ERROR)
			return ERROR;
	}

	parameters->delta = getKeyValueVector(kv, "delta");
	if (CheckParam(parameters->delta, nsp, "dispersal probability, delta") == ERROR)
		return ERROR;

	parameters->bbar = getKeyValueVector(kv, "bbar");
	if (CheckParam(parameters->bbar, nsp, "plasticity, bbar") == ERROR)
		return ERROR;

	/*** landscape ***/

	parameters->space_size = getKeyValueint(kv, "space_size");
     if (parameters->space_size == KV_INTERR || parameters->space_size < 1)
	{
		fprintf(stderr, "valid space_size not specified, using 100\n");
		parameters->space_size = 100;
	}

	parameters->opt_slope = getKeyValuedouble(kv, "opt_slope");
     if (parameters->opt_slope == KV_FLOATERR)
	{
		fprintf(stderr, "need to specify slope of optimum phenotype, "
			    "opt_slope\n");
		return -1;
	}

	/* TODO need different slopes for theta and epsilon, i.e., B and C */
	parameters->env_slope = 1;

	parameters->initial_num = getKeyValuestring(kv, "initial_num");
     if (parameters->initial_num == 0)
	{
		fprintf(stderr, "need to specify initial abundances, initial_num\n");
		return -1;
	}

	parameters->initial_abar = getKeyValuestring(kv, "initial_abar");
     if (parameters->initial_abar == 0)
	{
		fprintf(stderr, "need to specify initial breeding values, initial_abar\n");
		return -1;
	}

	/*** record-keeping ***/

	parameters->start_t = getKeyValueint(kv, "start_t");
     if (parameters->start_t == KV_INTERR)
	{
		fprintf(stderr, "valid start_t not specified, using 0\n");
		parameters->start_t = 0;
	}

	parameters->stop_t = getKeyValueint(kv, "stop_t");
     if (parameters->stop_t == KV_INTERR || 
	    parameters->stop_t < parameters->start_t)
	{
		fprintf(stderr, "valid stop_t not specified, using start_t + 1000\n");
		parameters->stop_t = parameters->start_t + 1000;
	}

	parameters->record_interval = getKeyValueint(kv, "record_interval");
     if (parameters->record_interval == KV_INTERR || 
	    parameters->record_interval <= 0)
	{
		fprintf(stderr, "valid record_interval not specified, using elapsed "
		                "time/10\n");
		parameters->record_interval = 
			(parameters->stop_t - parameters->start_t)/10;
	}

	return 0;		/* returns -1 elsewhere if there's an error */
}


Params *GetParams(int argc, char *argv[])
{
	struct KeyValue *kv;
	char k[100], v[100];
	Params *parameters;
	int i, j;

	/* open the user's file */

	if (argc < 2)
	{
		fprintf(stderr, "\n           *** This is EvoCompX ver 0.2 ***\n"
		     "  Input parameters are required.  Please see the README and\n"
		     "  example params.in file for usage examples and options.\n\n");
		exit(1);
	}

	parameters = NewParams();
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

	if (AcquireParams(kv, parameters) == -1)
	{
		fprintf(stderr, "unable to proceed -- "
				"not all required parameters specified\n");
		exit(ERROR);
	}

	/* convert any scalar parameters into a vector of length num_species */
	VectorizeParams(parameters);

	/*
	PrintParams(parameters);
	exit(ERROR); // XXX
	*/

	return parameters;
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

/*
	if (VectorSize(p->V_p) == 1)
		p->V_p = VecPar(p->V_p, nsp);

	if (VectorSize(p->V_u) == 1)
		p->V_u = VecPar(p->V_u, nsp);
*/

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
/*
	printf("V_p = ");   printVector(p->V_p);
	printf("V_u = ");   printVector(p->V_u);
*/
	printf("beta = ");  printVector(p->beta);
	printf("delta = "); printVector(p->delta);
	printf("bbar = ");  printVector(p->bbar);
}
/*
	const char *alpha_file;
	double alpha[MAX_NUM_SP][MAX_NUM_SP];

	int space_size;
	double opt_slope;
	double env_slope;
	const char *initial_num;
	const char *initial_abar;

	int start_t;
	int stop_t;
	int record_interval;
*/
