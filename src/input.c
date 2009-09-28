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

#include "input.h"
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
	free(params);
}


/*** read in the user's input parameter values ***/
int AcquireParams(struct KeyValue *kv, Params *parameters)
{
	parameters->num_sp = getKeyValueint(kv, "num_species");
     if (parameters->num_sp == KV_INTERR)
	{
		fprintf(stderr, "need to specify number of species, num_species\n");
		return -1;
	}
     if (parameters->num_sp < 2 || parameters->num_sp > MAX_NUM_SP)
	{
		fprintf(stderr, "need 2 <= num_species <= %d\n", MAX_NUM_SP);
		return -1;
	}

	/*** biology ***/

	parameters->r = getKeyValuedouble(kv, "r");
     if (parameters->r == KV_FLOATERR)
	{
		fprintf(stderr, "need to specify growth rate, r\n");
		return -1;
	}

	parameters->K = getKeyValuedouble(kv, "K");
     if (parameters->K == KV_FLOATERR)
	{
		fprintf(stderr, "need to specify carrying capacity, K\n");
		return -1;
	}

	parameters->h2 = getKeyValuedouble(kv, "h2");
     if (parameters->h2 == KV_FLOATERR)
	{
		fprintf(stderr, "need to specify heritability, h2\n");
		return -1;
	}

	parameters->V_s = getKeyValuedouble(kv, "V_s");
     if (parameters->V_s == KV_FLOATERR)
	{
		fprintf(stderr, "need to specify variance of stabilizing selection, "
		        "V_s\n");
		return -1;
	}

	parameters->V_p = getKeyValuedouble(kv, "V_p");
     if (parameters->V_p == KV_FLOATERR)
	{
		fprintf(stderr, "need to specify phenotypic variance, V_p\n");
		return -1;
	}

	parameters->V_u = getKeyValuedouble(kv, "V_u");
     if (parameters->V_u == KV_FLOATERR)
	{
		fprintf(stderr, "need to specify variance of competition function, "
		        "V_u\n");
		return -1;
	}

	parameters->alpha_file = getKeyValuestring(kv, "alpha_file");
     if (parameters->alpha_file == 0)
	{
		fprintf(stderr, "competition matrix not specified, using 1\n");
	}

	parameters->beta = getKeyValuedouble(kv, "beta");
     if (parameters->beta == KV_FLOATERR)
	{
		/* default is no hybridization */
		parameters->beta = 0;
	}

	parameters->delta = getKeyValuedouble(kv, "delta");
     if (parameters->delta == KV_FLOATERR)
	{
		fprintf(stderr, "need to specify dispersal probability, delta\n");
		return -1;
	}

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

	parameters->initial_num = getKeyValuestring(kv, "initial_num");
     if (parameters->initial_num == 0)
	{
		fprintf(stderr, "need to specify initial abundances, initial_num\n");
		return -1;
	}

	parameters->initial_zbar = getKeyValuestring(kv, "initial_zbar");
     if (parameters->initial_zbar == 0)
	{
		fprintf(stderr, "need to specify initial phenotypes, initial_zbar\n");
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
		exit(1);
	}

	/* overwrite parameter values with those specified on the command line */

     if (argc > 2) for(i = 2; i < argc; i++)
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
		exit(1);
	}

	return parameters;
}
