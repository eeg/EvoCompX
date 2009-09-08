#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "input.h"
#include "keyvalue.h"

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
		fprintf(stderr, "need to specify variance of stabilizing selection,"
			    " V_s\n");
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
		fprintf(stderr, "need to specify variance of competition function,"
			    " V_u\n");
		return -1;
	}

	parameters->beta = getKeyValuedouble(kv, "beta");
     if (parameters->beta == KV_FLOATERR)
	{
		// default is no hybridization
		parameters->beta = 0;
	}

	parameters->delta = getKeyValuedouble(kv, "delta");
	parameters->dispersal_file = getKeyValuestring(kv, "dispersal_file");
     if (parameters->delta == KV_FLOATERR)
	{
		// if none specified, use dispersal file
     	if (parameters->dispersal_file == 0)
		{
			fprintf(stderr, "need to specify dispersal, either delta or " 
					"dispersal_file\n");
			return -1;
		}
		else
		{
			fprintf(stderr, "dispersal file %s specified, but can't yet get" 
					" values from it; use delta instead\n",
					parameters->dispersal_file);
			//return -1;
		}
	}
	else
	{
		if (parameters->dispersal_file != 0)
			fprintf(stderr, "ignoring dispersal_file because delta was " 
					"specified\n");
	}

	/*** landscape ***/

	parameters->space_size = getKeyValueint(kv, "space_size");
     if (parameters->space_size == KV_INTERR)
	{
		parameters->space_size = 100000;
	}

	parameters->opt_slope = getKeyValuedouble(kv, "opt_slope");
     if (parameters->opt_slope == KV_FLOATERR)
	{
		fprintf(stderr, "need to specify slope of optimum phenotype,"
			    " opt_slope\n");
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
		parameters->start_t = 0;
	}

	parameters->stop_t = getKeyValueint(kv, "stop_t");
     if (parameters->stop_t == KV_INTERR)
	{
		parameters->stop_t = 100000;
	}

	parameters->record_interval = getKeyValueint(kv, "record_interval");
     if (parameters->record_interval == KV_INTERR)
	{
		parameters->record_interval = 
			(parameters->stop_t - parameters->start_t)/10;
	}

	return 0;		// return -1 elsewhere if there's an error
}
