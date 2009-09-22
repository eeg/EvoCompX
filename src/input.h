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
	double start_t;
	double stop_t;
	double record_interval;

} Params;

Params *NewParams();
void FreeParams(Params *params);
int AcquireParams(struct KeyValue *kv, Params *parameters);
Params *GetParams(int argc, char *argv[]);

#endif
