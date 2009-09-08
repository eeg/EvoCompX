#ifndef __INPUT_H__
#define __INPUT_H__ 

#define TINY 1E-16

#include "keyvalue.h"

typedef struct
{
	/* biology */
	double r;
	double K;
	double h2;
	double V_s;
	double V_p;
	double V_u;
	double beta;
	double delta;
	const char *dispersal_file;

	/* landscape */
	int space_size;
	double opt_slope;
	const char *initial_num;
	const char *initial_zbar;

	/* record-keeping */
	double start_t;
	double stop_t;
	double record_interval;

} Params;

Params *NewParams();
void FreeParams(Params *params);
int AcquireParams(struct KeyValue *kv, Params *parameters);


#endif
