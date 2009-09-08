/***
  * local adaptation
  * keep track of numbers and mean phenotypes of each of two species
  * density dependent growth (equal intra- and inter-specific competition)
  *
  * no dispersal barriers
  *
  * one-dimensional space
  * discrete time
  * allow fancy dispersal kernels (may be different for the two species)
  *
  *** kernel file should look something like this:
	  note: odd number of lines, sums to zero
		0.1
		0.2
		-0.4	it's now okay to put any number here; normalization is automatic
		0.1
		0.0
	  indicates proportion of this cell added to neighboring cells (or subtracted from own)
***/


#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "dispersal.h"
#include "landscape.h"
#include "optimum.h"
#include "keyvalue.h"
#include "input.h"


int main(int argc, char *argv[])
{

/* might want to alter these */

	const int MAX_KERN_WIDTH = 3;
	char *kernel_file[2] = {"kernel.in", "kernel.in"};

	/*--------------------------------------------------
	* const int REAL_COLS = 100; 
	* const double start_t = 0;
	* const double stop_t = 100000;
	* double record_interval = (stop_t-start_t)/10;		// in units of t
	* / *** population growth parameters *** /
	* const double r = 0.1;
	* const double K = 10;
	* / *** natural selection parameters *** /
	* const double h2 = 0.6;
	* const double V_s = 200;
	* const double V_p = 1;
	* const double V_u = 40;
	* / *** parameter for frequency-independent hybridization (0=none) *** /
	* const double beta = 0.00;
	*--------------------------------------------------*/
	struct KeyValue *kv;
	char k[100], v[100];
	Params *parameters;

	double r;
	double K;
	double h2;
	double V_s; 
	double V_p;
	double V_u;
	double beta;
	double delta = 0.;    // if non-zero in param file, don't use kernel


/* no need to alter these */

	/*--------------------------------------------------
	* const int TOTAL_COLS = (REAL_COLS + MAX_KERN_WIDTH-1); 
	* const int REAL_START = ((MAX_KERN_WIDTH-1)/2);
	* const int REAL_STOP = (REAL_START + REAL_COLS -1);
	*--------------------------------------------------*/
	int TOTAL_COLS;
	int REAL_START;
	int REAL_STOP; 

	const int UNDEF = -9999;		// use this for undefined mean phenotypes

	// keep two copies of the space so events can happen simultaneously
	//Cell space[TOTAL_COLS][2];  FIXME: declared below, but should be up here
	int old, new;

	int t_steps;
	int recorded;

	double *base_kernel[2];
	int kernel_width[2], kernel_halfwidth[2];
	int kern_check = -99;

	FILE *fp_time, *fp_num1, *fp_num2, *fp_zbar1, *fp_zbar2;

	int i, j, t, sp;

	double n_old, n_other, n_new;
	double zbar_old, zbar_other, zbar_new;
	double w_bar, opt;

	/*** read in the parameter file ***/

	if (argc < 2)
	{
		fprintf(stderr, "need to specify a parameters file\n");
		return -1;
	}

	parameters = NewParams();
	kv = loadKeyValue(argv[1]);
	if (kv == 0)
	{
		fprintf(stderr, "unable to load parameters file\n");
		return -1;
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
		return -1;
	}

	r     = parameters->r;
	K     = parameters->K;
	h2    = parameters->h2;
	V_s   = parameters->V_s; 
	V_p   = parameters->V_p;
	V_u   = parameters->V_u;
	beta  = parameters->beta;
	delta = parameters->delta;

	TOTAL_COLS = (parameters->space_size + MAX_KERN_WIDTH-1); 
	REAL_START = ((MAX_KERN_WIDTH-1)/2);
	REAL_STOP = (REAL_START + parameters->space_size - 1);
	Cell space[TOTAL_COLS][2];  // FIXME: bad form to declare this down here

	t_steps = parameters->stop_t - parameters->start_t + 1;

	/*** if delta is specified, set up a simple kernel ***/
	if (delta > 0)
	{
		for (sp=0; sp<2; sp++)
		{
			kernel_halfwidth[sp] = 1;
			kernel_width[sp] = 3;

			base_kernel[sp] = allocate_kernel(kernel_width[sp]);
			base_kernel[sp][0] = delta;
			base_kernel[sp][1] = -2*delta;
			base_kernel[sp][2] = delta;
		}
	}

	/*** otherwise, read in the dispersal kernel to base_kernel ***/
	else
	{
		for (sp=0; sp<2; sp++)
		{
			kern_check = read_kernel(kernel_file[sp], &base_kernel[sp], 
					&kernel_width[sp], MAX_KERN_WIDTH);
			if (kern_check != 0)
				return -1;
			kernel_halfwidth[sp] = (kernel_width[sp]-1)/2;
		}
	}

	/*** assign dispersal kernels to each cell ***/
	for (sp=0; sp<2; sp++)
		assign_kernels(space, base_kernel[sp], kernel_width[sp], 
				MAX_KERN_WIDTH, TOTAL_COLS, REAL_START, REAL_STOP, sp);

	/*** clear the landscape and place initial individuals ***/
	initialize_landscape(space, parameters, TOTAL_COLS, REAL_START, REAL_STOP, UNDEF);

	old = 0;
	new = 1;

	/*** write the initial condition to the files ***/
	fp_time = fopen("time.dat", "w");
	fp_num1 = fopen("num1.dat", "w");
	fp_num2 = fopen("num2.dat", "w");
	fp_zbar1 = fopen("zbar1.dat", "w");
	fp_zbar2 = fopen("zbar2.dat", "w");
	recorded = 0;
	fprintf(fp_time, "%d\n", 0);
	record_landscape(fp_num1, fp_num2, fp_zbar1, fp_zbar2, space, old, 
			REAL_START, REAL_STOP);
	recorded++;


	printf("t = 0\n");

	for (t=1; t<t_steps; t++)
	{

		/*** competition and hybridization ***/
		
		for (i=REAL_START; i<=REAL_STOP; i++)
		{
			/* get the optimal phenotype for this cell */
			opt = get_optimum(i-REAL_START, parameters->opt_slope);

			/* effects of competition and stabilizing selection */

			for (sp=0; sp<2; sp++)
			{
				n_old = space[i][old].num[sp];

				if (n_old > 0)
				{
					zbar_old = space[i][old].zbar[sp];
					zbar_other = space[i][old].zbar[(sp+1)%2];
					n_other = space[i][old].num[(sp+1)%2];
					
					/* competition/logistic growth changes numbers */
					w_bar = r - (r/K)*sqrt(V_u/(V_p+V_u))*(n_old + n_other*exp((-1/(4*(V_p+V_u)))*pow(zbar_old-zbar_other, 2))) - V_p/(2*V_s) - pow(opt-zbar_old, 2)/(2*V_s);
					n_new = exp(w_bar) * n_old;

					/* competition/selection changes mean phenotype */
					zbar_new = zbar_old + h2*V_p*((opt-zbar_old)/V_s + (r/K)*sqrt(V_u/(V_p+V_u))*n_other*(zbar_old-zbar_other)/(2*(V_p+V_u))*exp((-1/(4*(V_p+V_u)))*pow(zbar_old-zbar_other, 2)));

					/* hybridization changes numbers */
					n_new = n_new * n_old/(n_old + beta*n_other);
				}
				else 
				{
					n_new = n_old;
					zbar_new = space[i][old].zbar[sp];
				}

				/* put the post-competition values into "new" (not "old" because then sp1's changes would affect sp2 immediately) */
				space[i][new].num[sp] = n_new;
				space[i][new].zbar[sp] = zbar_new;
				space[i][new].ztotal[sp] = zbar_new * n_new;
			} 
		}

		/* swap new and old if they didn't get swapped at the end of the spatial loop */
		new = old;
		old = (new+1)%2;


		/*** dispersal ***/
		empty_border(space, REAL_START, REAL_STOP, TOTAL_COLS, UNDEF); // empties old and new
		dispersal_happens(space, kernel_halfwidth, old, REAL_START, REAL_STOP, UNDEF);


		/*** record the new state, if it's time to ***/
		if (t == recorded * parameters->record_interval)
		{
			fprintf(fp_time, "%d\n", t);
			record_landscape(fp_num1, fp_num2, fp_zbar1, fp_zbar2, space, new, REAL_START, REAL_STOP);
			recorded++;
			printf("t = %d\n", t);
//printf("n1 = %f, n2 = %f, z1bar = %f, z2bar = %f\n", space[20+REAL_START][old].num[0], space[20+REAL_START][old].num[1], space[20+REAL_START][old].zbar[0], space[20+REAL_START][old].zbar[1]);
		}

		/*** the new state becomes the old state ***/
		new = old;
		old = (new+1)%2;
	}


	fclose(fp_time);
	fclose(fp_num1);
	fclose(fp_num2);
	fclose(fp_zbar1);
	fclose(fp_zbar2);


	/*** record the final state ***/
	fp_time = fopen("final.dat", "w");
	for (i=REAL_START; i<=REAL_STOP; i++)
		fprintf(fp_time, "%e\t%e\t%e\t%e\n", space[i][old].num[0], space[i][old].num[1], space[i][old].zbar[0], space[i][old].zbar[1]);
	fclose(fp_time);

	FreeParams(parameters);


	return 0;
}
