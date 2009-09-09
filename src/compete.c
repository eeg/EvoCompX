/***
  * keep track of numbers and mean phenotypes of each of two species
  * local adaptation
  * density dependent growth (equal intra- and inter-specific competition)
  * one-dimensional space
  * discrete time
  * nearest-neighbor dispersal (though could be repeated within a generation)
  * no dispersal barriers
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
	Params *parameters;

	/* all of these are specified in the input parameter file */
	double r;
	double K;
	double h2;
	double V_s; 
	double V_p;
	double V_u;
	double beta;
	double delta;
	int space_size;

	/* use this for undefined mean phenotypes */
	const int UNDEF = -9999;

	/* keep two copies of the space so events can happen simultaneously */
	Cell space[10000][2];  /* MAX_SPACE_SIZE */
	int old, new;

	/* output files */
	FILE *fp_time, *fp_num1, *fp_num2, *fp_zbar1, *fp_zbar2;

	int i, t, sp;

	int t_steps;
	int recorded;

	double n_old, n_other, n_new;
	double zbar_old, zbar_other, zbar_new;
	double w_bar, opt;


	/*** read in the parameter file ***/

	parameters = GetParams(argc, argv);

	r     = parameters->r;
	K     = parameters->K;
	h2    = parameters->h2;
	V_s   = parameters->V_s; 
	V_p   = parameters->V_p;
	V_u   = parameters->V_u;
	beta  = parameters->beta;
	delta = parameters->delta;

	space_size = parameters->space_size;
	t_steps = parameters->stop_t - parameters->start_t + 1;


	/*** clear the landscape and place initial individuals ***/

	initialize_landscape(space, parameters, space_size, UNDEF);

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
	                 space_size);
	recorded++;


	printf("t = 0\n");

	for (t=1; t<t_steps; t++)
	{
		/*** competition and hybridization ***/
		
		for (i=0; i<space_size; i++)
		{
			/* get the optimal phenotype for this cell */
			opt = get_optimum(i, parameters->opt_slope);

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

				/* put the post-competition values into "new" (not "old"
 				 *   because then sp1's changes would affect sp2
				 *   immediately) */
				space[i][new].num[sp] = n_new;
				space[i][new].zbar[sp] = zbar_new;
				space[i][new].ztotal[sp] = zbar_new * n_new;
			} 
		}

		/* swap new and old, since dispersal happens after competition */
		new = old;
		old = (new+1)%2;


		/*** dispersal ***/
		dispersal_happens(space, old, space_size, delta, UNDEF);


		/*** record the new state, if it's time to ***/
		if (t == recorded * parameters->record_interval)
		{
			fprintf(fp_time, "%d\n", t);
			record_landscape(fp_num1, fp_num2, fp_zbar1, fp_zbar2, space,
			                 new, space_size);
			recorded++;
			printf("t = %d\n", t);
		}

		/* swap new and old, in preparation for the next generation */
		new = old;
		old = (new+1)%2;
	}

	/* got all the results now */
	fclose(fp_time);
	fclose(fp_num1);
	fclose(fp_num2);
	fclose(fp_zbar1);
	fclose(fp_zbar2);
	FreeParams(parameters);

	/*** record the final state ***/
	fp_time = fopen("final.dat", "w");
	for (i=0; i<space_size; i++)
		fprintf(fp_time, "%e\t%e\t%e\t%e\n", space[i][old].num[0],
		        space[i][old].num[1], space[i][old].zbar[0], 
		        space[i][old].zbar[1]);
	fclose(fp_time);


	return 0;
}
