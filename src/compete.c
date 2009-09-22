/***
  * keep track of numbers and mean phenotypes of each of several species
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
	/* FIXME: don't use these? */
	double r;
	double K;
	double h2;
	double V_s; 
	double V_p;
	double V_u;
	double beta;
	double delta;
	int space_size;
	int num_sp;

	/* keep two copies of the space so events can happen simultaneously */
	Cell space[MAX_SPACE_SIZE][2];
	int old, new;

	/* output files */
	FILE *fp_time, *fp_num[MAX_NUM_SP], *fp_zbar[MAX_NUM_SP];
	char str[1000];

	int i, t, sp, osp;

	int t_steps;
	int recorded;

	double n_old, n_other, n_new;
	double zbar_old, zbar_other, zbar_new;
	double w_bar, opt, w_temp, z_temp;


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

	num_sp = parameters->num_sp;
	space_size = parameters->space_size;
	t_steps = parameters->stop_t - parameters->start_t + 1;


	/*** clear the landscape and place initial individuals ***/

	initialize_landscape(space, parameters);

	old = 0;
	new = 1;

	/*** write the initial condition to the output files ***/

	fp_time = fopen("time.dat", "w");
	for (sp=0; sp<num_sp; sp++)
	{
		sprintf(str, "num%d.dat", sp+1);
		fp_num[sp] = fopen(str, "w");

		sprintf(str, "zbar%d.dat", sp+1);
		fp_zbar[sp] = fopen(str, "w");
	}

	recorded = 0;
	fprintf(fp_time, "%d\n", 0);
	record_landscape(fp_num, fp_zbar, space, parameters, old);
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

			for (sp=0; sp<num_sp; sp++)
			{
				n_old = space[i][old].num[sp];
				zbar_old = space[i][old].zbar[sp];

				/* FIXME functionalize ?  */
				if (n_old > 0)
				{

					/* competition/selection change num and zbar */
					w_temp = 0;
					z_temp = 0;
					n_other = 0;
					for (osp=0; osp<num_sp; osp++) /* "other" species */
					{
						if (osp != sp)
						{
							zbar_other = space[i][old].zbar[osp];
							n_other = space[i][old].num[osp];

							w_temp += n_other * exp((-1/(4*(V_p+V_u)))*pow(zbar_old-zbar_other, 2));
							z_temp += n_other * (zbar_old - zbar_other) * exp((-1/(4*(V_p+V_u)))*pow(zbar_old-zbar_other, 2));
						}
					}

					w_bar = r - pow(opt-zbar_old, 2)/(2*V_s) - V_p/(2*V_s) - (r/K)*sqrt(V_u/(V_p+V_u)) * (n_old + w_temp);

					n_new = exp(w_bar) * n_old;

					zbar_new = zbar_old + h2*V_p*((opt-zbar_old)/V_s + r/(2*K)*sqrt(V_u)/pow(V_p+V_u, 1.5) * z_temp);
							
					/* hybridization changes numbers */
					n_new = n_new * n_old/(n_old + beta*n_other);
				}
				else 
				{
					n_new = n_old;
					zbar_new = zbar_old;
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
		dispersal_happens(space, old, parameters);


		/*** record the new state, if it's time to ***/
		if (t == recorded * parameters->record_interval)
		{
			fprintf(fp_time, "%d\n", t);
			record_landscape(fp_num, fp_zbar, space, parameters, new);
			recorded++;
			printf("t = %d\n", t);
		}

		/* swap new and old, in preparation for the next generation */
		new = old;
		old = (new+1)%2;
	}

	/* got all the results now */
	fclose(fp_time);
	for (sp=0; sp<num_sp; sp++)
	{
		fclose(fp_num[sp]);
		fclose(fp_zbar[sp]);
	}
	FreeParams(parameters);

	/*** record the final state ***/
	fp_time = fopen("final.dat", "w");
	for (i=0; i<space_size; i++)
	{
		for (sp=0; sp<num_sp; sp++)
			fprintf(fp_time, "%e\t", space[i][old].num[sp]); 
		for (sp=0; sp<num_sp; sp++)
			fprintf(fp_time, "%e\t", space[i][old].zbar[sp]); 
		fprintf(fp_time, "\n");
	}
	fclose(fp_time);


	return 0;
}
