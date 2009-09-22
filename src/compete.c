#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "compete.h"
#include "optimum.h"


/*** 
 * competition, stabilizing selection, hybridization 
 ***/


void competition_happens(Cell space[][2], int old, Params *params)
{
	double nz_new[2];           /* n_new, zbar_new */
	double n_old, n_other;
	int new = (old+1)%2;
	double opt;
	int i, sp, osp;

	for (i=0; i<params->space_size; i++)
	{
		/* get the optimal phenotype for this cell */
		opt = get_optimum(i, params->opt_slope);

		/*** effects of competition, stabilizing selection, hybridization ***/

		for (sp=0; sp<params->num_sp; sp++)
		{
			n_old = space[i][old].num[sp];

			if (n_old > 0)  /* if the species is present, do stuff */
			{
				/* competition and selection update nz_new[] */
				comp_sel(nz_new, sp, i, opt, space, old, params);

				/* hybridization changes numbers */
				if (params->beta > 0)
				{
					n_other = 0;
					for (osp=0; osp<params->num_sp; osp++)
						if (osp != sp)
							n_other += space[i][old].num[osp];
					nz_new[0] = nz_new[0] * 
							  n_old/(n_old + params->beta*n_other);
				}
			}
			else            /* if the species is absent, nothing changes */
			{
				nz_new[0] = n_old;
				nz_new[1] = space[i][old].zbar[sp];
			}

			/* put the post-competition values into "new" (not "old"
 			 *   because then sp1's changes would affect sp2 immediately) */
			space[i][new].num[sp] = nz_new[0];
			space[i][new].zbar[sp] = nz_new[1];
			space[i][new].ztotal[sp] = nz_new[1] * nz_new[0];
		} 
	}
}

/*** competition and selection change num and zbar ***/
void comp_sel(double nz_new[2], int sp, int i, double opt, Cell space[][2], 
              int old, Params *p)
{
	double n_old, n_other;
	double zbar_old, zbar_other;
	double w_temp, z_temp, w_bar;
	int osp;

	w_temp = 0;
	z_temp = 0;
	n_other = 0;

	n_old = space[i][old].num[sp];
	zbar_old = space[i][old].zbar[sp];

	for (osp=0; osp<p->num_sp; osp++) /* "other" species */
	{
		if (osp != sp)
		{
			zbar_other = space[i][old].zbar[osp];
			n_other = space[i][old].num[osp];

			w_temp += n_other * exp((-1/(4*(p->V_p+p->V_u))) *
			                         pow(zbar_old-zbar_other, 2));
			z_temp += n_other * (zbar_old - zbar_other) * 
			          exp((-1/(4*(p->V_p+p->V_u))) *
			              pow(zbar_old-zbar_other, 2));
		}
	}

	w_bar = p->r - pow(opt-zbar_old, 2)/(2*p->V_s) - p->V_p/(2*p->V_s) - 
	        (p->r/p->K) * sqrt(p->V_u/(p->V_p+p->V_u)) * (n_old + w_temp);

	nz_new[0] = exp(w_bar) * n_old;

	nz_new[1] = zbar_old + p->h2*p->V_p*((opt-zbar_old)/p->V_s + 
	            p->r/(2*p->K) * sqrt(p->V_u)/pow(p->V_p+p->V_u, 1.5) * z_temp);
}
