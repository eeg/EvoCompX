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
#include <math.h>
#include <stdlib.h>

#include "compete.h"
#include "optimum.h"


/*** 
 * competition, stabilizing selection, hybridization 
 ***/


/* construct the matrix of competition coefficients */
void make_alpha(Params *params)
{
	FILE *fp;
	int i, j;

	/* if the alpha's were provided */
	if (params->alpha_file != 0)
	{
		fp = fopen(params->alpha_file, "r");
		for (i=0; i<params->num_sp; i++)
		{
			for (j=0; j<params->num_sp; j++)
			{
				if (fscanf(fp, "%lf", &params->alpha[i][j]) != 1)
				{
					fprintf(stderr, "Error: invalid input in "
					                "alpha_file\n");
					exit(1);
				}
			}
		}
		fclose(fp);
	}

	/* if the alpha's weren't provided, set them all to 1, as in CT2000 */
	else
	{
		for (i=0; i<params->num_sp; i++)
			for (j=0; j<params->num_sp; j++)
				params->alpha[i][j] = 1;
	}
}


void competition_happens(Cell space[][2], int old, Params *params)
{
	double nza_new[3];           /* n_new, zbar_new, abar_new */
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
				/* competition and selection update nza_new[] */
				comp_sel(nza_new, sp, i, opt, space, old, params);

				/* hybridization changes numbers */
				if (params->beta[sp] > 0)
				{
					n_other = 0;
					for (osp=0; osp<params->num_sp; osp++)
					{
						if (osp != sp)
						{
							n_other += space[i][old].num[osp];
						}
					}
					nza_new[0] = nza_new[0] * 
							  n_old/(n_old + params->beta[sp]*n_other);
				}
			}
			else            /* if the species is absent, nothing changes */
			{
				nza_new[0] = n_old;
				nza_new[1] = space[i][old].zbar[sp];
				nza_new[2] = space[i][old].abar[sp];
			}

			/* put the post-competition values into "new" (not "old" 
			 * because then sp1's changes would affect sp2 immediately) */
			space[i][new].num[sp] = nza_new[0];
			space[i][new].abar[sp] = nza_new[2];
			space[i][new].atotal[sp] = nza_new[2] * nza_new[0];
			space[i][new].zbar[sp] = nza_new[1];
			/* note that mean phenotype is empty now; it will be
			 * re-determined during development
			 * space[i][new].ztotal[sp] = nza_new[1] * nza_new[0]; */
		} 
	}
}

/*** competition and selection change num and zbar ***/
void comp_sel(double nza_new[3], int sp, int i, double opt, Cell space[][2], 
              int old, Params *p)
{
	double n_old, n_other;
	double zbar_old, zbar_other;
	double var_temp, w_temp, z_temp, wz_temp, w_bar;
	int osp;

	n_old = space[i][old].num[sp];
	zbar_old = space[i][old].zbar[sp];

	w_temp = 0;
	z_temp = 0;

	/* sum up the effects from each "other" species */
	for (osp=0; osp<p->num_sp; osp++)
	{
		if (osp != sp)
		{
			zbar_other = space[i][old].zbar[osp];
			n_other = space[i][old].num[osp];

			var_temp = (p->V_p[sp] + p->V_p[osp] + p->V_u[sp] + p->V_u[osp]);

			wz_temp = p->alpha[sp][osp] * n_other * 
					sqrt((p->V_u[sp] + p->V_u[osp]) / var_temp) * 
					exp(-pow(zbar_old - zbar_other, 2) / (2 * var_temp));

			w_temp += wz_temp;
			z_temp += wz_temp * (zbar_old - zbar_other) / var_temp;
		}
	}

	/* mean fitness: pop growth, stabilizing selection, inter- and intraspecific competition */
	w_bar = p->r[sp] - (p->V_p[sp] + pow(opt - zbar_old, 2)) / ( 2 *p->V_s[sp])
		-(p->r[sp] / p->K[sp]) * 
			( w_temp +
				p->alpha[sp][sp] * n_old * sqrt(p->V_u[sp] / (p->V_u[sp] + p->V_p[sp]))
			);

	/* population size changes */
	nza_new[0] = exp(w_bar) * n_old;

	/* mean breeding value changes */
	nza_new[2] = space[i][old].abar[sp] + p->h2[sp] * p->V_p[sp] * (
		(opt - zbar_old) / p->V_s[sp] + p->r[sp] / p->K[sp] * z_temp );

	/* mean phenotype will be set during development, which should occur
	 * immediately after competition/reproduction */
	nza_new[1] = UNDEF_PHEN;
}
