#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "landscape.h"

/***
 * Specify the optimum phenotype over space.
 ***/


/* linear now, but could be anything */
double get_optimum(int cell_num, double slope)
{
	double opt;

	opt = slope*cell_num;

	return opt;
}
