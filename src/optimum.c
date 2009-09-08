#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "landscape.h"

/***
 * Specify the optimum phenotype over space.
 ***/


double get_optimum(int cell_num, double slope)
{
	//double slope = 0.5;
	double opt;

	opt = slope*cell_num;

	return opt;
}
