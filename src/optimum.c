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

#include "landscape.h"


/***
 * specify the optimum phenotype over space
 ***/


/* linear now, but could be anything */
double get_optimum(int cell_num, double slope)
{
	double opt;

	opt = slope*cell_num;

	return opt;
}

/*--------------------------------------------------
* // TODO: error checking for kink_stop < space_size
* double get_optimum(int cell_num, double slope)
* {
* 	double opt;
* 
* 	int kink_start = 5;
* 	int kink_stop = 10;
* 	double kink_steep = 6;    / * how much steeper than the main slope * /
* 
* 	if (cell_num < kink_start)
* 	{
* 		opt = cell_num * slope;
* 	}
* 	
* 	else if (cell_num > kink_stop)
* 	{
* 		opt = slope * kink_start + 
* 		      slope * kink_steep * (kink_stop - kink_start) + 
* 		      slope * (cell_num - kink_stop);
* 	}
* 	
* 	else
* 	{
* 		opt = slope * kink_start +
* 		      slope * kink_steep * (cell_num - kink_start);
* 	}
* 
* 	return opt;
* }
*--------------------------------------------------*/
