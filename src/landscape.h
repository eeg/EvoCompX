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


#ifndef __LANDSCAPE_H__
#define __LANDSCAPE_H__

#include "input.h"

#define MAX_NUM_SP 20
#define MAX_SPACE_SIZE 1000

/*** spatial units on the landscape ***/
typedef struct cell
{
	double num[MAX_NUM_SP];      /* number of individuals of each species */
	double zbar[MAX_NUM_SP];     /* mean phenotype for each species       */
	double ztotal[MAX_NUM_SP];   /* num * zbar, for each species          */
} Cell;

void record_landscape(FILE **fp_num, FILE **fp_zbar, Cell space[][2],
                      Params *params, int old_new);

void initialize_landscape(Cell space[][2], Params *params);

#endif
