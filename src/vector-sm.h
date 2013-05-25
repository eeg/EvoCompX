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


/*******************************************************************************
 * These Vector structures are for use with keyvalue.
 * They were written by Walter Brisken.
 ******************************************************************************/


#ifndef __VECTOR__
#define __VECTOR__

typedef double*  Vector;
typedef double   VectorType;

#define VectorSize(v)   	(((int *)(v))[-1])

Vector newVector(int n);
Vector newVectorfromstring(const char *str);
int stringtoVector(const char *str, Vector v);
void deleteVector(Vector v);

#endif
