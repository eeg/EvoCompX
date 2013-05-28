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


#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <glib.h>
#include "vector-sm.h"


/***
 * written by Walter Brisken
 * used in keyvalue.c
 ***/

static int totalvectordata = 0;

Vector newVector(int n)
{
	Vector v;

	v = g_new(VectorType, n+1);
	v++;
	VectorSize(v) = n;

	totalvectordata += n;

	return v;
}

Vector newVectorfromstring(const char *str)
{
	Vector v;
	int n;
	
	n = stringtoVector(str, 0);
	if(n == 0) return 0;
	
	v = newVector(n);
	stringtoVector(str, v);

	return v;
}

int stringtoVector(const char *str, Vector v)
{
	int n=0, l, i=0;
	double d;
	gchar *s;

	s = g_strdup(str);

	for(i = 0; s[i]; i++) 
	{
		if(s[i] >= '0' && s[i] <= '9') continue;
		if(s[i] == '.' || s[i] == '+' || s[i] == '-') continue;
		if((s[i] == 'e' || s[i] == 'E') &&
		   (s[i+1] == '+' || s[i+1] == '-')) continue;
		s[i] = ' ';
	}
	i = 0;

	while(sscanf(s+n, "%lf%n", &d, &l) > 0)
	{
		if(v) if(i < VectorSize(v)) v[i] = d;
		i++;
		n+=l;
	}

	g_free(s);

	return i;
}

void deleteVector(Vector v)
{
	totalvectordata -= VectorSize(v);
	g_free(v-1);
}

void printVector(const Vector v)
{
	int i, m;

	m = VectorSize(v);

	for(i = 0; i < m; i++)
	{
		if(i != 0) printf(" ");
		printf("%e", v[i]);
	}
	printf("\n");
}
