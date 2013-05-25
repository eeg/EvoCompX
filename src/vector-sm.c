#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <glib.h>
#include "vector-sm.h"

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

