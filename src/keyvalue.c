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
#include <string.h>
#include <stdlib.h>
#include <glib.h>

#include "keyvalue.h"


/***
 * written by Walter Brisken
 * used in input.c
 ***/


struct KeyValue *newKeyValue()
{
	struct KeyValue *p;

	p = g_new(struct KeyValue, 1);
	p->nparms = 0;

	return p;
}

void deleteKeyValue(struct KeyValue *p)
{
	int i;
	
	if(!p) return;
	
	for(i = 0; i < p->nparms; i++)
	{
		g_free(p->key[i]);
		g_free(p->value[i]);
	}
	g_free(p);
}

void KeyValueaddparm(struct KeyValue *p, const char *key, const char *value)
{
	p->key[p->nparms]   = g_strdup(key);
	p->value[p->nparms] = g_strdup(value);
	p->nparms++;
}

void KeyValueupdateparm(struct KeyValue *p, const char *key, const char *value)
{
	int i;

	i = KeyValuekeyindex(p, key);
	if(i >= 0)
	{
		g_free(p->value[i]);
		p->value[i] = g_strdup(value);
	}
	else
	{
		p->key[p->nparms]   = g_strdup(key);
		p->value[p->nparms] = g_strdup(value);
		p->nparms++;
	}
}

void KeyValueupdateparmdouble(struct KeyValue *p, const char *key, double value)
{
	int i;
	char v[100];

	sprintf(v, "%e", value);

	i = KeyValuekeyindex(p, key);
	if(i >= 0)
	{
		g_free(p->value[i]);
		p->value[i] = g_strdup(v);
	}
	else
	{
		p->key[p->nparms]   = g_strdup(key);
		p->value[p->nparms] = g_strdup(v);
		p->nparms++;
	}
}

struct KeyValue *loadKeyValue(const char *filename)
{
	struct KeyValue *p;
	FILE *in;
	int i;
	char str[1000], K[500], V[500];

	if(strcmp(filename, "-") == 0) in = stdin;
	else in = fopen(filename, "r");

	if(!in)
	{
		fprintf(stderr, "loadKeyValue: no file: %s\n", filename);
		return 0;
	}

	p = newKeyValue();

	while(1)
	{
		int q;

		if(fgets(str, 999, in) == 0);
		if(feof(in)) break;
		for(i = 0; str[i] != 0; i++)
		{
			if(str[i] == '#') str[i] = 0;
			else if(str[i] == '=') str[i] = ' ';
		}
		for(i--; str[i] <= ' '; i--);
		
		/* remove trailing whitespace */
		str[i+1] = 0;
		if(str[0] == 0) continue;

		/* allow space separating vector elements */
		if(sscanf(str, "%s %n%s", K, &q, V) != 2) continue;
		while(str[q] == ' ')
		{
			++q;
		}

		KeyValueupdateparm(p, K, str+q);
	}

	if(strcmp(filename, "-") != 0) fclose(in);
	
	return p;
}

int saveKeyValue(const struct KeyValue *kv, const char *filename)
{
	int i;
	FILE *out;

	out = fopen(filename, "w");
	if(!out) return 0;

	for(i = 0; i < kv->nparms; i++)
		fprintf(out, "%s = %s\n", kv->key[i], kv->value[i]);
	
	fclose(out);

	return 1;
}

/* return -1 if not in the list of keys */
int KeyValuekeyindex(const struct KeyValue *p, const char *key)
{
	int i;

	for(i = 0; i < p->nparms; i++)
		if(strcmp(p->key[i], key) == 0) return i;
	
	return -1;
}

long long int getKeyValuelonglongint(const struct KeyValue *p, const char *key)
{
	int i = KeyValuekeyindex(p, key);

	if(i < 0) return KV_INTERR;
	return atoll(p->value[i]);
}

int getKeyValueint(const struct KeyValue *p, const char *key)
{
	int i = KeyValuekeyindex(p, key);

	if(i < 0) return KV_INTERR;
	return atoi(p->value[i]);
}

double getKeyValuedouble(const struct KeyValue *p, const char *key)
{
	int i = KeyValuekeyindex(p, key);

	if(i < 0) return KV_FLOATERR;
	return atof(p->value[i]);
}

Vector getKeyValueVector(const struct KeyValue *p, const char *key)
{
	int i = KeyValuekeyindex(p, key);

	if(i < 0) return 0;
	return newVectorfromstring(p->value[i]);
}
 
const char *getKeyValuestring(const struct KeyValue *p, const char *key)
{
	int i = KeyValuekeyindex(p, key);

	if(i < 0) return 0;
	return p->value[i];
}

void printKeyValue(const struct KeyValue *p)
{
	int i;

	for(i = 0; i < p->nparms; i++)
		printf("  %s = %s\n", p->key[i], p->value[i]);
}
