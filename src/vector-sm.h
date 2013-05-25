#ifndef __VECTOR__
#define __VECTOR__

typedef double*  Vector;
typedef double   VectorType;

#define VectorSize(v)   	(((int *)(v))[-1])


Vector newVector(int n);
Vector newVectorfromstring(const char *str);
int stringtoVector(const char *str, Vector v);


#endif
