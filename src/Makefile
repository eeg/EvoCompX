SOURCE  = EvoCompX.c compete.c converge.c develop.c disperse.c landscape.c optimum.c keyvalue.c input.c vector-sm.c
HEADERS = compete.h converge.h develop.h disperse.h landscape.h optimum.h keyvalue.h input.h vector-sm.h
OBJECTS = ${SOURCE:.c=.o}
TARGETS = EvoCompX

.PREFIXES = .c .o 
CFLAGS = -Wall -g -O3 -mpc64 `pkg-config --cflags glib-2.0`
LDFLAGS = `pkg-config --libs glib-2.0` -lm
CC = gcc

%.o: %.c $(HEADERS)
	$(CC) -c -o $@ $< $(CFLAGS)

EvoCompX: ${OBJECTS}
	${CC} -o $@ ${OBJECTS} ${LDFLAGS}
	
all: ${TARGETS}

clean:
	rm -f ${TARGETS} *.o
