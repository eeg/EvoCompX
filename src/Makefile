SOURCE  = EvoCompX.c compete.c develop.c dispersal.c landscape.c optimum.c keyvalue.c input.c
HEADERS = compete.h develop.h dispersal.h landscape.h optimum.h keyvalue.h input.h
OBJECTS = ${SOURCE:.c=.o}
TARGETS = EvoCompX

.PREFIXES = .c .o 
CFLAGS = -Wall -g -O3
LDFLAGS = -lm
CC = gcc

%.o: %.c $(HEADERS)
	$(CC) -c -o $@ $< $(CFLAGS)

EvoCompX: ${OBJECTS}
	${CC} ${OBJECTS} ${LDFLAGS} -o $@
	
all: ${TARGETS}

clean:
	rm -f ${TARGETS} *.o
