SOURCE  = compete.c dispersal.c landscape.c optimum.c keyvalue.c input.c
HEADERS = dispersal.h landscape.h optimum.h keyvalue.h input.h
OBJECTS = ${SOURCE:.c=.o}
TARGETS = compete

.PREFIXES = .c .o 
CFLAGS = -Wall -g
LDFLAGS = -lm
CC = gcc

.c.o: 
	${CC} -c ${CFLAGS} $< 

compete: ${OBJECTS}
	${CC} ${LDFLAGS} ${OBJECTS} -o $@
	
all: ${TARGETS}

clean:
	rm ${TARGETS} *.o
