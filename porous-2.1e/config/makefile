CC = mpicc
COMP = gnu
LDIR = /global/lib/checs

#CFLAGS = -O3 -no-prec-div
CFLAGS = -O3 -ffast-math
#CFLAGS = -Wall -g

INC = proto.h parms.h
OBJ = lognorm.o io.o

LGSL = -L$(LDIR)/gsl-$(COMP)/lib -lgsl -lgslcblas
LFFT = -L$(LDIR)/fftw-$(COMP)/lib -lfftw3
#LIBS = $(LGSL) $(LFFT) -limf -lm
LIBS = $(LGSL) $(LFFT) -lm

IGSL = -I$(LDIR)/gsl-$(COMP)/include
IFFT = -I$(LDIR)/fftw-$(COMP)/include
INCS = $(IGSL) $(IFFT)

lognorm: lognorm.o ./io.c
	$(CC) $(CFLAGS) -c ./io.c
	$(CC) $(CFLAGS) -o ../bin/lognorm lognorm.o io.o $(LIBS)

random: random.o ./io.c
	$(CC) $(CFLAGS) -c ./io.c
	$(CC) $(CFLAGS) -o ../bin/random random.o io.o $(LIBS)

.c.o:
	$(CC) $(CFLAGS) $(INCS) -c $<

$(OBJ): $(INC)

clean:	
	rm -rf *.o

all:
	make clean
	make lognorm
	make random
