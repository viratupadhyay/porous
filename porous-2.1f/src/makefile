CC        = mpicc
HYPRE_DIR = /home/ladd/viratu/hypre-fnode/src/hypre

COPTS     = -g -Wall #-wd981
CINCLUDES = -I$(HYPRE_DIR)/include
CDEFS     = -DHAVE_CONFIG_H -DHYPRE_TIMING
CFLAGS    = $(COPTS) $(CINCLUDES) $(CDEFS)

INC = proto.h parms.h
OBJ = porous.o io.o functions.o hypresolve.o flux.o models.o utils.o

LINKOPTS  = $(COPTS)
LIBS      = -L$(HYPRE_DIR)/lib -lHYPRE -lm
LFLAGS    = $(LINKOPTS) $(LIBS) -lstdc++

.c.o:
	$(CC) $(CFLAGS) -c $<

porous: $(OBJ)
	$(CC) $(CFLAGS) -o ../bin/porous $(OBJ) $(LIBS)
	$(CC) -o $@ $^ $(LFLAGS)
	rm porous

clean:
	rm -f *.o

all:
	mkdir -p ../bin
	make clean
	make porous
