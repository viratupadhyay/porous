#include <stdlib.h>
#include <stdio.h>
#include "proto.h"
#include "parms.h"

extern struct parms avgp (double *c, double *r, double *q, double *h, char *fname, struct parms parms)
{
    double havg, delh;
    FILE *fileptr;

    havg = avgf(data, 0, parms);
    delh = havg - parms.havg;
