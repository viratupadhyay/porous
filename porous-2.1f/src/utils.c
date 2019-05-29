/* UTILS; Utility routines
        warning            (Prints warning messages on processor)
        error            (Prints error message - sets flag for exit)
        wtime            (Wall clock time)
        winit            (Initialize wall clock) */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <sys/time.h>
#include <time.h>
#include "proto.h"

#define MAX_W 10

void warning (char *msg)
{
    static int n_warn=0;

    fprintf (stdout, "warning: %s.\n", msg);
    n_warn++;
    if (n_warn >= MAX_W){
        sprintf(msg, "maximum number of warnings exceeded %d", MAX_W);
        error(msg);
    }
    fflush (stdout);
}

void error (char *msg)
{
        fprintf (stdout, "FATAL ERROR: %s\n", msg);
        fflush (stdout);
        exit(1);
}

void winit(struct timeval *t)
{
    gettimeofday(t, (struct timezone*)0);
}

double wtime(struct timeval *t1)
{
    struct timeval t2;
    double time;

    gettimeofday(&t2, (struct timezone*)0);
    time = ((t2.tv_sec-t1->tv_sec)*1000000
           + t2.tv_usec-t1->tv_usec)*1.0e-6;

    return time;
}
