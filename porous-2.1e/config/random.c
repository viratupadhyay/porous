#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <fftw3.h>
#include <gsl/gsl_randist.h>
#include "mpi.h"
#include "proto.h"

int main(int argc, char ** argv)
{
    fftw_complex *z;
    fftw_plan plan;
    gsl_rng *r;
    const  gsl_rng_type *T;
    struct parms parms;
    double *phi, avg, ran, std;
    double *philocal;
    int Nx, Ny, Nz, size, rank, xnn;
    int i, j, n1, n2, irank, i0, j0;
    char filename[64];
    FILE *file;
    void   initrng(void);
    double uniform(void);
    MPI_Status status;

    Nx = atoi(argv[1]);
    Ny = atoi(argv[2]);
    Nz = atoi(argv[3]);
    avg = atof(argv[4]);
    std = atof(argv[5]);
    parms.Nx = Nx;  parms.Ny = Ny;  parms.Nz = Nz;
    //  MPI
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size (MPI_COMM_WORLD, &size);
    if (rank == 0){
        //  RNG setup
        gsl_rng_env_setup();
        T = gsl_rng_default;
        r = gsl_rng_alloc (T);
        printf ("r is a '%s' generator\n", gsl_rng_name (r));
        for (i=0; i < 1000000; i++)
            ran = avg + gsl_ran_gaussian(r,std); 
        xnn = Nx*Ny*Nz;
        phi = (double *) calloc(Nx*Ny*Nz, sizeof(double));
        for (i = 0; i < Nx*Ny*Nz; i++){
            ran = avg + gsl_ran_gaussian(r,std);
            phi[i] = (ran);
        }
    }


    parms.Nx = Nx;  parms.Ny = Ny;    parms.Nz = Nz;
    parms.rank = rank-1; parms.size = size-1;        //Slave COMM
    parms = mapping(parms);
    philocal = (double *) calloc(parms.nx*parms.ny, sizeof(double));
    
    if (rank == 0){
        for (irank = 0; irank < parms.size; irank++){
            parms.rank = irank;
            parms = mapping(parms);
            for (i = 0; i < parms.nx; i++)
            for (j = 0; j < parms.ny; j++){
                i0 = parms.x0 + i;        // Global coords
                j0 = parms.y0 + j;
                n1 = i*parms.ny + j;        //Local index
                n2 = i0*parms.Ny + j0;        //Global index
                philocal[n1] = phi[n2];
            }
            n1 = parms.nx*parms.ny;
            MPI_Send(philocal, n1, MPI_DOUBLE, irank+1, 0, MPI_COMM_WORLD);
        }
    }
    else{
        n1 = parms.nx*parms.ny;        //Fill local height
        MPI_Recv(philocal, n1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status);

        sprintf(filename, "%s.%02d", "phi.dat", parms.rank);
        file = fopen(filename, "w");
        fprintf(file, "%5d %5d %5d %5d %5d    0    0.0\n", parms.Nx, parms.Ny, parms.nx, parms.ny, parms.Nz);
        for (i = 0; i < parms.nx; i++)
        for (j = 0; j < parms.ny; j++){
            fprintf(file, " % .3e", philocal[i*parms.ny+j]);
            if ((i*parms.ny+j+1)%10 == 0) fprintf(file, "\n");
        }
    fclose(file);
    }
    MPI_Finalize();
    return 0;
}
