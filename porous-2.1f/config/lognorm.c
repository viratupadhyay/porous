/* Lognorm - log-normal aperture distribution - Master/Slave */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "mpi.h"
#include <fftw3.h>
#include <gsl/gsl_rng.h>
#include "proto.h"

int main(int argc, char **argv)
{
    fftw_complex *z;
    fftw_plan plan;
    gsl_rng *r;
    const  gsl_rng_type *T;
    void   initrng(void);
    struct parms parms;
    double uniform(void);
    double *h, avg, rgh, var;
    double *hlocal;
    double kx, ky, k2, lambda, lambdax, lambday;
    double norm, ph, radial, xnn, ran, Pi;
    int    Nx, Ny, Nz, Nx2, Ny2;
    int    i, j, i0, j0, n1, n2;
    int    rank, size, irank;
    char   filename[64];
    FILE   *file;
    MPI_Status status;
    
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    Nx = atoi(argv[1]);
    Ny = atoi(argv[2]);
    Nz = atoi(argv[3]);
    lambda = atof(argv[4]);
    rgh    = atof(argv[5]);
    if (rank == 0){
        lambdax = lambda;
        lambday = lambda;

/*  RNG setup  */
        gsl_rng_env_setup();
        T = gsl_rng_default;
        r = gsl_rng_alloc (T);
        printf ("r is a '%s' generator\n", gsl_rng_name (r));
        for (i=0; i < 1000000; i++)
            ran = gsl_rng_uniform(r);

/*  Fill z  */
        z = (fftw_complex*) fftw_malloc(Nx*Ny*sizeof(fftw_complex));
    
        Nx2 = Nx/2;
        Ny2 = Ny/2;
        Pi  = acos(-1.0);
        norm= sqrt(log(1.0+rgh*rgh));
        norm = norm*sqrt(2.0*Pi*lambdax*lambday);
    
        for (i = 0; i <= Nx2; i++)
        for (j = 0; j <= Ny2; j++){
            kx=2.0*Pi*i*lambdax/Nx;
            ky=2.0*Pi*j*lambday/Ny;
            k2=kx*kx+ky*ky;
            if(j==0 && i==0)  radial = 0.0;
            else              radial = norm/pow((1.0+k2),0.75);

            if (i==0) i0=0;
            else      i0=Nx-i;
            if (j==0) j0=0;
            else      j0=Ny-j;

            if (i%Nx2 + j%Ny2 > 0){
                ran=gsl_rng_uniform(r);
                ph=2*Pi*ran;
            }
            else ph = 0.0;

            n1 =  i*Ny + j;
            n2 = i0*Ny + j0;
            z[n1][0] = radial*cos(ph);
            z[n1][1] = radial*sin(ph);
            z[n2][0] = radial*cos(ph);
            z[n2][1] =-radial*sin(ph);
        }
    
        for (i = 1; i < Nx2; i++)
        for (j = 1; j < Ny2; j++){
            kx=2.0*Pi*i*lambdax/Nx;
            ky=2.0*Pi*j*lambday/Ny;
            k2=kx*kx+ky*ky;
            radial=norm/pow((1.0+k2),0.75);

            i0=Nx-i;
            j0=Ny-j;

            ran=gsl_rng_uniform(r);
            ph=2*Pi*ran;
            n1 =  i*Ny + j0;
            n2 = i0*Ny + j;
            z[n1][0] = radial*cos(ph);
            z[n1][1] = radial*sin(ph);
            z[n2][0] = radial*cos(ph);
            z[n2][1] =-radial*sin(ph);
        }

/*  Compute transform  */
        plan =  fftw_plan_dft_2d(Nx,Ny, z,z, FFTW_BACKWARD, FFTW_ESTIMATE);
        fftw_execute(plan);
        fftw_destroy_plan(plan);

/* Calculate global height array  */
        h = (double *) calloc(Nx*Ny, sizeof(double));
        xnn = Nx*Ny;
        norm =sqrt(xnn);
        for (i = 0; i < Nx*Ny; i++)
            h[i] = exp((z[i][0])/norm);
        avg = 0.0;
        for (i = 0; i < Nx*Ny; i++)
            avg = avg + h[i];
        xnn = Nx*Ny;
        avg = avg/xnn;
        for (i = 0; i < Nx*Ny; i++)
            h[i] = h[i]/avg;
        avg = 0.0;
        var = 0.0;
        for (i = 0; i < Nx*Ny; i++){
            avg = avg + h[i];
            var = var + h[i]*h[i];
        }
        avg = avg/xnn;
        var = var/xnn;
        var = sqrt(var - avg*avg);
        printf("Mean % .3e, Variance % .3e\n", avg, var);
    }

/*  Distribute global height data to slaves  */
    parms.Nx = Nx;  parms.Ny = Ny;    parms.Nz = Nz;
    parms.rank = rank-1; parms.size = size-1;        /* Slave COMM */
    parms = mapping(parms);
    hlocal = (double *) calloc(parms.nx*parms.ny, sizeof(double));

    if (rank == 0){                        /* rank == 0 */
        for (irank = 0; irank < parms.size; irank++){
            parms.rank = irank;
            parms = mapping(parms);
            for (i = 0; i < parms.nx; i++)
            for (j = 0; j < parms.ny; j++){
                i0 = parms.x0 + i;        /* Global coords */
                j0 = parms.y0 + j;
                n1 = i*parms.ny + j;        /* Local index */
                n2 = i0*parms.Ny + j0;        /* Global index */
                hlocal[n1] = h[n2];
            }
            n1 = parms.nx*parms.ny;
            MPI_Send(hlocal, n1, MPI_DOUBLE, irank+1, 0, MPI_COMM_WORLD);
        }
    }

/* Slaves receive data from master node */
    else{                        /* rank != 0 */
        n1 = parms.nx*parms.ny;            /* Fill local height */
        MPI_Recv(hlocal, n1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status);

        sprintf(filename, "%s.%02d", "phi.dat", parms.rank);
        file = fopen(filename, "w");
        fprintf(file, "%5d %5d %5d %5d %5d    0    0.0\n", parms.Nx, parms.Ny, parms.nx, parms.ny, parms.Nz);
        for (i = 0; i < parms.nx; i++)
        for (j = 0; j < parms.ny; j++){
            fprintf(file, " % .5e", hlocal[i*parms.ny+j]);
            if ((i*parms.ny+j+1)%10 == 0) fprintf(file, "\n");
        }
    fclose(file);
    }
    
    MPI_Finalize();
    return 0;
}
