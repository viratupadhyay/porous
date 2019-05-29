#include "includes.h"
#include "parms.h"
#include <string.h>

int main (int argc, char *argv[])
{
    struct parms parms;
    struct timeval start;
    int i, rank, size, ncy, cyc=0;
    int  freq[5]={PHI,P,V,C};
    double *data, *phi0, *r, *dr, *phit, *dphi;
    double telp, time, tcyc, havg, vavg, delh;
    char *wkdir=0;
    char input[64]={"init/input.dat"};
    char phichek[64]={"init/phi.chk"};
    char phii[64] = {"data/phi."};
    char p[64] = {"data/p."};
    char vx[64] = {"data/vx."};
    char vy[64] = {"data/vy."};
    char c[64] = {"data/c."};
    char avname[64]={"data/avg."};
    char fname[64];
    FILE *fileptr;

    if (argc > 1) cyc = atoi(argv[1]);
    if (argc > 2) wkdir = argv[2];

    winit(&start);
    // MPI Initialization
    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    parms.rank = rank; parms.size = size;

    parms.V0 = 0;
    parms.Nvar = (3+Dim);
    parms.nsten = (2*Dim)+1;
    nameinput(fname, input, wkdir);
    fileptr = fopen(fname, "r");            // Read input data
    fscanf(fileptr, "%*s %*s %*s %*s %*s %*s %*s");
    fscanf(fileptr, "%d %lf %lf %lf %lf %lf %lf", &ncy, &parms.Re, &parms.Pe, &parms.Da, &parms.pow, &parms.dx, &parms.dt);
    fclose(fileptr);

    // Reading the initial porosity profile
#ifdef SEED
    char phseed[64] = {"init/phi.dat-seed"};
    //sprintf(phseed, "seed%f/phi.dat", scale);
    filename(fname, phseed, wkdir, parms);
#endif
#ifndef SEED
    char phinit[64] = {"init/phi.dat"};
    filename(fname, phinit, wkdir, parms);
#endif
    fileptr = fopen(fname, "r");
    fscanf(fileptr, "%5d %5d %5d %5d %5d %5d %lf", &parms.Nx, &parms.Ny, &parms.nx, &parms.ny, &parms.Nz, &parms.cyc, &parms.t);
    fclose(fileptr);

    // 2D Domain decomposition
    parms = mapping(parms);
    parms.N = (parms.nx*parms.ny*parms.Nz);
    
    // Initialize pointers
    phi0 = (double *) calloc(parms.N, sizeof(double));
    data = (double *) calloc(parms.N*parms.Nvar, sizeof(double));
    dphi = (double *) calloc(parms.N, sizeof(double));
    phit = (double *) calloc(parms.N, sizeof(double));
    r = (double *) calloc(parms.N, sizeof(double));
    dr = (double *) calloc(parms.N, sizeof(double));
    
    // Create Spiral indexing
    int *sp = spiral(parms);
    
    readf(fname, phi0, parms);                  // Read from the input porosity profile
    writedata(phi0, data, 0, 0, parms);        // Write the input porosity to pointer "data" @ 0

    for(i = 0; i < parms.N; i++){
        data[(i*parms.Nvar)+4] = 1.0;          // Initialize c = 1.0
    }
    erode(data, r, dr, sp, parms);              // Intial fields
    parms.V0 = avgf(data, 2, parms);

    if (cyc > 0){
        parms.cyc = cyc;
        filename(fname, phichek, wkdir, parms);
        readb(fname, phi0, parms);
    }

    writedata(phi0, data, 0, 0, parms);        // Write the input porosity to pointer "data" @ 0
    telp = wtime(&start);
    time = telp;
    if (parms.rank == 0){
        fprintf(stdout, "Initialization; elapsed time: % .3e\n\n", telp);
    }
    while (parms.cyc <= ncy){

        erode(data, r, dr, sp, parms);

        //if (parms.rank == 0){
        havg = avgf(data, 0, parms);
        vavg = avgf(data, 2, parms);
        vavg = vavg*parms.Ny;
        delh = havg - parms.havg;
        parms.havg = havg;
        if (parms.rank == 0){
            filename(fname, avname, wkdir, parms);
            fileptr = fopen(fname, "w");
            if (parms.cyc == 0) fprintf(fileptr, " cycle   time       havg       qavg       delh\n");
            fprintf(fileptr, "%5d % .3e % .3e % .3e % .3e\n", parms.cyc, parms.t, havg, vavg, delh);
            fclose(fileptr);
        }
        filename(fname, phii, wkdir, parms);    // Write data files
        writef(fname, data, 0, parms, freq[0]);
        filename(fname, p, wkdir, parms);
        writef(fname, data, 1, parms, freq[1]);
        filename(fname, vx, wkdir, parms);
        writef(fname, data, 2, parms, freq[2]);
        filename(fname, vy, wkdir, parms);
        writef(fname, data, 3, parms, freq[2]);
        filename(fname, c, wkdir, parms);
        writef(fname, data, 4, parms, freq[3]);

        // Calculate erosion during dt
        double dt = parms.dt/parms.Da;
#ifdef EULER
        for (i = 0; i < parms.N; i++)
            dphi[i] = r[i]*dt;
#endif

#ifdef MIDPT
        for (i = 0; i < parms.N; i++){
            phit[i] = phi0[i] + 0.5*r[i]*dt;
        }
        writedata(phit, data, 0, 0, parms);
        erode (data, r, dr, sp, parms);
        for (i = 0; i < parms.N; i++)
            dphi[i] = r[i]*dt;
#endif

#ifdef RK4
        for (i = 0; i < parms.N; i++){
            phit[i] = phi0[i] + 0.5*r[i]*dt;
            dphi[i] = r[i]*dt/6.0;
        }
        writedata(phit, data, 0, 0, parms);
        erode (data, r, dr, sp, parms);
        for (i = 0; i < parms.N; i++){
            phit[i] = phi0[i] + 0.5*r[i]*dt;
            dphi[i]+= r[i]*dt/3.0;
        }
        writedata(phit, data, 0, 0, parms);
        erode (data, r, dr, sp, parms);
        for (i = 0; i < parms.N; i++){
            phit[i] = phi0[i] + r[i]*dt;
            dphi[i]+= r[i]*dt/3.0;
        }
        writedata(phit, data, 0, 0, parms);
        erode (data, r, dr, sp, parms);
        for (i = 0; i < parms.N; i++)
            dphi[i]+= r[i]*dt/6.0;
#endif

        // Update height profile
        for (i = 0; i < parms.N; i++)
            phi0[i] += dphi[i];
        writedata(phi0, data, 0, 0, parms);

        telp = wtime(&start);
        tcyc = time;
        time += telp;
        tcyc = time - tcyc;
        if (parms.rank == 0){
            fprintf(stdout, "Cycle %6d;   cycle time: % .3e  total time: % .3e\n\n", parms.cyc, tcyc, time);
        }
        parms.cyc++;
        parms.t += parms.dt;
    
#if CHEK > 0
        if ((parms.cyc-1)%CHEK == 0){
            filename(fname, phichek, wkdir, parms);
            writeb(fname, data, 0, parms);            /* Checkpoint aperture */
        }
#endif
    }
    // Free memory
    free(phit); free(r); free(dr);
    free (dphi); free(phi0);
    free(data); free(sp);
    // Finalize MPI
    MPI_Finalize();
    return (0);
}
