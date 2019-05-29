/* IO:    Read/Write routines
        mapping        (parms.Nprocessor mapping)
        readp        (Read parameters from h.ini - formatted)
        readf        (Read height profile - formatted)
        readb        (Read height profile - binary)
        writef        (Write field - formatted)
        writeb        (Write field - binary)
        filename    (Constructs file name) */

#include <string.h>
#include <math.h>
#include "includes.h"

char msg[128];

extern struct parms readf(char *fname, double *a, struct parms parms)
{
    int  x, y, z;
    FILE *file;

    file = fopen(fname, "r");
    if (file == 0){
        sprintf(msg, "file %s is missing: EXIT", fname);
    }
    else {
        fscanf(file, "%5d %5d %5d %5d %5d %5d %lf", &parms.Nx, &parms.Ny, &parms.nx, &parms.ny, &parms.Nz, &parms.cyc, &parms.t);
        for (x = 0; x < parms.nx; x++){
            for (y = 0; y <  parms.ny; y++){
                for (z = 0; z <  parms.Nz; z++){
                    fscanf(file, "%le", &a[(x*parms.ny+y)*parms.Nz+z]);
                }
            }
        }
        fclose(file);
    }

    return parms;
}


extern struct parms readb(char *fname, double *a, struct parms parms)
{
    FILE *file;

    file = fopen(fname, "rb");
    if (file == 0){
        sprintf(msg, "file %s is missing: EXIT", fname);
    }
    else {
        fread(&parms, sizeof(parms), 1, file);		
        fread(a, sizeof(double), parms.nx*parms.ny*parms.Nz, file);
        fclose(file);
    }

    return parms;
}


void filename(char *fname, char *name, char *workdir, struct parms parms)
{
    int i;
    char id[4];
    strcpy(fname, name);
    if (workdir != '\0')
    {
        for (i = strlen(fname)+1; i > 0; i--)
            fname[i+strlen(workdir)]=fname[i-1];
        fname[strlen(workdir)] = '/';
        for (i = 0; i < strlen(workdir); i++)
            fname[i] = workdir[i];
    }
    if (name[strlen(name)-1] == '.'){
        sprintf(id, "%04d.%02d", parms.cyc, parms.rank);
        strcat(fname, id);
    }
    else{
        sprintf(id, ".%02d", parms.rank);
        strcat(fname, id);
    }
}


void nameinput(char *fname, char *name, char *workdir)  // Input file without rank in the name
{
    int i;
    strcpy(fname, name);
    if (workdir != '\0')
    {
        for (i = strlen(fname)+1; i > 0; i--)
            fname[i+strlen(workdir)]=fname[i-1];
        fname[strlen(workdir)] = '/';
        for (i = 0; i < strlen(workdir); i++)
            fname[i] = workdir[i];
    }
}


void writef(char *fname, double *data, int dtype, struct parms parms, int freq)
{
    FILE *file;
    int i, j, k;

    if (freq == 0)  return;               // No output
    if (parms.cyc%freq != 0) return;      // No output this cycle

    file = fopen(fname, "w");
    fprintf(file, "%5d %5f %5d %5d %5d %5d %5d %5d %5d\n", parms.cyc, parms.t, parms.nx, parms.ny, parms.Nz, parms.Np, parms.Nq, parms.p, parms.q);
    for (i = 0; i < parms.nx; i++){
        for (j = 0; j < parms.ny; j++){
            for (k = 0; k < parms.Nz; k++){
                fprintf(file, " % .5e", data[((i*parms.ny+j)*parms.Nz+k)*parms.Nvar+dtype]);
                if (((i*parms.ny+j)*parms.Nz+k+1)%10 == 0) fprintf(file, "\n");
            }
        }
    }
    fclose(file);
}


void writeb(char *fname, double *data, int dtype, struct parms parms)
{
    FILE *file;
    int i;
    double *a;
    
    file = fopen(fname, "wb");
    fwrite(&parms, sizeof(parms), 1, file);
    a = (double *) calloc(parms.N, sizeof(double));
    for (i = 0; i < parms.N; i++){
        a[i] = data[(i*parms.Nvar)+dtype];
    }
    fwrite(a, sizeof(double), parms.N, file);
    
    free(a);
    fclose(file);
}
