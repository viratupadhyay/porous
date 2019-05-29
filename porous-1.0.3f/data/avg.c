/* Average profiles - Use ../src/parms.h for list */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "mpi.h"
#include "../src/proto.h"
#include "../src/parms.h"

#define MAX 10000

int main(int argc, char ** argv)
{
	double havg[MAX]={MAX*0};
	double pavg[MAX]={MAX*0};
	double cavg[MAX]={MAX*0};
	double dh[MAX]={MAX*0};
	double dp[MAX]={MAX*0};
	double dc[MAX]={MAX*0};
	double *f;
	double hmax, pmax, cmax;
	int  x, y, z, n, xh, xp, xc;
	char hname[64]={"data/h."};
	char pname[64]={"data/p."};
	char cname[64]={"data/c."};
	char aname[64]={"avg."};
	char mname[64]={"max."};
	char fname[64], *wkdir;
	FILE *fileptr;
    
	if (argc > 1) parms.cyc = atoi(argv[1]);
	if (argc > 2) wkdir = argv[2];

#if (H > 0)
	if (parms.cyc%H == 0){
		filename(fname, hname, wkdir, parms);
		parms = readp(fname, parms);
		n = parms.nx*parms.ny*parms.Nz;
		f = (double *) calloc(n, sizeof(double));
		readf(fname, f, parms);
		for (x = 0; x <= parms.nx; x++)
		for (y = 0; y <   parms.ny; y++){
			n = x*parms.ny + y;
			havg[x] += f[n]/(double)(parms.ny);
			dh[x] += f[n]*f[n]/(double)(parms.ny);
		}
	}
	free(f);
#endif

#if (P > 0)
	if (parms.cyc%P == 0){
		filename(fname, pname, wkdir, parms);
		parms = readp(fname, parms);
		n = parms.nx*parms.ny*parms.Nz;
		f = (double *) calloc(n, sizeof(double));
		readf(fname, f, parms);
		for (x = 0; x < parms.nx; x++)
		for (y = 0; y < parms.ny; y++){
            for (z = 0; z < parms.Nz; z++){
			    n = ((x*parms.ny+y)*parms.Nz+z);
			    pavg[x] += f[n]/(double)(parms.ny*parms.Nz);
			    dp[x] += f[n]*f[n]/(double)(parms.ny*parms.Nz);
		    }
        }
	
    }
	free(f);
#endif

#if (C > 0)
	if (parms.cyc%C == 0){
		filename(fname, cname, wkdir, parms);
		parms = readp(fname, parms);
		n = (parms.nx+1)*parms.ny;
		f = (double *) calloc(n, sizeof(double));
		readf(fname, f, parms);
		for (x = 0; x <= parms.nx; x++)
		for (y = 0; y <  parms.ny; y++){
			n = x*parms.ny + y;
			cavg[x] += f[n]/(double)(parms.ny);
			dc[x] += f[n]*f[n]/(double)(parms.ny);
		}
	}
	free(f);
#endif

	filename(fname, aname, wkdir, parms);
	fileptr = fopen(fname, "w");
	fprintf(fileptr, "    x     havg       pavg       cavg        dh         dp         dc\n");
	for (x = 0; x <= parms.nx; x++){
		dh[x] -= havg[x]*havg[x];
		dp[x] -= pavg[x]*pavg[x];
		dc[x] -= cavg[x]*cavg[x];
		fprintf(fileptr, "%6d % .3e % .3e % .3e % .3e % .3e % .3e\n", x, havg[x], pavg[x], cavg[x], dh[x], dp[x], dc[x]);
	}
	fclose(fileptr);

	hmax = pmax = cmax = 0.0;  xh = xp = xc = 0;			/* Find max fluctuations */
	for (x = 0; x <= parms.nx; x++){
		if (dh[x] > hmax){
			hmax = dh[x]; xh = x;
		}
		if (dp[x] > pmax){
			pmax = dp[x]; xp = x;
		}
		if (dc[x] > cmax){
			cmax = dc[x]; xc = x;
		}
	}
	
	filename(fname, mname, wkdir, parms);
	fileptr = fopen(fname, "w");
	fprintf(fileptr, "   time        xh    hmax        xp    pmax        xc    cmax\n");
	fprintf(fileptr, "% .3e %6d % .3e %6d % .3e %6d % .3e\n", parms.t, xh, hmax, xp, pmax, xc, cmax);
	fclose(fileptr);

	return 0;
}
