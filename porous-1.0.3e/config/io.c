/* IO:	Read/Write routines
		mapping		(parms.Nprocessor mapping)
		readp		(Read parameters from h.ini - formatted)
		readf		(Read height profile - formatted)
		readb		(Read height profile - binary)
		writef		(Write field - formatted)
		writeb		(Write field - binary)
		filename	(Constructs file name) */

#include <string.h>
#include <math.h>
#include "includes.h"

char msg[128];

extern struct parms mapping(struct parms parms)
{
	int nx, ny, psave, p1, p2;
	int p, q;	// parms.Nproc coordinates in domain
	double fit;
	fit = fabs(parms.Ny - parms.size*parms.Nx);  psave = 1;
	for (p1 = 2; p1 <= parms.size; p1++){			// Check layout
		if (parms.size%p1 == 0){
			p2 = p1*p1;				// Valid layout
			if (fabs(parms.Ny - parms.size*parms.Nx/p2) < fit){
				psave = p1;
				fit = fabs(parms.Ny - parms.size*parms.Nx/p2);
			}
		}
	}
	parms.Np = psave; parms.Nq = parms.size/psave;		// Best layout

	nx = (parms.Nx-1)/parms.Np + 1;		// Domain size
	ny = (parms.Ny-1)/parms.Nq + 1;
	if (nx <= 1)  sprintf(msg, "Invalid proc layout %d %d\n", nx, ny);
	if (ny <= 1)  sprintf(msg, "Invalid proc layout %d %d\n", nx, ny);
	if ((parms.rank+1)%parms.Np == 0)  nx = parms.Nx - (parms.Np-1)*nx;
	if ((parms.rank+1)%parms.Nq == 0)  ny = parms.Ny - (parms.Nq-1)*ny;
	parms.nx = nx;  parms.ny = ny;
	parms.N = parms.nx*parms.ny*parms.Nz;
	
	// Calculate the rank of neighbors using the proc coordinates
	p = parms.rank/parms.Nq;
	q = parms.rank%parms.Nq;
	if(p == 0 && q == 0){
		parms.L = (parms.Np-1)*parms.Nq+q;
		parms.D = p*parms.Nq+parms.Nq-1;
	}
	else if (p == 0 && q != 0){
		parms.L = (parms.Np-1)*parms.Nq+q;
		parms.D = p*parms.Nq+q-1;
	}
	else if(p != 0 && q == 0){
		parms.L = (p-1)*parms.Nq+q;
		parms.D = p*parms.Nq+parms.Nq-1;
	}
        else{
		parms.L = (p-1)*parms.Nq+q;
		parms.D = p*parms.Nq+q-1;
	}
	
	parms.R = ((p+1)%parms.Np)*parms.Nq+q;

	parms.U = p*parms.Nq+(q+1)%parms.Nq;

	//Calculate the extents of the domain
	parms.x0 = (parms.rank/parms.Nq)*nx;	// Domain origin
	parms.y0 = (parms.rank%parms.Nq)*ny;
	parms.x1 = parms.x0 + nx-1;
	parms.y1 = parms.y0 + ny-1;
	return parms;
}


extern struct parms readf(char *name, double *a, struct parms parms)
{
	int  x, y;
	FILE *fileptr;

	fileptr = fopen(name, "r");
	if (fileptr == 0){
		sprintf(msg, "file %s is missing: EXIT", name);
	}
	else {
		fscanf(fileptr, "%5d %5d", &parms.nx, &parms.ny);
		for (x = 0; x < parms.nx; x++)
		for (y = 0; y <  parms.ny; y++)
			fscanf(fileptr, "%le", &a[x*parms.ny+y]);
		fclose(fileptr);
	}

	return parms;
}


char *filename (char *name, struct parms parms)
{
	char id[4];							// Max # digits
	sprintf(id, ".%02d", parms.rank);
	strcat(name, id);
	return(name);
}

