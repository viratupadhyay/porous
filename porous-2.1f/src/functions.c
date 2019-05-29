#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "proto.h"
#include "mpi.h"
#include "parms.h"

char msg[128];

extern struct parms mapping(struct parms parms)
{
    int nx, ny, psave, p1, p2;
    int p, q;    // parms.Nproc coordinates in domain
    double fit;
    fit = fabs(parms.Ny - parms.size*parms.Nx);  psave = 1;
    for (p1 = 2; p1 <= parms.size; p1++){            // Check layout
        if (parms.size%p1 == 0){
            p2 = p1*p1;                // Valid layout
            if (fabs(parms.Ny - parms.size*parms.Nx/p2) < fit){
                psave = p1;
                fit = fabs(parms.Ny - parms.size*parms.Nx/p2);
            }
        }
    }
    parms.Np = psave; parms.Nq = parms.size/psave;        // Best layout

    nx = (parms.Nx-1)/parms.Np + 1;        // Domain size
    ny = (parms.Ny-1)/parms.Nq + 1;
    if ((parms.rank+1)%parms.Np == 0)  nx = parms.Nx - (parms.Np-1)*nx;
    if ((parms.rank+1)%parms.Nq == 0)  ny = parms.Ny - (parms.Nq-1)*ny;
    if (nx <= 1)  sprintf(msg, "Invalid proc layout %d %d\n", nx, ny);
    if (ny <= 1)  sprintf(msg, "Invalid proc layout %d %d\n", nx, ny);
    parms.nx = nx;  parms.ny = ny;
    parms.N = parms.nx*parms.ny*parms.Nz;
    
    // Calculate the rank of neighbors using the proc coordinates
    p = parms.rank/parms.Nq;
    q = parms.rank%parms.Nq;
    parms.p = p;
    parms.q = q;

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
    
    
    //parms.L = ((p-1)%parms.Np)*parms.Nq+q;

    //parms.D = p*parms.Nq+(q-1)%parms.Nq;
    
    parms.R = ((p+1)%parms.Np)*parms.Nq+q;

    parms.U = p*parms.Nq+(q+1)%parms.Nq;

    //Calculate the extents of the domain
    parms.x0 = (parms.rank/parms.Nq)*nx;    // Domain origin
    parms.y0 = (parms.rank%parms.Nq)*ny;
    parms.x1 = parms.x0 + nx-1;
    parms.y1 = parms.y0 + ny-1;
    return parms;
}

// Create spiral indexing
#define valid(i, j) 0 <= i && i < m && 0 <= j && j < n && !s[i][j]      
int *spiral(struct parms parms)
{
    int i, j, k, m, n;
    m = parms.nx;
    n = parms.ny;
    int **s = calloc(1, sizeof(int *) * m + sizeof(int) * m * n);
    int *sp = malloc(sizeof(int)*((m+2)*(n+2)*parms.Nz));
    s[0] = (int*)(s + m);
    for (i = 1; i < m; i++) s[i] = s[i - 1] + n;
    int dx = 1, dy = 0, val = 0, t;
    for (i = j = 0; valid(i, j); i += dy, j += dx ) {
        for (; valid(i, j); j += dx, i += dy)
            s[i][j] = ++val;
        j -= dx; i -= dy;
        t = dy; dy = dx; dx = -t;
    }
    for(i = 1; i < m+1; i++){
        for(j = 1; j < n+1; j++){
            for(k = 0; k < parms.Nz; k++){
                sp[(i*(n+2)+j)*parms.Nz+k] = (m*n-s[i-1][j-1])*parms.Nz+k;
            }
        }
    }
    for (j = 0; j < n+2; j++){
        for(k = 0; k < parms.Nz; k++){
            sp[(j)*parms.Nz+k] = ((m*n)+j)*parms.Nz+k;
            sp[((m+1)*(n+2)+j)*parms.Nz+k] = ((m+2)*(n+1)+1-j)*parms.Nz+k;
        }
    }
    for (i = 1; i < m+2; i++){
        for(k = 0; k < parms.Nz; k++){
            sp[(i*(n+2))*parms.Nz+k] = ((m+2)*(n+2)-i)*parms.Nz+k;
            sp[((i+1)*(n+2)-1)*parms.Nz+k] = ((m+1)*(n)+1+i)*parms.Nz+k;
        }
    }
    free(s);
    return (sp);
}

// Create padding for a particular variable in data
double *padding(double *data, int *sp, int dtype, struct parms parms)
{
    MPI_Status status;
    int i, j, k, spij;
    int m, n, p, inn;
    m = parms.nx;
    n = parms.ny;
    p = parms.Nz;
    double *temp;
    temp = (double *) calloc((parms.nx+2)*(parms.ny+2)*parms.Nz, sizeof(double));
    // Building temp array with padding in spiral indexing format
    for (i = 0; i < parms.nx; i++){
        for (j = 0; j < parms.ny; j++){
            for (k = 0 ; k < parms.Nz; k++){
                spij = sp[((i+1)*(parms.ny+2)+j+1)*parms.Nz+k];
                temp[spij] = data[((i*parms.ny+j)*parms.Nz+k)*parms.Nvar+dtype];
            }
        }
    }
    
    // MPI communication calls

    // Order of calls is Down-Right-Up-Left-Odd (D.R.U.L.O.)

    inn = sp[(2*(n+2)+1)*p];
    double *add_send[5] = {&temp[inn], &temp[inn+(m-2)*p], &temp[inn+(m+n-3)*p],&temp[inn+(2*m+n-4)*p],&temp[inn+(2*(m+n)-5)*p]};
    double *add_recv[5] = {&temp[(n*(m+1)+3)*p], &temp[(m*n+1)*p], &temp[((m+1)*(n+1)+n+3)*p], &temp[((m+1)*(n+1)+2)*p], &temp[(n*(m+1)+2)*p]};

    // Down
    MPI_Sendrecv(add_send[0], (m-1)*p, MPI_DOUBLE, parms.D, 1, add_recv[0], (m-1)*p, MPI_DOUBLE, parms.U, 1, MPI_COMM_WORLD, &status);

    // Right
    MPI_Sendrecv(add_send[1], n*p, MPI_DOUBLE, parms.R, 2, add_recv[1], n*p, MPI_DOUBLE, parms.L, 2, MPI_COMM_WORLD, &status);

    // Up
    MPI_Sendrecv(add_send[2], m*p, MPI_DOUBLE, parms.U, 3, add_recv[2], m*p, MPI_DOUBLE, parms.D, 3, MPI_COMM_WORLD, &status);

    // Left
    MPI_Sendrecv(add_send[3], n*p, MPI_DOUBLE, parms.L, 4, add_recv[3], n*p, MPI_DOUBLE, parms.R, 4, MPI_COMM_WORLD, &status);

    // Odd
    MPI_Sendrecv(add_send[4], p, MPI_DOUBLE, parms.D, 5, add_recv[4], p, MPI_DOUBLE, parms.U, 5, MPI_COMM_WORLD, &status);

    return(temp);
}

double *matP(double *data, int *sp, struct parms parms)
{
    int i, j, k, l, stencil[parms.nsten], n, no;
    int Nmat = parms.N*(parms.nsten+1);
    double *tphi = padding(data, sp, 0, parms);        // Padding porosity field
    double *Ab = (double *) calloc(Nmat, sizeof(double));
    double p_in = 0.0;
    double p_out = 1.0*parms.Nx-1;
    for(i = 0; i < parms.nx; i++){
        for(j = 0; j < parms.ny; j++){
            for(k = 0; k < parms.Nz; k++){
                no = ((k*parms.ny+j)*parms.nx)+i;
                n = (((i+1)*(parms.ny+2)+j+1)*parms.Nz)+k;
                stencil[0] = sp[n];
                stencil[1] = sp[n-((parms.ny+2)*parms.Nz)];    // Left neighbor
                stencil[2] = sp[n-parms.Nz];                // Down neighbor
                stencil[3] = sp[n+parms.Nz];                // Up neighbor
                stencil[4] = sp[n+((parms.ny+2)*parms.Nz)];    // Right neighbor
#if Dim == 3
                stencil[5] = sp[n-k+(k+parms.Nz-1)%parms.Nz];
                stencil[6] = sp[n-k+(k+parms.Nz+1)%parms.Nz];
#endif
                if (parms.x0+i == 0){
                    Ab[no*parms.nsten] = 1;
                    Ab[parms.N*parms.nsten+no] = p_in;        // RHS
                }
                else if(parms.x0+i == parms.Nx-1){
                    Ab[no*parms.nsten] = 1;
                    Ab[parms.N*parms.nsten+no] = p_out;        // RHS
                }
                else{
                    for(l = 1; l < parms.nsten; l++){
                        Ab[no*parms.nsten+l] = 0.5*(K(tphi[stencil[0]])+K(tphi[stencil[l]]));
                        Ab[no*parms.nsten] -= Ab[no*parms.nsten+l];
                        Ab[parms.N*parms.nsten+no] = 0.0;    // RHS
                    }
                }
            }
        }
    }

    free(tphi);
    return (Ab);
}


double *matC(double *data, int *sp, double *r, double *dr, struct parms parms)
{
    int i, j, k, l, n, no, np, stencil[parms.nsten];
    int Nmat = parms.N*(parms.nsten+1);
    double *tphi = padding(data, sp, 0, parms);
    double *Ab = (double *) calloc(Nmat, sizeof(double));
    double c_in = 1.0, Vx, Vy;
    for(i = 0; i < parms.nx; i++){
        for(j = 0; j < parms.ny; j++){
            for(k = 0; k < parms.Nz; k++){
                n = ((i*parms.ny+j)*parms.Nz)+k;                    // Forward index
                no = ((k*parms.ny+j)*parms.nx)+i;                    // Reverse index
                np = (((i+1)*(parms.ny+2)+j+1)*parms.Nz)+k;            // Padding index
                Ab[parms.N*parms.nsten+no] = dr[n]*data[n*parms.Nvar+Dim+2] - r[n];
                Vx = data[n*parms.Nvar+2]/parms.dx;
                Vy = data[n*parms.Nvar+3]/parms.dx;
                if (parms.x0+i == 0){
                    Ab[no*parms.nsten] = 1;
                    Ab[parms.N*parms.nsten+no] = c_in;
                }
                else if (parms.x0+i == parms.Nx-1){
                    Ab[no*parms.nsten] = Vx + dr[n];
                    Ab[no*parms.nsten+1] = -Vx;                        // Left neighbor
                }
                else{
                    Ab[no*parms.nsten] = dr[n];
                    if (Vx < 0.0){
                        Ab[no*parms.nsten+4] = Vx;
                        Ab[no*parms.nsten] -= Ab[no*parms.nsten+4];
                    }
                    else{
                        Ab[no*parms.nsten+1] = -Vx;
                        Ab[no*parms.nsten] -= Ab[no*parms.nsten+1];
                    }
                    if (Vy < 0.0){
                        Ab[no*parms.nsten+3] = Vy;
                        Ab[no*parms.nsten] -= Ab[no*parms.nsten+3];
                    }
                    else{
                        Ab[no*parms.nsten+2] = -Vy;
                        Ab[no*parms.nsten] -= Ab[no*parms.nsten+2];
                    }
#if Dim == 3
                    double Vz = data[n*parms.Nvar+4]/parms.dx;
                    if (Vz < 0.0){
                        Ab[no*parms.nsten+6] = Vz;
                        Ab[no*parms.nsten] -= Ab[no*parms.nsten+6];
                    }
                    else{
                        Ab[no*parms.nsten+5] = -Vz;
                        Ab[no*parms.nsten] -= Ab[no*parms.nsten+5];
                    }
#endif
                    stencil[0] = sp[np];
                    stencil[1] = sp[np-((parms.ny+2)*parms.Nz)];     // Left neighbor
                    stencil[2] = sp[np-parms.Nz];                    // Down neighbor
                    stencil[3] = sp[np+parms.Nz];                    // Up neighbor
                    stencil[4] = sp[np+((parms.ny+2)*parms.Nz)];     // Right neighbor
#if Dim == 3
                    stencil[5] = sp[np-k+(k+parms.Nz-1)%parms.Nz];
                    stencil[6] = sp[np-k+(k+parms.Nz+1)%parms.Nz];
#endif
                    
                    //parms.DC = D(data, n, parms);                     // Model for Dispersion
                    
                    parms.DC = 1.0/(parms.Pe*parms.dx*parms.dx);
                    for(l = 1; l < parms.nsten; l++){
                        Ab[no*parms.nsten+l] -= parms.DC*(tphi[stencil[0]] + tphi[stencil[l]])/2.0;
                        Ab[no*parms.nsten] += parms.DC*(tphi[stencil[0]] + tphi[stencil[l]])/2.0;
                    }
                }
            }
        }
    }
    free(tphi);
    return(Ab);
}


// Write a set of values to a particular data variable
void writedata(double *values, double *data, int dtype, int order, struct parms parms)
{
    int i, j, k, n;
    for(i = 0; i < parms.nx; i++){
        for(j = 0; j < parms.ny; j++){
            for(k = 0; k < parms.Nz; k++){
                n = ((i*parms.ny+j)*parms.Nz+k);
                if(order == 0){
                    data[(n*parms.Nvar+dtype)] = values[n];
                }
                else if(order == 1){
                    data[(n*parms.Nvar+dtype)] = values[(((k*parms.ny)+j)*parms.nx+i)];
                }
            }
        }
    }
}


double sumf(double *data, int dtype, struct parms parms)
{
    double sum = 0;
    int i, j, k, n;

    for (i = 0; i < parms.nx; i++)
    for (j = 0; j < parms.ny; j++){
        for(k = 0; k < parms.Nz; k++){
            n = ((i*parms.ny+j)*parms.Nz+k);
            sum += data[(n*parms.Nvar+dtype)];
        }
    }
    return sum;
}


double avgf(double *data, int dtype, struct parms parms)
{
    double localsum, globalsum, avg;
    localsum = sumf(data, dtype, parms);
    MPI_Allreduce(&localsum, &globalsum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    avg = globalsum/(double)(parms.Nx*parms.Ny*parms.Nz);
    return avg;
}


void erode(double *data, double *kr, double *r, double *dr, int *sp, struct parms parms)
{
    // Solve for pressure field using HYPRE
    double *AbP = matP(data, sp, parms);        // Create matrix, vector for HYPRE for solving for p
#ifdef DEBUG
    if (parms.rank == 0){
        fprintf(stdout, "Pressure solve: \n");
    }
#endif
    double *xp = solve(AbP, SLN, parms);             // Solve for pressure fields using HYPRE
    writedata(xp, data, 1, 1, parms);           // Write pressure to data @ 1 in reverse order

    // Solve for velocity fields
    parms = fluxV(data, sp, parms);

    // Solve for concentration field using HYPRE
    double *AbC = matC(data, sp, r, dr, parms); // Create matrix, vector for HYPRE for solving for c
#ifdef DEBUG
    if (parms.rank == 0){
        fprintf(stdout, "Concentration solve: \n");
    }
#endif
    double *xc = solve(AbC, SLN, parms);             // Solve for concentration fields using HYPRE
    writedata(xc, data, Dim+2, 1, parms);           // Write concentration to data @ 1 in reverse order
    fluxC(data, kr, r, dr, parms);                  // Update reaction rates

    free(xp); free(xc);
}
