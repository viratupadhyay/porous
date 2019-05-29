#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "proto.h"
#include "parms.h"
#include "mpi.h"

extern struct parms fluxV(double *data, int *sp, struct parms parms)
{
    int i, j, k, d, n, nc, np, nm;
    double Kp, Km, *v[Dim];
    double *tphi = padding(data, sp, 0, parms);
    double *tp = padding(data, sp, 1, parms);

    for(d = 0; d < Dim; d++){
        v[d] = (double *) calloc(parms.N, sizeof(double));
    }

/*  Interior nodes 1 <= x < Nx-1 */
    for(i = 0; i < parms.nx; i++){
        for(j = 0; j < parms.ny; j++){
            for(k = 0; k < parms.Nz; k++){
                for(d = 0; d < Dim; d++){
                    n = (((i+1)*(parms.ny+2)+j+1)*parms.Nz)+k;
                    nc = sp[n];                                          // Center
                    switch (d){
                        case 0: 
                            nm = sp[(n-((parms.ny+2)*parms.Nz))];        // Right/Left
                            np = sp[(n+((parms.ny+2)*parms.Nz))];
                            break;
                        case 1:
                            nm = sp[(n-parms.Nz)];                        // Up/Down
                            np = sp[(n+parms.Nz)];
                            break;
                        case 2:
                            nm = sp[(n-1)];
                            np = sp[(n+1)];
                            break;
                    }
                    if (parms.x0+i == 0 || parms.x0+i == parms.Nx-1) continue;
                    Km = 0.5*(K(tphi[nm])+K(tphi[nc]));
                    Kp = 0.5*(K(tphi[np])+K(tphi[nc]));
                    v[d][((i*parms.ny+j)*parms.Nz)+k] = 0.5*((Km*(tp[nc] - tp[nm]))+(Kp*(tp[np] - tp[nc])));
                }
            }
        }
    }
    

/*  Boundary nodes x = 1, x = Nx-1 */
    if (parms.x0 == 0){
        for(j = 0; j < parms.ny; j++){
            for(k = 0; k < parms.Nz; k++){
                n = ((1*(parms.ny+2)+j+1)*parms.Nz)+k;
                nc = sp[n];
                np = sp[(n+((parms.ny+2)*parms.Nz))];
                Kp = 0.5*(K(tphi[np])+K(tphi[nc]));
                v[0][((0*parms.ny+j)*parms.Nz)+k] = Kp*(tp[np] - tp[nc]);
                for(d = 1; d < Dim; d++){
                    v[d][((0*parms.ny+j)*parms.Nz)+k] = 0;
                }
            }
        }
    }
    if (parms.x1 == parms.Nx-1){
        for(j = 0; j < parms.ny; j++){
            for(k = 0; k < parms.Nz; k++){
                n = (((parms.nx)*(parms.ny+2)+j+1)*parms.Nz)+k;
                nc = sp[n];
                nm = sp[(n-((parms.ny+2)*parms.Nz))];
                Km = 0.5*(K(tphi[nm])+K(tphi[nc]));
                v[0][(((parms.nx-1)*parms.ny+j)*parms.Nz)+k] = Km*(tp[nc] - tp[nm]);
                for(d = 1; d < Dim; d++){
                    v[d][(((parms.nx-1)*parms.ny+j)*parms.Nz)+k] = 0;
                }
            }
        }
    }

    for(d = 0; d < Dim; d++){
        writedata(v[d], data, d+2, 0, parms);
        free (v[d]);
    }
    double Vavg;
    if (parms.V0 == 0){
        Vavg = 1.0;
    }
    else{
        Vavg = parms.V0;
#ifdef QCON
        parms.avgV = avgf(data, 2, parms);
        Vavg = parms.avgV;
#endif
    }
    for(i = 0; i < parms.nx; i++)
    for(j = 0; j < parms.ny; j++){
        for(k = 0; k < parms.Nz; k++){
            n = ((i*parms.ny+j)*parms.Nz)+k;
            for(d = 0; d < Dim; d++){
                data[n*parms.Nvar+2+d] /= Vavg;
            }
        }
    }
    free(tphi); free(tp);
    return parms;
}


void fluxC(double *data, double *kr, double *r, double *dr, struct parms parms)
{
    int i, j, k, n;
    double Da, c0;
    double Dsh = SH/parms.Pe;
    for(i = 0; i < parms.nx; i++){
        for(j = 0; j < parms.ny; j++){
            for(k = 0; k < parms.Nz; k++){
                n = ((i*parms.ny+j)*parms.Nz)+k;
#ifdef MIXED
                Da = parms.Da*kr[n]*data[n*parms.Nvar]/Dsh;
                c0 = data[n*parms.Nvar+4]/(1.0 + Da);
                r[n] = Dsh*S(data[n*parms.Nvar])*(data[n*parms.Nvar+4] - c0)/data[n*parms.Nvar];
                dr[n] = Dsh*S(data[n*parms.Nvar])*(1.0 - 1.0/(1.0 + Da))/data[n*parms.Nvar];    // Mixed case
#endif
#ifdef TLIM
                r[n] = Dsh*S(data[n*parms.Nvar])*data[n*parms.Nvar+Dim+2]/data[n*parms.Nvar];
                dr[n] = Dsh*S(data[n*parms.Nvar])/data[n*parms.Nvar];
#endif
#ifdef RLIM
                r[n] = parms.Da*kr[n]*S(data[n*parms.Nvar])*data[n*parms.Nvar+Dim+2];
                dr[n] = parms.Da*kr[n]*S(data[n*parms.Nvar]);
#endif
            }
        }
    }
}
