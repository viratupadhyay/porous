#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "proto.h"
#include <math.h>
#include "parms.h"

//#define Dm = 1.0
//#define lambdaL = 0.5
//#define lambdaT = 0.1
//#define alp0 = 0.5
//#define v0 = 2.0

// Permeability model

double K(double phi)
{
    return (phi*phi*phi);
}

// Area model
double S(double phi)
{
    return 1.0;
    //return (pow((AMAX-phi+0.01),SAF));
    //return (pow(((AMAX-phi)/(AMAX+1)),SAF));
}

// Dispersion model

double *D(double *data, int n, struct parms parms)
{
    int i, j;
    double *DC, De, Dl, Dt, modv = 0;
    DC = (double *) calloc(Dim*Dim, sizeof(double));
    double Dm = 1.0, lambdaL = 0.5, lambdaT = 0.1, alp0 = 0.5, v0 = 2.0;
    De = alp0*Dm*data[n*parms.Nvar];
    Dl = lambdaL*v0*data[n*parms.Nvar];
    Dt = lambdaT*v0*data[n*parms.Nvar];
    // Scale velocities
    //for(i = 0; i < Dim; i++){
    //    data[n*parms.Nvar+2+i] = data[n*parms.Nvar+2+i]/v0;
    //    modv += data[n*parms.Nvar+2+i]*data[n*parms.Nvar+2+i];
    //}
    modv = pow(modv,0.5);

    // Dispersion tensor
    for(i = 0; i < Dim; i++){
        for(j = 0; j < Dim; j++){
            if(i != j){
                DC[i*Dim+j] = 0;
                //DC[i*Dim+j] = (Dl - Dt)*data[n*parms.Nvar+2+i]*data[n*parms.Nvar+2+j]/modv;
            }
            else{
                DC[0] = De + Dl;
                DC[Dim+1] = De + Dt;
                //DC[0] = De + Dl*data[n*parms.Nvar+2];
                //DC[Dim+1] = De + Dt*data[n*parms.Nvar+2];
            }
        }
    }
    //for(i = 0; i < Dim*Dim; i++){
        //        printf("%4f\n",DC[i]);
        //}
    return(DC);
}
