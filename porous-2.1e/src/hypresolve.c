#include <stdlib.h>
#include <stdio.h>
#include "proto.h"
#include "mpi.h"
#include "_hypre_utilities.h"
#include "HYPRE_struct_ls.h"
#include "parms.h"

double *solve(double *Ab, int solver_id, struct parms parms)
{
    int i, j;
    double final_res_norm;
    int time_index, n_pre, n_post, num_iterations;
    n_pre  = 1; n_post = 1;
    double *A_val, *b_val;
    A_val = (double *) calloc(parms.N*parms.nsten, sizeof(double));
    b_val = (double *) calloc(parms.N, sizeof(double));
                    
    for (i = 0; i < (parms.N*parms.nsten); i++){
        A_val[i] = Ab[i];
    }
    for (i = 0; i < parms.N; i++){
        b_val[i] = Ab[i+parms.N*parms.nsten];
    }

    // HYPRE //
    HYPRE_StructGrid     grid;
    HYPRE_StructStencil  stencil;
    HYPRE_StructMatrix   A;
    HYPRE_StructVector   b;
    HYPRE_StructVector   x;
    HYPRE_StructSolver   solver;
    HYPRE_StructSolver   precond;

#if Dim == 2
    HYPRE_Int ilower[2] = {parms.x0, parms.y0};
    HYPRE_Int iupper[2] = {parms.x1, parms.y1};
#endif

#if Dim == 3
    HYPRE_Int ilower[3] = {parms.x0, parms.y0, 0};
    HYPRE_Int iupper[3] = {parms.x1, parms.y1, parms.Nz-1};
#endif
    {
    // Create an empty 2D grid object
        HYPRE_StructGridCreate(MPI_COMM_WORLD, Dim, &grid);

    // Add a new box to the grid
        HYPRE_StructGridSetExtents(grid, ilower, iupper);

    // 1. Set up periodic boundary condition in y-direction and create the grid 
        int pset[3]; 
        pset[0] = 0; pset[1] = parms.Ny; pset[2] = 0;
#if Dim == 3
        pset[2] = parms.Nz;
#endif
    //HYPRE_StructGridSetNumGhost(grid,pset)
        HYPRE_StructGridSetPeriodic(grid, pset);
        HYPRE_StructGridAssemble(grid);
    }

    // 2. Define the discretization stencil
    {
        if (Dim == 2){

        // Create an empty 2D, 5-pt stencil object
            HYPRE_StructStencilCreate(2, parms.nsten, &stencil);

        // Define the geometry of the stencil
            {
                int offsets[5][2] = {{0,0}, {-1,0}, {0,-1}, {0,1}, {1,0}};
                for (i = 0; i < parms.nsten; i++)
                    HYPRE_StructStencilSetElement(stencil, i, offsets[i]);
            }
        }
        else
        {
            HYPRE_StructStencilCreate(3, parms.nsten, &stencil);

            // Define the geometry of the 3D stencil
            {
                int offsets[7][3] = {{0,0,0}, {-1,0,0}, {0,-1,0}, {0,1,0}, {1,0,0}, {0,0,-1}, {0,0,1}};

                for (i = 0; i < parms.nsten; i++)
                    HYPRE_StructStencilSetElement(stencil, i, offsets[i]);
            }
        }
    }
    // 3. Set up a Struct Matrix A from Aval
    {
        HYPRE_Int stencil_indices[parms.nsten];

        // Create an empty matrix object
        HYPRE_StructMatrixCreate(MPI_COMM_WORLD, grid, stencil, &A);

        // Indicate that the matrix coefficients are ready to be set
        HYPRE_StructMatrixInitialize(A);

        for (j = 0; j < parms.nsten; j++)
            stencil_indices[j] = j;

        HYPRE_StructMatrixSetBoxValues(A, ilower, iupper, parms.nsten, stencil_indices, A_val);

        free(A_val);
    }

    // 4. Set up Struct Vectors for b from b_val and set x = 0
    {
        double *values;

        HYPRE_StructVectorCreate(MPI_COMM_WORLD, grid, &b);
        HYPRE_StructVectorCreate(MPI_COMM_WORLD, grid, &x);

        HYPRE_StructVectorInitialize(b);
        HYPRE_StructVectorInitialize(x);

        values = calloc((parms.N), sizeof(double));

        for (i = 0; i < (parms.N); i++)
            values[i] = 0.0;
        HYPRE_StructVectorSetBoxValues(x, ilower, iupper, values);
        HYPRE_StructVectorSetBoxValues(b, ilower, iupper, b_val);

        free(b_val);
        free(values);
    }

    //Finalize the vector and matrix assembly

    HYPRE_StructMatrixAssemble(A);
    HYPRE_StructVectorAssemble(b);
    HYPRE_StructVectorAssemble(x);
#ifdef DEBUG
    //HYPRE_StructMatrixPrint("./poisson.matrix", A, 0);
    //HYPRE_StructVectorPrint("./poisson.rhs", b, 0);
#endif
    // 6. Set up and use a solver (SMG)
    if (solver_id == 0)
    {
        time_index = hypre_InitializeTiming("SMG Setup");
        hypre_BeginTiming(time_index);
        HYPRE_StructSMGCreate(MPI_COMM_WORLD, &solver);
        HYPRE_StructSMGSetMemoryUse(solver, 0);
        HYPRE_StructSMGSetMaxIter(solver, 100);
        HYPRE_StructSMGSetTol(solver, 1.0e-12);
        HYPRE_StructSMGSetRelChange(solver, 0);
        HYPRE_StructSMGSetNumPreRelax(solver, n_pre);
        HYPRE_StructSMGSetNumPostRelax(solver, n_post);
        // Logging must be on to get iterations and residual norm info below
        HYPRE_StructSMGSetLogging(solver, 1);

        // Setup and print setup timings
        HYPRE_StructSMGSetup(solver, A, b, x);
        hypre_EndTiming(time_index);
//#ifdef DEBUG
        hypre_PrintTiming("Setup phase times", MPI_COMM_WORLD);
//#endif
        hypre_FinalizeTiming(time_index);
        hypre_ClearTiming();

        // Solve and print solve timings
        time_index = hypre_InitializeTiming("SMG Solve");
        hypre_BeginTiming(time_index);
        HYPRE_StructSMGSolve(solver, A, b, x);
        hypre_EndTiming(time_index);
//#ifdef DEBUG
        hypre_PrintTiming("Solve phase times", MPI_COMM_WORLD);
//#endif
        hypre_FinalizeTiming(time_index);
        hypre_ClearTiming();

        // Get some info on the run
        HYPRE_StructSMGGetNumIterations(solver, &num_iterations);
        HYPRE_StructSMGGetFinalRelativeResidualNorm(solver, &final_res_norm);
//#ifdef DEBUG
        if (parms.rank == 0){
            fprintf(stdout, "Number of Iterations = %4d ; Final Relative Residual Norm = %e\n\n", num_iterations, final_res_norm);
        }
//#endif
        // Clean up 
        HYPRE_StructSMGDestroy(solver);
    }

    // 6. Set up and use a solver (PCG) with SMG Preconditioner
    if (solver_id == 1)
    {
        HYPRE_StructPCGCreate(MPI_COMM_WORLD, &solver);
        //HYPRE_StructPCGSetMemoryUse(solver, 0);
        HYPRE_StructPCGSetMaxIter(solver, 100);
        HYPRE_StructPCGSetTol(solver, 1.0e-12);
        HYPRE_StructPCGSetTwoNorm(solver, 1);
        HYPRE_StructPCGSetRelChange(solver, 0);
        //HYPRE_StructPCGSetPrintLevel(solver, 2 ); /* print each CG iteration */
        HYPRE_StructPCGSetLogging(solver, 1);
       
        /* Use symmetric SMG as preconditioner */
        HYPRE_StructSMGCreate(MPI_COMM_WORLD, &precond);
        HYPRE_StructSMGSetMemoryUse(precond, 0);
        HYPRE_StructSMGSetMaxIter(precond, 32);
        HYPRE_StructSMGSetTol(precond, 0.0);
        HYPRE_StructSMGSetZeroGuess(precond);
        HYPRE_StructSMGSetNumPreRelax(precond, 1);
        HYPRE_StructSMGSetNumPostRelax(precond, 1);
 
        /* Set the preconditioner and solve */
        HYPRE_StructPCGSetPrecond(solver, HYPRE_StructSMGSolve, HYPRE_StructSMGSetup, precond);
        HYPRE_StructPCGSetup(solver, A, b, x);
        HYPRE_StructPCGSolve(solver, A, b, x);
 
        /* Get some info on the run */
        HYPRE_StructPCGGetNumIterations(solver, &num_iterations);
        HYPRE_StructPCGGetFinalRelativeResidualNorm(solver, &final_res_norm);
#ifdef DEBUG
        if (parms.rank == 0){
            fprintf(stdout, "Number of Iterations = %4d ; Final Relative Residual Norm = %e\n\n", num_iterations, final_res_norm);
        }
#endif

        /* Clean up */
        HYPRE_StructSMGDestroy(precond);
        HYPRE_StructPCGDestroy(solver);
    }

    // get the local solution
    double *values = calloc(parms.N, sizeof(double));
    HYPRE_StructVectorGetBoxValues(x, ilower, iupper, values);

    // Free memory
    HYPRE_StructGridDestroy(grid);
    HYPRE_StructStencilDestroy(stencil);
    HYPRE_StructMatrixDestroy(A);
    HYPRE_StructVectorDestroy(b);
    HYPRE_StructVectorDestroy(x);
    free(Ab);
    return(values);
}
