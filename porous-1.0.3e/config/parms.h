/* General */
#define MAXITR		100	/* Max iterations for NR */
#define TOL		1e-6	/* Tolerance iterative solvers */
#define DEBUG		1	/* Debug level (0=none; 1=min; 2=med; 3=max) */

/* dissol.c */
#undef	VARk			/* Variable reaction rate - needs init/k.dat */
#undef	EULER			/* Euler time update */
#undef	MIDPT			/* Midpoint time update */
#define	RK4			/* RK4 time update */

#define PHI		1	/* Output frequency - 0 = none */
#define K		1
#define P		1
#define V		1
#define C		1
#define CHEK		10	/* Checkpoint frequency - 0 = none */

#define Nmob		1	/* Number of mobile reactant species in Complex kinetics */
#define Nimob		0	/* Number of immobile species */
