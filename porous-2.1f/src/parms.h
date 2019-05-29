#undef DEBUG

#undef	EULER		// Euler time update
#define	MIDPT		// Midpoint time update
#undef	RK4			// RK4 time update

#define PHI		1   // Output frequency - 0 = none
#define P		0
#define V		0
#define C		0
#define CHEK	5	// Checkpoint frequency - 0 = none

#define Dim    2    // Dimensions

#define Nmob	1	// Number of mobile reactant species in Complex kinetics
#define Nimob	0	// Number of immobile species

#undef APLIM        // Aperture limits

// Constants
#define SH      8.0 // Sherwood number

// Kinetics
#define AMAX    2.0 // Maximum Aperture
#undef QCON
#undef RLIM		    // Reaction limited kinetics
#undef TLIM			// Transport limited kinetics
#undef MIXED        // Mixed kinetics

// Solver
#define SLN     0   // 0 - SMG, 1 - PCG

// Constitutive models
#undef SAF     0.667   // Surface area model factor
