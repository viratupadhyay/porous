// dissol.c

#define DEBUG

#undef	EULER		// Euler time update
#define	MIDPT		// Midpoint time update
#undef	RK4			// RK4 time update

#define PHI		1	// Output frequency - 0 = none
#define P		0
#define V		0
#define C		1
#define CHEK	5	// Checkpoint frequency - 0 = none

#define Dim     2

#define Nmob	1	// Number of mobile reactant species in Complex kinetics
#define Nimob	0	// Number of immobile species

#define SH      8.0 // Sherwood number

// Kinetics
#undef QCON
#undef RLIM		// Reaction limited kinetics
#undef TLIM			// Transport limited kinetics

#define SLN   0       // 0 - SMG, 1 - PCG
