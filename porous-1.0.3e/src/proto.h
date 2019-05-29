/* Structures */
struct parms{int Nx, Ny, Nz, nsten, Nvar, N, cyc; 
		int nx, ny, Np, Nq, rank, size, p, q;
		double Re, Pe, Da, DC, k, pow, t, dt, dx, avgV, V0, havg;
        int x0, y0, x1, y1;
		int L, R, U, D;};	// Rank of neighbors

/* I/O */
struct parms readb(char *, double *, struct parms);
struct parms readf(char *, double *, struct parms);
void filename(char *, char *, char *, struct parms);
void nameinput(char *, char *, char *);
void writef(char *, double *, int, struct parms, int);
void writeb(char *, double *, int, struct parms);

/* Functions */
struct parms mapping(struct parms);
int *spiral(struct parms);
double *padding(double *, int *, int, struct parms);
double *matP(double *, int *, struct parms);
double *matC(double *, int *, double *, double *, struct parms);
void writedata(double *, double *, int, int, struct parms);
double sumf(double *, int, struct parms);
double avgf(double *, int, struct parms);
struct parms fluxV(double *, int *, struct parms);
void fluxC(double *, double *, double *, struct parms);
void erode(double *, double *, double *, int *, struct parms);

/* Models */
double K(double phi);
double S(double phi);
double *D(double *, int, struct parms);

/* Hypre */
double *solve(double *, int, struct parms);

/* Utilities */
double wtime(struct timeval *);
void winit(struct timeval *);
void warning(char *);
void error(char *);
