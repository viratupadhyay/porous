/* Structures */
struct parms{int Nx, Ny, Nz, nx, ny, nsten, Nvar, N; 
		int Np, Nq, rank, size;
		int x0, y0, x1, y1;
		int L, R, U, D;};	// Rank of neighbors

/* I/O */
struct parms mapping(struct parms);
//struct parms readp(char *, struct parms);
struct parms readf(char *, double *, struct parms);
//struct parms readb(char *, double *, struct parms);
//void writef(char *, double *, struct parms, int);
//void writeb(char *, double *, struct parms); */
char *filename(char *, struct parms);

/* Functions */
int *spiral(struct parms);
double *padding(double *, int, struct parms);
double *matvec(double *, int, struct parms);
double perm(double, int);
//double coeff(struct parms);

/* Hypre */
void solve(double *, struct parms);
