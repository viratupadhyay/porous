# porous
A code in C runnable in parallel to simulate dissolution in thin fractures with a weak acid.
============================

13th February, 2013
v1.0.1:

============================
make files
.version with the version number e.g. 1.0.1

.workdir with the location of work directory e.g. /work/viratu/porous


-----------------
    x-scripts
-----------------

The sequence of x-scripts is same as dissol.

1) First creat initial condition using xconf

2) Then initalize it on the node by running xinit.

3) Finally xrun to start the simulations

-- xconf has 1 added parameter at the end i.e. (Number of processors+1)
e.g. xconf test 1 128 128 1 5    #(for 4 processors)

-- xrun has 1 parameter added i.e. Number of processors
e.g. xrun test 0 4               #(for 4 processors)

-- xinit and xclean are same too.
e.g. xinit test
     xclean test

============================

24th December, 2012
v1.0.0:

============================

-- First working parallel version of dissol named porous

-- Parameters and models based on dissol-1.3.2

-- models.c can be modified to add more models

============================

Working of the code

============================
-----------------
	input.dat 
-----------------

-- Nx, Ny, Nz

-- Np (# of processors for 2D domain decomposition)

-- Nvar (# of variables in data pointer)

-- Remaining parameters in input.dat are like dissol

-----------------
   fracs.dat
-----------------

-- non-linear kinetics lines removed

-----------------
    x-scripts
-----------------

-- xconf has 1 added parameter at the end Nz
e.g. xconf test 1 128 128 1

-- xrun is same as dissol
e.g. xrun test 0

-- xinit and xclean are same too.
e.g. xinit test
     xclean test

-----------------
  data types
-----------------

 phi = porosity
 p   = pressure
 vx  = velocity in x-direction
 vy  = velocity in y-direction
 c   = concentration

-----------------
   File name
-----------------

e.g. phi.00.0010

"00" is the rank of processor
"0010" is the cyc
