#ver 3.0
-----------

Attempted changes

1) Try to change the perm model

2) Try to modeify the Surface area model

3) Try to implement the Dispersion model from Golfier(2002b)

4) Try to apply non-linear kinetics from Dreybrodt(1990)

#ver 2.1
-----------

1) Change in mapping function
   parms.D and parms.L calculated using modulus instead of  if loops

2) Fixed the bug in v2.0.1 for 2D codes. 
   In hypresolve added #if for ilower, iupper and pset 

3) AMAX instead of MAXAP

4) Surface area model modified in models.c

5) Mixed kinetics in flux.c

6) Fixed a bug in averaging velocity in z-direction in flux.c

7) Fixed a bug in function.c line 269 and 270 in matC. Changed to np instead of n.

#ver 2.0.1
-----------

1) The code now a 3D version. Works for both 2D and 3D.

2) Parameter in parms.h as Dim 2/3. 

3) Changes made in all functions accordingly

4) The code now has MAXAP and APLIM for the 3D case since its a porous rock.

5) Surface area model is turned ON too.

6) Rate constant is now a pointer of parms.N values instead of single constant

#ver 1.0.3
-----------

1) Put DEBUG #ifdef anywhere I need to debug e.g. No. of iterations in hypre 

2) checkpoint for aperture/porosity CHEK in parms.h ...write phi.chk in porous.c

3) read/write binary file in io.c

4) phi.dat.0000 ---> phi.dat.00 in io.c nameinit

5)******* hfile.txt individual for each simulation in init directory ********
    --- xconf xrun and xinit modified accordingly

6) function nameinit2 changed to nameinput and nameinit removed and merged with filename

7) changes in 'ida.py' and 'bplt.py' ... 
    --new function 'pget' to copy from clients by reading from 'hfile.txt' and then 'paradata' now calls 'pget' and reads those files and combines them to for single 'data' and in 'bplt.py' 
    --reading in each function 'plotc', 'plotx' and 'ploty' removed and read just once in 'Bplot' and makes any plot
    --directories 'imgs' and 'plots' copied to local dir after making all plots
=============================================================================
#ver 1.0.2

datafile name changed to type.cyc.rank e.g. c.0030.01

Bug in velocity fixed. 
-Averaging was being done only in local domains. But now adding all velocities
-Added lines of code in flux.c with an MPI command and functions.c
-Added #ifdef QCON whether to scale the velocity or not

=============================================================================
#ver 1.0.1

Now only 1 "input.dat" for all ranks. Rank removed from the name in function "nameinit2"

datafile name still type.rank.cyc (first rank then cyc) e.g. c.01.0030

Bug in velocity averaging
=============================================================================

#ver 1.0.0

"input.dat" has rank in its name e.g. input.dat.0000

Bug in velocity averaging
