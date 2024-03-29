# Compressible Falkner-Skan-Cooke solver

Solves the Falkner-Skan-Cooke equations by computing the compressible similarity solutions for a swept-wing  boundary layer at Pr=1 using the Stewartson transformation.

This is updated to build and run with gfortran on Darwin and Linux. 

## Building

to run a test problem simply build using

    ln -s gcc.mak Makefile
    make all

## Running

To run FSC enter the following:

    ./fsc < test.inp

or for the attachment line problem

    ./attach < at.inp

The Blasius Tollmein-Schlichting case from Collis Ph.D. thesis
Chapter 4 is run using:

    ./fsc < blasius.inp

See `stab/thesis` for more details including linear stability analysis for this flow. 

## Notes:

1. The gcc.mak uses gfortran and should be rather
   portable.  There are other, old, `*.mak` that you
   may find helpful but they have not been updated.
2. One may, optionally, use Numerical Recipes integration
   if you have the code (not supplied) and are 
   licensed to use it. To do so, build by using:

   `make USE_NR=1`

   this assumes that you have placed the needed NR
   routines in `nr_runge.f90` (again not supplied).
3. Note that you do not need NR software as there
   is not equivalent public domain software that
   we now use.
4. The primary output files are `cprofile.dat` and
   `sprofile.dat` these are to be compared with 
   `cprofile.ref` and `sprofile.ref` after running
   
   `./fsc < test.inp`

Also contains publically available software PPPack
from `NetLib.org`.

## Additional Information

Here are some typical parameters for this code.

### Case 1:
    M_inf = 0.1, R=80000
    lambda = 45 deg
    BetaH = 1
    Tw = 1, T0 = 1, muw = 1
    T0/TN0 = 1.001
    f'', g'
    1.2332556981943E+00  5.7053820816947E-01`

### Case 2: (Blasius TS case Collis, Ch. 4, Ph.D. thesis)
    M_inf = 0.3, R=1000
    lambda = 0
    BetaH = 0
    Tw = 1, T0=1, muw =1
    f'', g'
    4.6959998836720E-01  4.6959998883669E-01`

Or equivalently:  $M_\infty = 0.3$, $R_\delta{_1} = 1000$, $\lambda = 0$, $\beta_h = 0$, $T_w = 1$, $T_0 = 1$, $\mu_w = 1$
with initial guesses $f'' = 4.6959998836720\times10^{-1}$ and $g' = 4.6959998883669\times10^{-1}$.

S. Scott Collis\
flow.physics.simulation@gmail.com
