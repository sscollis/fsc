# Compressible Falkner-Skan-Cooke solver

This is updated to build and run with gfortran on Darwin

## Building

to run a test problem simply build using

    ln -s gcc.mak Makefile
    make all

## Running

To run FSC enter the following:

    ./fsc < test.inp

or for the attachment line problem

    ./attach < at.inp

Notes:

1. The gcc.mak uses gfortran and should be rather
   portable.  There are other, old, `*.mak` that you
   may find helpful but they have not been updated.
2. One optionally use Numerical Recipes integration
   if you have the code (not supplied) and are 
   licensed to use it. To do so, build by using:

   `make USE_NR=1`

   this assumes that you have placed the needed NR
   routines in `numrec.f90` (again not supplied).
3. The primary output files are `cprofile.dat` and
   `sprofile.dat`

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

### Case 2:
    M_inf = 0.3, R=1000
    lambda = 0
    BetaH = 0
    Tw = 1, T0=1, muw =1
    f'', g'
    4.6959998836720E-01  4.6959998883669E-01`

S. Scott Collis \
Thu Feb 20 07:25:02 MST 2020
