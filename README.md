# phonon-enumeration

This set of codes is used to determine the unique arrangements of
atoms that have been displaced from a lattice. It is still in
development and not ready for wide spread use. If you desire to use it
anyway then in order to have it find the arrangements of your system
you must modify the source code in phonon_enumeration.py before
running. At line 20 you will need to change the elements of trans to
be the translation group of the lattice you are using. At line 21 you
will need to change the elements of rots to be [[group operation on
lattice],[effect on arrows]]. At present the code is only 2D so the
effect on the arrows should be a list permutation of numbers 0 to 3
(where 0 is up, 1 in left, 2 is down, and 3 is right) that map the
arrows to their new orientation. At line 19 you will need to change
the col to match your system, this is a fixed concentration code so
the values you enter will be the values used the first number is the
color to be used, the second number indicates if it is displaced (-1
indicates no displacement and 0-3 indicates displacement in one of the
directions used for the arrows).

Good luck.

## Python Packages Used

The phonon_enumeration.py code require the following python packages to run:
-copy
-numpy

## Running the program

The program is run using:
```
python phonon_enumeration.py
```
