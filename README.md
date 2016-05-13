# phonon-enumeration

This code is used to enumerate all the derivative structures of a
system within a crystalographic system within specified concentration
and cell size ranges. The code uses a newly developed group theoretic
approach that is extremely efficient and can include the enumeration
of displacement directions, or arrow directions, within a system.

## Prerequisites

The code currently requires a modified version of the previous
enumeration code, available at https://github.com/msg-byu/enumlib, to
run. To make this modified code do the following, get the symlib library:

```
git clone https://github.com/msg-byu/symlib.git
cd symlib/src/
make F90=(your compiler, gfortran or ifort)
cd ../../
```

Then get a copy of enumlib:

```
git clone https://github.com/msg-byu/enumlib.git
```

Now copy the Makefile, derivative_structure_generator.f90, and
wrapper.f90 from the phonon-enumeration/support directory to the
enumlib/src/ directory. Now we can make the enum.x executable and
place it in our path:

```
cd enumlib/src/
make F90=(your compiler) enum.x
cp enum.x /bin/.
```

In order for enum.x to run you will need to have its input folder
struct_enum.in, an example of which can be found in the input folder,
for the system you desire to model. You may then choose to run enum.x
yourself to generate the needed input files by typing:

```
enum.x
```

This will now generate a number of files titled cell_# where # is the
cell size. These files contain the information needed to run the new
enumeration code. If you do not run enum.x the enumeration.py code
will execute it for you as long as its in your path. The input files
are setup so that each HNF with it's SNF and left transform (as
described in http://msg.byu.edu/papers/multi.pdf and
http://msg.byu.edu/papers/GLWHart_enumeration.pdf) are listed in a
file titeled matrices:

```
  #n	SNF		   HNF			          left transform
   1  1  1  4    1  0  1  0  0  4      1    0    0    0    1    0    0    0    1
   1  1  1  4    1  0  1  0  1  4      1    0    0    0    1    0    0   -1    1   
```

The first digit indicates which of the group.n files contains the
symmetry group for that system. As can be seen only the diagonals of
the SNF and lower traingular entries of the HNF should be included in
this file. The group.n files contain the permutations of the sites on
the lattice that constitute the symmtery group.

## Installing the code

To install the code use the following command in the
phonon-enumeration directory:

```
python setup.py install
```

## Running the code

You now have everything you need to run the new enumeration code. You
have two options for how to proceed. First the algorithm can use the
burnside polya algorithm to predict the number of unique structures
that exist for each HNF and symmetry group produced. This mode is run
as follows:

```
enumeration.py -polya
```

and expects a file called lattice.in an example of which can be found
in the input folder. This mode produces a file for each cell size that
lists the number of unique configurations for each HNF at every
possible concentration range for the cell size. This data can be very
useful when modeling large systems as it will allow the user to select
an appropriate distribution of structures to use given the number of
each type available.

The second option is the actual enumeration of derivative
structures. This mode is run using:

```
enumeration.py -enum
```

and expects a file called enum.in which can also be found in the input
folder. The enum.in folder should contain a list of the desired HNFs,
their concentration ranges, and the number of arrangements for the HNF
concetrtaion range pair the user would like. For example:
```
# HNF                           Conc.       Number
  1 0 1 0 2 11                  8 3         2
  1 0 1 3 4 8                   4 4         1
  1 0 1 1 4 11                  6 5         3
  1 0 1 0 0 10                  6 4         2
  1 0 1 1 5 10                  8 2         1
  1 0 1 1 2 10                  7 3         1
  1 0 1 0 3 11                  7 4         3
  1 0 1 0 2 9                   5 4         2
```

## Python Packages Used

The enumeration.py code require the following python packages to run:

- numpy

- pyparsing

- termcolor
