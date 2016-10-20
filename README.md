[![PyPI](https://img.shields.io/pypi/v/phenum.svg)](https://pypi.python.org/pypi/phenum/) [![Build Status](https://travis-ci.org/wsmorgan/phonon-enumeration.svg?branch=master)](https://travis-ci.org/wsmorgan/phonon-enumeration)[![Coverage Status](https://coveralls.io/repos/github/wsmorgan/phonon-enumeration/badge.svg?branch=phenum-dev)](https://coveralls.io/github/wsmorgan/phonon-enumeration?branch=phenum-dev)

# phonon-enumeration

This code is used to enumerate all the derivative structures of a
system within a crystalographic system within specified concentration
and cell size ranges. The code uses a newly developed group theoretic
approach that is extremely efficient and can include the enumeration
of displacement directions, or arrow directions, within a system.

Full API Documentation available at: [github pages](https://wsmorgan.github.io/phenum/).

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

# Enumerating a system

You now have everything you need to run the new enumeration code. The
work flow for this code is as follows. Samples of all input files can
be found in the input folder. First find the polya distribution for the
system described in your lattice.in file:

```
enumeration.py -polya
```

Next we need to build an enum.in file. You may either build this by
hand or have the code build it for you using the `-distribution`
option which takes two arguments, the type of distribution and the
number of structures we want in the results.

```
enumeration.py -distribution all all
enumeration.py -distribution all 100
```

If any option other than 'all' is passed into the first argument then
the code will not produce an enum.in file that will be useful for the
actual enumeration. The options of 'HNF', 'shape', and 'conc' are
simply for the user's viewing purposes.

Once an enum.in file has been constructed we can enumerate the entire
set of unique configurations:

```
enumeration.py -enum
```

This will make an enum.out file listing the unique configurations.

# Making POSCARS

Phenum contains a second executable for making POSCARS. To make a
POSCAR first select a structure number, or range of structures, for
the POSCARs to be constructed for from the enum.out file. Then run:

```
makeStr.py 10
```

This would make the POSCAR for the 10th structure. For a range of
structures use:

```
makeStr.py 20 30
```

To make POSCARs for the 20th to 30th structures. The POSCARS are saved
as vasp.* files. To have the code calculate the lattice parameter as
well use:

```
makeStr.py 10 -species Al Ni
```

Where the Al and Ni are replaced with the elements in the system being
modeled.

# Visualization

At times it is useful to construct the distribution based off the
shapes of the supercells. For the non-expert user an option has been
added to the code to make pictures of the supercells. To do this first
we need to make a distribution of only the supercells:

```
enumeration.py -distribution shape all -savedist
```

This created an enum.in file that lists only the supercells and the
number of unique arrangements within each supercell. We can now
visualize each of the supercells:

```
enumeration.py -visualize -shape
```

This creates a pdf file for each of the supercells. The `-shape`
option forces the code to include the lines that define the cell in
the pdfs.

## Python Packages Used

The enumeration.py code require the following python packages to run:

- numpy

- pyparsing

- termcolor

- matplotlib
