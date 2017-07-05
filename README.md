[![PyPI](https://img.shields.io/pypi/v/phenum.svg)](https://pypi.python.org/pypi/phenum/) [![Build Status](https://travis-ci.org/wsmorgan/phonon-enumeration.svg?branch=master)](https://travis-ci.org/wsmorgan/phonon-enumeration)[![codecov](https://codecov.io/gh/wsmorgan/phonon-enumeration/branch/master/graph/badge.svg)](https://codecov.io/gh/wsmorgan/phonon-enumeration)[![Code Issues](https://www.quantifiedcode.com/api/v1/project/675d6268247c4c2cb80669f832bef46c/badge.svg)](https://www.quantifiedcode.com/app/project/675d6268247c4c2cb80669f832bef46c)

# phonon-enumeration

This code is used to enumerate all the derivative structures of a
system within a crystalographic system within specified concentration
and cell size ranges. The code uses a newly developed group theoretic
approach that is extremely efficient and can include the enumeration
of displacement directions, or arrow directions, within a system.

Full API Documentation available at: [github pages](https://wsmorgan.github.io/phenum/).

## Installing the code

To install the code use the following command in the
phonon-enumeration directory:

```
pip intsall phenum
```

Alternatively you can clone this repository and use:

```
python setup.py install
```

from within the phonon-enumeration directory.

## Running the code

### Complete enumeration of a system.

If you want to enumerate every possible structure specified in your
'struct_enum.in' file, a sample of which can be found in the input
folder, then use the commend:

```
enumeration-py -enum
```

### Enumerating a subset

If the total number of arrangements in your system is huge you may not
want to enumerate them all. In that case it is possible to enumerate a
subset of the system as follows. Examples of all input files can be
found in the input folder. First have the code find the polya
distribution for the system described in your lattice.in file:

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
