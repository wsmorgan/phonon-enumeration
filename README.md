# phonon-enumeration

This code is used to enumerate all the derivative structures of a
system within a crystalographic system within specified concentration
and cell size ranges. The code uses a newly developed group theoretic
approach that is extremely efficient and can include the enumeration
of displacement directions, or arrow directions, within a system.

## Prerequisites

The code currently requires a modified version of the previous
enumeration code, available at https://github.com/msg-byu/enumlib, to
run. To make this modified code first get a copy of the original code
then copy the files 'derivative_structure_generator.f90', wrapper.f90
and 'Makefile' from the support folder into the old codes
directory. Then you can follow the instruction to compile the code
found on that site.

Once the old code is compiled you will need to make a file called
struct_enum.in, an example of which can be found in the input folder,
for the system you desire to model. Then run the compiled enum.x from
the enumlib code. This will now generate a number of files titled
cell_# where # is the cell size. These files contain the information
needed to run the new enumeration code such as the symmetry group and
the HNFs that exist for that cell size.

## Running the code

You now have everything you need to run the new enumeration code. You
have two options for how to proceed. First the algorithm can use the
burnside polya algorithm to predict the number of unique structures
that exist for each HNF and symmetry group produced. This mode is run
as follows:

```
python phonon_enumeration.py -polya
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
python phonon_enumeration.py -enum
```

and expects a file called enum.in which can also be found in the input
folder. Depending on what is found in the input folder the code will
either do a full enumeration or enumerate a subset of the structures
and unique configurations. If a subset is desired then the last
entries in the input file should list the HNF, concentration, and
number of configurations desired of each system in that order on
seperate lines.

## Python Packages Used

The phonon_enumeration.py code require the following python packages to run:
-copy
-numpy
-random

