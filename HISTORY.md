# Revision History for "phonon-enumeration"

## Revision 1.4.0

- Refactored the code so that fewer subroutines are located in the driver.
- Fixed issue #12.
- Updated the code so that it doesn't need the struct_enum.template file anymore.
- Added the creation of the enum.in files as described in issue # 10.
- Added unit tests to all global subroutines excluding those in io_utils.py.
- Reformated enum.out slightly so that it still works with makestr.x from enumlib.
- Moved all the hash functions to phonons.py and tree.py and renamed radix.py to numerics.py since it only contains numeric subroutines now.
- Fixed lots of random bugs that occured during the refactoring.

## Revision 1.3.0

- Added the path to the enum.x code as a command line option as described in issue # 6.
- Updated derivative_structure_generator.f90 so that it no longer produces unwanted debugging files and other output.
- Updated the format of enum.out so that it properly displays the arrows now as described in issue $4.
- Fixed a bug in tree.py and in enumeration.py that occured when only some of a color were being displaced by an arrow.
- Also fixed some errors that occured when grouptheory.py was added to the code to find the arrow group.
- Fixed issue number 11.

## Revision 1.2.0

- Moved the enum folder to phenum as a better naming convention

- Added grouptheory.py, also moved contents of grouper.py into this
  file, and symmetry.py in order to generate the arrow group.

- Implemented the finding of the arrow group in the code in cases
  where arrow enumerations are wanted. Otherwise the code still reads
  in the symmetry group from file.

- Changed the arrow concentration input to be fractions that represent
  the percentage of each atomic species that will have arrows on their
  sites.

## Revision 1.1.1

- Fixed issue #1. As a result removed polya.py from repo as it is no
 longer needed.

- Fixed issue #2.

- Removed struct_enum.out from repo as it isn't needed.


## Revision 1.1.0

Refactored the code to bring it into best practices standards. All
tabs are now 4 spaces instead of 8. Other changes include:

- Moved binomial_coefficient.py and Inverse_radix_number.py into
radix_num_generator.py and renamed the file radix.py.

- Moved phonon_brancher.py to phonons.py.

- Moved arrow_group.py to grouper.py.

- Moved branch_method.py to tree.py.

- Removed the src/python folder and moved all of its contents into the
src folder.

- Added a .gitignore file.

- Removed the old phonon_out.txt output file that got added to the repo
by mistake.

If anything breaks then I made a mistake while refactoring the code.

## Revision 1.0.0

Changed the driver to make it easeir to run the code and added two
different modes for running, polya and enum. The first, polya, is used
to determine the distribution of unique arrangements for the different
cell sizes and concentrations. The second, enum, actually enumerates
the requested number of configurations for each of the desired cell
sizes and concentrations.

Added structures.py, io_utils.py, and polya.py to the src/python
folder. Also added the input folder that contains sample input and the
support folder that has the needed altered files for enumlib.

Also updated the documentation in all the code files so that they have
better formatting and are accurate to the current code.

## Revision 0.3.0

Changed from the python function to a hash function of my making so
that we could compare results easily with MatLab. Also added an
inverse hash function for the arrows. These two changes enabled a
major shift in how duplicates are found. Now instead of finding each
configuration by permuting the arrows we simply loop over the possible
hash numbers.

The output is now also saved to file for comparison to test cases.

I also fixed a bug that occured because of the change to hash
functions where the arrows weren't being permuted properly.

## Revision 0.2.1

Changed how the code determines duplicate arrows. It now does it using
a python function int(str(),n) where string is a string of the arrows
and n is the number of cardinal directions used. This number is then
compared between permuted arrows, if the permuted arrows have a
smaller numeral then the current string is not unique and won't be
saved.

The code is hard-coded for 4 cardinal directions at the moment but will
soon be switched to handle 2, 4, or 6 directions for more general use.

## Revision 0.1.1

Found two bugs in the code. The first was, if the code was handed only
a single color in the col list, i.e. 4 reds, all with or without
arrows it would fail, because, the brancher_method.py wasn't set up to
handle a single color. The error was fixed by having the code check to
see how many colors there are in the list at the start. If there is
only one then it goes to the code that adds arrows immediately instead
of starting the branch_method code. The arrow code will return the
original configuration if there are no arrows on the lattice and the
correct arrow configurations if there are.

The other error was that if two atoms of the same color were both
displaced with an arrow then the code would alter both arrows at the
same time instead of incrementing only one at a time. For example
[[1,0],[1,0]] would become [[1,1],[1,1]] and not [[1,0],[1,1]] as was
intended. The result is that the arrows always had the same
orientation. The bug was caused by pythons coping linking the two
colors/arrows together so that whatever happened to one would happen
to the other. It was fixed by making each implementation a 'deepcopy'
instead.

## Revision 0.1.0

Changed the file format for a more efficient algorithm. The code is
now run from phnon_enumeration.py instead of phonon_brancher.py. This
allows us to handle the arrows as each unique color configuration is
found, which enables us to use the stabilizers for the last color
configuration for the phonon excitations. This saves a great deal of
computation time.

Also removed the generator.py program from the repository since it has
been replaced by arrow_group.py code.

I also cleaned up a lot of the code to remove unused subroutines that
were creating a lot of clutter so the code is a lot cleaner. The
code's documentation has also been improved.

Lastly I merger radix_num_generator.py and radix_num_generator2.py so
that the code only depends on radix_num_generator.py and removed
radix_num_generator2.py from the repository.

## Initial Repository: Revision 0.0.0

First commits in initial development. Code currently works for the 2D
case when the symmetry group and the concentrations are entered into
the code as described in the README.
