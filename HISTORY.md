# Revision History for "phonon-enumeration"

## Revision 2.1.3
- Fixed the hash function in `makeStr.py`.

## Revision 2.1.2
- Fixed bugs in `write_POSCAR` and `write_config` that were introduced
  by the merge and by the need to maintain element names in the first
  line of the file.

## Revision 2.1.1
- Fixed `makeStr.py` so that when making a config file output so that
  it always puts the atoms in order of number and not in order of
  species. If the atom numbers are out of order MTP can't read the
  config file.

- Fixed `makeStr.py` so that the default behavior is for it not to
  remove zeros from the concentration string. To remove zeros the
  `-remove_zero=t` flag must be used.

## Revision 2.1.0
- Changed makeStr.py so that it can optionally output on MTP config
  file instead of VASP POSCAR files. Also added a 'species_mapping'
  option to makeStr.py so that a binary enumeration can get used the
  edge of a ternary phase diagram for any edge of the ternary system
  by alternating the species_mapping.

## Revision 2.0.11
-Fixed a bug found by feifzhou in
 https://github.com/msg-byu/enumlib/issues/54 in which the code tries
 to read in the structures before `system["nD"]` has been set.

## Revision 2.0.10
-Fixed the calculation of the lattice parameter so that the lattice
 parameter for each atom will be determined by it's atomic volume with
 respect to the volume per atom of the parent cell.

## Revision 2.0.9
- Reversed the order of the tite written to the VASP style POSCARS by
  makeStr.py so that the atomic species get written first then the
  system title as given in enum.out. This is so that the POSCARS can
  be easily be read in by the quippy Atoms object.

## Revision 2.0.8
- Added the option to create subsets straight from the command line
  using `-enum -distribution all #`. This allows for the enumeration
  of subsets without having to run polya and distribution seperately
  first.

## Revision 2.0.7
- Fixed the bug reported in [issue
  #62](https://github.com/wsmorgan/phonon-enumeration/issues/62).

## Revision 2.0.6
- Fixed the issue reported in [issue
  #61](https://github.com/wsmorgan/phonon-enumeration/issues/61).
- Added a check to SmithNormalForm to ensure that the input matrix is
  integer.

## Revision 2.0.5
- Fixed the issue reported in [issue
  #58](https://github.com/wsmorgan/phonon-enumeration/issues/58).
- Removed a lot of unused variables.

## Revision 2.0.4
- Fixed the issue reported in issue #56 by changing M=list(HNF) to
  M=deepcopy(list(HNF)) in SmithNormalForm.
- Finished fixing the documentation.

## Revision 2.0.3
- Fixed final 4 issues found by quantifiedcode, this involved minor
  changes to phenum/symmetry.py, phenum/grouptheory.py,
  phenum/vector_utils.py, and phenum/HNFs.py.
- Added CONTRIBUTING.md to repo.

## Revision 2.0.2
- Removed deprecated code from phenum/structures.py.
- Fixed lots of minor style issues found by quantified code.
- Finished improving sphinx documentation.
- Moved phenum.symmetry._get_lattice_pointGroup to phenum.symmetry.get_lattice_pointGroup.

## Revision 2.0.1
- Started changing numpy style documentation to google style to be consitent.
- Switched to codecov and quantified code.
- Changed the nested for loops in io_utils.py and symmetry.py to itertools iterators.

## Revision 2.0.0
- Added HNFs.py to the repo which finds the unique HNFs.
- Updated the driver so that the polya mode is now stand alone.
- Updated the driver so that the distribution mode is now stand alone.
- Updated the driver so that the enumeration mode is now stand alone.
- Implemented new unit tests.
- Fixed bug in makeStr.py, the lattice vectors weren't being
  transposed after being read in.
- Fixed bugs revealed by study if hexagonal crystals.
- Fixed some integer issues for python 3.
- Finished unit testing. Now at 100% test coverage.

## Revision 1.8.3
- Removed an error message from vector_utils.py that was causing the
  code to die in cases where the code was behaving as it should.
- Fixed the data being written out to file so that the number of unit
  cells rather than the total number of atoms are written out for the
  idx value in 'enum.out'.
- Changed the largest allowed Coefficient of the
  guess_and_check_brancher subroutine. The previous limit was causing
  -1's to show up in the labeling.

## Revision 1.8.2
- Fixed the bug in tree.py that was causing the enumretation to hang. See Issue #48.

## Revision 1.8.1
- Forgot to import the new guess_and_check_brancher into the enum_sys
  subgroutine. Fixed now.

## Revision 1.8.0
- Updated the unit test outputs and test cases to reflect the changes
  in the code introduced by including the ability to eliminate
  superperiodic structures.
- Added the ability to include/exclude superperiodic structures from
  the list of unique configurations.
- Updated the fortran makefile to reflect the recent changes in enumlib.
- Changed the code so that when we are looking for a very small number
  of configurations within a very large number of possibilities the
  code will no longer do the full tree search but will instead do a
  'guess and check' method that picks a new configuration then checks
  to make sure it's not a duplicate of those present in the survivors
  list.

## Revision 1.7.1
- Refactored the drivers (enumeration.py and makeStr.py) for ease of unit testing.
- Changed some error messages to raise value erros instead.
- Fixed a number of bugs in grouptheory.py related to improper
  translation of original fortran codes matrix operations.
- Removed depricated binomial calculator from numerics.py.
- Cleaned up the calls to numpy and removed unneeded calls to
  numpy.array.
-Fixed an indexing error in polyaburnside.py.

## Revision 1.7.0

- Added the displaying of HNFs for the non-expert user as discussed in #22.
- Added the rattle option as described in #33.
- Fixed the bug reported in #41.
- Fixed the sheild in the README.md.
- Added lots of unit tests.

## Revision 1.6.0

- Added an option to save the distribution as described in #20.
- Fixed the bug reported in #14.
- Added the filter option for finer control over the enumerations as described in #21.
- Added unit tests for structures.py and the -distribution option of enumeration.py.

## Revision 1.5.10

- Fixed the bugs reported in #35, #36, and #38.
- Added the -cellsdir option as described in #31 in order to allow the
  user to specify where the cells.{} folders are


## Revision 1.5.9

- Fixed the bug reported in #34.
- Updated the -help info for the -displacement flag.

## Revision 1.5.8

- Added options to makeStr.py to allow the lattice parameter to be
  calculated and override the default title as descrbide in #32.

## Revision 1.5.7

- Fixed the bug reported in issue #15.
- Added a different workflow for the purely arrow enumerations where
  only a small number of unique arrangements are wanted.
  
## Revision 1.5.6

- Optimization for the arrow enumeration.
- Added progress bars for fun (they are useful...)
- Added profiling with `vprof`. 

## Revision 1.5.5

- Added support for `python 3`.

## Revision 1.5.4

- Added the actualy displacement to the atomic positions in the output POSCAR.
- Fixed a rounding error found in vector_utils.py that was causing the
  minkouski_reduction to fail.

## Revision 1.5.3

- Fixed the single structure inputs which broke when introducing the
  range option.
- Added the all option to the structures parameter as suggested in #28.

## Rivision 1.5.2

- Changed the input argumants as mentioned in #28

## Revision 1.5.1

- Had to change revision number to push up a new pyPI package because
  I maed an error.

## Revision 1.5.0

- Added makestr.py to the repo as suggested in #7.
- Added vector_utils.py to enable makestr.py.

## Revision 1.4.4

- Fixed the output file error reported in issue # 23.
- Fixed the output file error reported in issue # 24.
- Added an error message to indicate of the code ever brakes in a major way.

## Revision 1.4.3

- Added `-sizes` option to `enumeration.py` so that distribution can be size limited.
- Added sorting by cell size, then value, then HNF to the distribution print out.

## Revision 1.4.2

- Added counting by cell size to the enumeration output.
- Added sorting of output by `labeling` in the `enum.out`file.

## Revision 1.4.1

- Minor fixes for the pip install.


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
