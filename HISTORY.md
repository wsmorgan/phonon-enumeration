# Revision History for "phonon-enumeration"

## Revision 0.1.1

Found two bugs in the code. One was that if the code was handed only a
single color in the col list, i.e. 4 reds, with or without arrows it
would fail because the brancher_method.py wasn't set up to handle that
a single color. The error was fixed by having the code check to see
how many colors it has at the start. If there is only one then it goes
to the code that adds arrows immediately insead of starting into the
branch_method code. The arrow code will return the original
configuration if there are no arrows on the lattice and the correct
arrow configurations of there are.

The other error was that if two atoms of the same color were both
displaced with an arrow then the code would alter both arrows at the
same time instead of incrementing only one at a time. For example
[[1,0],[1,0]] would next become [[1,1],[1,1]] and not [[1,0],[1,1]] as
was inteded. The result is that the arrows always had the same
orientation. The bug was caused by pythons copy linking the two
colors/arrows together so that whatever happened to one would happen
to the other. It was fixed by making each implementation a deepcopy
instead.


## Revision 0.1.0

Changed the file format for a more efficient algorithm. The code is
now run from phnon_enumeration.py instead of phonon_brancher.py. This
allows us to handle the arraws as each unique color configuration is
found, which enables us to use the stabalizers for the last color
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
