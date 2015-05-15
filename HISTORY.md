# Revision History for "phonon-enumeration"

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
