"""Methods for testing the subroutines in the enumeration module."""
import unittest as ut

class TestEnumIn(ut.TestCase):
    """Tests of the _enum_in subroutine."""

    def _compare_files(self,file1,file2):
        """Compares the contents of two enum.in type files.
        :arg file1: The test generated output file.
        :arg file2: The control output file.
        """
        out1 = []
        with open(file1,"r") as f1:
            for line in f1:
                out1.append(line.strip())
        out2 = []
        with open(file2,"r") as f2:
            for line in f2:
                out2.append(line.strip())
        test = len(out1) == len(out2) and sorted(out1) == sorted(out2)
        self.assertEqual(True,test)

    def test_EnumIn1(self):
        from phenum.enumeration import _enum_in
        from os import system
        args = {'profile': None, 'savedist': False, 'verbose': None, 'exec': 'enum.x',
                'outfile': 'test_enum.in', 'enum': False, 'cellsdir': 'tests/enumeration/fcc_1/',
                'lattice': 'lattice.in', 'dataformat': 'cells.{}', 'acceptrate': None,
                'examples': False, 'sizes': [2,6], 'debug': False, 'input': 'enum.in',
                'polya': True, 'super': False, 'distribution': ['all','all'],'seed':None}
        _enum_in(args)
        self._compare_files('test_enum.in','tests/enumeration/fcc_1/enum.in_2_6')
        system("rm test_enum.in")

    def test_EnumIn2(self):
        from phenum.enumeration import _enum_in
        from os import system
        args = {'profile': None, 'savedist': False, 'verbose': None, 'exec': 'enum.x',
                'outfile': 'test_enum.in', 'enum': False, 'cellsdir': 'tests/enumeration/fcc_2/',
                'lattice': 'lattice.in', 'dataformat': 'cells.{}', 'acceptrate': None,
                'examples': False, 'sizes': [3,4], 'debug': False, 'input': 'enum.in',
                'polya': True, 'super': False, 'distribution': ['all','all'],'seed':None}
        _enum_in(args)
        self._compare_files('test_enum.in','tests/enumeration/fcc_2/enum.in_3_4')
        system("rm test_enum.in")

    def test_EnumIn3(self):
        from phenum.enumeration import _enum_in
        from os import system
        import sys
        args = {'profile': None, 'savedist': False, 'verbose': None, 'exec': 'enum.x',
                'outfile': 'test_enum.in', 'enum': False, 'cellsdir': 'tests/enumeration/fcc_1/',
                'lattice': 'lattice.in', 'dataformat': 'cells.{}', 'acceptrate': None,
                'examples': False, 'sizes': None, 'debug': False, 'input': 'enum.in',
                'polya': True, 'super': False, 'distribution': ['all','100'],'seed':999}
        _enum_in(args)
        if sys.version_info[0] < 3:
            self._compare_files('test_enum.in','tests/enumeration/fcc_1/enum.in_100_p2')
        else:
            self._compare_files('test_enum.in','tests/enumeration/fcc_1/enum.in_100_p3')
        system("rm test_enum.in")

    def test_EnumIn4(self):
        from phenum.enumeration import _enum_in
        from os import system
        import sys
        args = {'profile': None, 'savedist': False, 'verbose': None, 'exec': 'enum.x',
                'outfile': 'test_enum.in', 'enum': False, 'cellsdir': 'tests/enumeration/fcc_2/',
                'lattice': 'lattice.in', 'dataformat': 'cells.{}', 'acceptrate': None,
                'examples': False, 'sizes': None, 'debug': False, 'input': 'enum.in',
                'polya': True, 'super': False, 'distribution': ['all','50'],'seed':2461}
        _enum_in(args)
        if sys.version_info[0] < 3:
            self._compare_files('test_enum.in','tests/enumeration/fcc_2/enum.in_50_p2')
        else:
            self._compare_files('test_enum.in','tests/enumeration/fcc_2/enum.in_50_p3')
        system("rm test_enum.in")

    def test_EnumIn5(self):
        from phenum.enumeration import _enum_in
        from os import system
        args = {'profile': None, 'savedist': False, 'verbose': None, 'exec': 'enum.x',
                'outfile': 'test_enum.in', 'enum': False, 'cellsdir': 'tests/enumeration/sc_1/',
                'lattice': 'lattice.in', 'dataformat': 'cells.{}', 'acceptrate': None,
                'examples': False, 'sizes': None, 'debug': False, 'input': 'enum.in',
                'polya': True, 'super': False, 'distribution': ['all','100'],'seed':None}
        _enum_in(args)
        self._compare_files('test_enum.in','tests/enumeration/sc_1/enum.in_100')
        system("rm test_enum.in")
