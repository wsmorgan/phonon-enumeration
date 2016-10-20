"""Methods for testing the subroutines in the enumeration module."""
import unittest as ut
import os

import pytest

def get_sargs(args):
    """Returns the list of arguments parsed from sys.argv.
    """
    import sys
    sys.argv = args
    from phenum.enumeration import _parser_options
    return _parser_options()    

def test_examples():
    """Makes sure the script examples work properly.
    """
    argv = ["py.test", "-examples"]
    assert get_sargs(argv) is None

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

    def _check_total(self,test_file,number):
        """Compares the number of structures in an enum.in file and compares
        it to the given number.
        :arg test_file: The test output file.
        :arg number: The correct number of structures in the file.
        """
        total = 0
        with open(test_file,"r") as of:
            for line in of:
                if "#" in line:
                    pass
                else:
                    total += int(line.strip().split()[-1])

        self.assertEquals(total,number)

    def test_EnumIn1(self):
        from phenum.enumeration import _enum_in
        from os import system
        args = {'profile': None, 'savedist': False, 'verbose': None, 'exec': 'enum.x',
                'outfile': 'test_enum.in', 'enum': False, 'cellsdir': 'tests/enumeration/fcc_1/',
                'lattice': 'lattice.in', 'dataformat': 'cells.{}', 'acceptrate': None,
                'examples': False, 'sizes': [2,6], 'debug': False, 'input': 'enum.in',
                'polya': False, 'super': False, 'distribution': ['all','all'],'seed':None,
                'filter':None}
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
                'polya': False, 'super': False, 'distribution': ['all','all'],'seed':None,
                'filter':None}
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
                'polya': False, 'super': False, 'distribution': ['all','100'],'seed':999,
                'filter':None}
        _enum_in(args)
        if sys.version_info[0] < 3:
            self._compare_files('test_enum.in','tests/enumeration/fcc_1/enum.in_100_p2')
        else:
            self._compare_files('test_enum.in','tests/enumeration/fcc_1/enum.in_100_p3')
        system("rm test_enum.in")

    def test_EnumIn4(self):
        from phenum.enumeration import _enum_in
        from os import system, getcwd, chdir
        direct = getcwd()
        import sys
        args = {'profile': None, 'savedist': False, 'verbose': None, 'exec': 'enum.x',
                'outfile': 'test_enum.in', 'enum': False, 'cellsdir': 'tests/enumeration/fcc_2/',
                'lattice': 'lattice.in', 'dataformat': 'cells.{}', 'acceptrate': None,
                'examples': False, 'sizes': None, 'debug': False, 'input': 'enum.in',
                'polya': False, 'super': False, 'distribution': ['all','50'],'seed':2461,
                'filter':None}
        _enum_in(args)
        if sys.version_info[0] < 3:
            self._compare_files('test_enum.in','tests/enumeration/fcc_2/enum.in_50_p2')
        else:
            self._compare_files('test_enum.in','tests/enumeration/fcc_2/enum.in_50_p3')
        system("rm test_enum.in")
        chdir(direct)

    def test_EnumIn5(self):
        from phenum.enumeration import _enum_in
        from os import system, getcwd, chdir
        direct = getcwd()
        args = {'profile': None, 'savedist': False, 'verbose': None, 'exec': 'enum.x',
                'outfile': 'test_enum.in', 'enum': False, 'cellsdir': 'tests/enumeration/sc_1/',
                'lattice': 'lattice.in', 'dataformat': 'cells.{}', 'acceptrate': None,
                'examples': False, 'sizes': None, 'debug': False, 'input': 'enum.in',
                'polya': False, 'super': False, 'distribution': ['all','100'],'seed':None,
                'filter':None}
        _enum_in(args)
        self._compare_files('test_enum.in','tests/enumeration/sc_1/enum.in_100')
        system("rm test_enum.in")
        chdir(direct)

    def test_EnumIn6(self):
        from phenum.enumeration import _enum_in
        from os import system
        args = {'profile': None, 'savedist': True, 'verbose': None, 'exec': 'enum.x',
                'outfile': 'test_enum.in', 'enum': False, 'cellsdir': 'tests/enumeration/sc_1/',
                'lattice': 'lattice.in', 'dataformat': 'cells.{}', 'acceptrate': None,
                'examples': False, 'sizes': None, 'debug': False, 'input': 'enum.in',
                'polya': False, 'super': False, 'distribution': ['shape','all'],'seed':None,
                'filter':None}
        _enum_in(args)
        self._compare_files('test_enum.in','tests/enumeration/sc_1/enum.in_shape')
        system("rm test_enum.in")

    def test_EnumIn7(self):
        from phenum.enumeration import _enum_in
        from os import system
        args = {'profile': None, 'savedist': True, 'verbose': None, 'exec': 'enum.x',
                'outfile': 'test_enum.in', 'enum': False, 'cellsdir': 'tests/enumeration/sc_1/',
                'lattice': 'lattice.in', 'dataformat': 'cells.{}', 'acceptrate': None,
                'examples': False, 'sizes': None, 'debug': False, 'input': 'enum.in',
                'polya': False, 'super': False, 'distribution': ['size','all'],'seed':None,
                'filter':None}
        _enum_in(args)
        self._compare_files('test_enum.in','tests/enumeration/sc_1/enum.in_size')
        system("rm test_enum.in")

    def test_EnumIn8(self):
        from phenum.enumeration import _enum_in
        from os import system
        args = {'profile': None, 'savedist': True, 'verbose': None, 'exec': 'enum.x',
                'outfile': 'test_enum.in', 'enum': False, 'cellsdir': 'tests/enumeration/sc_1/',
                'lattice': 'lattice.in', 'dataformat': 'cells.{}', 'acceptrate': None,
                'examples': False, 'sizes': None, 'debug': False, 'input': 'enum.in',
                'polya': False, 'super': False, 'distribution': ['conc','all'],'seed':None,
                'filter':None}
        _enum_in(args)
        self._compare_files('test_enum.in','tests/enumeration/sc_1/enum.in_conc')
        system("rm test_enum.in")

    def test_EnumIn9(self):
        from phenum.enumeration import _enum_in
        from os import system
        args = {'profile': None, 'savedist': False, 'verbose': None, 'exec': 'enum.x',
                'outfile': 'test_enum.in', 'enum': False, 'cellsdir': 'tests/enumeration/sc_1/',
                'lattice': 'lattice.in', 'dataformat': 'cells.{}', 'acceptrate': None,
                'examples': False, 'sizes': None, 'debug': False, 'input': 'enum.in',
                'polya': False, 'super': False, 'distribution': ['all','all'],'seed':None,
                'filter':['shape','enum.in_shape']}
        _enum_in(args)
        self._compare_files('test_enum.in','tests/enumeration/sc_1/enum.in_100')
        system("rm test_enum.in")

    def test_EnumIn10(self):
        from phenum.enumeration import _enum_in
        from os import system
        args = {'profile': None, 'savedist': False, 'verbose': None, 'exec': 'enum.x',
                'outfile': 'test_enum.in', 'enum': False, 'cellsdir': 'tests/enumeration/sc_1/',
                'lattice': 'lattice.in', 'dataformat': 'cells.{}', 'acceptrate': None,
                'examples': False, 'sizes': None, 'debug': False, 'input': 'enum.in',
                'polya': False, 'super': False, 'distribution': ['all','all'],'seed':None,
                'filter':['conc','enum.in_conc']}
        _enum_in(args)
        self._compare_files('test_enum.in','tests/enumeration/sc_1/enum.in_100')
        system("rm test_enum.in")

    def test_EnumIn11(self):
        from phenum.enumeration import _enum_in
        from os import system
        import sys
        args = {'profile': None, 'savedist': False, 'verbose': None, 'exec': 'enum.x',
                'outfile': 'test_enum.in', 'enum': False, 'cellsdir': 'tests/enumeration/sc_1/',
                'lattice': 'lattice.in', 'dataformat': 'cells.{}', 'acceptrate': None,
                'examples': False, 'sizes': None, 'debug': False, 'input': 'enum.in',
                'polya': False, 'super': False, 'distribution': ['all','10'],'seed':1010,
                'filter':['conc','enum.in_conc']}
        _enum_in(args)
        if sys.version_info[0] < 3:
            self._compare_files('test_enum.in','tests/enumeration/sc_1/enum.in_conc_10_p2')
        else:
            self._compare_files('test_enum.in','tests/enumeration/sc_1/enum.in_conc_10_p3')
        system("rm test_enum.in")

    def test_EnumIn12(self):
        from phenum.enumeration import _enum_in
        from os import system
        import sys
        args = {'profile': None, 'savedist': False, 'verbose': None, 'exec': 'enum.x',
                'outfile': 'test_enum.in', 'enum': False, 'cellsdir': 'tests/enumeration/sc_1/',
                'lattice': 'lattice.in', 'dataformat': 'cells.{}', 'acceptrate': None,
                'examples': False, 'sizes': None, 'debug': False, 'input': 'enum.in',
                'polya': False, 'super': False, 'distribution': ['all','10'],'seed':1010,
                'filter':['shape','enum.in_shape']}
        _enum_in(args)
        if sys.version_info[0] < 3:
            self._compare_files('test_enum.in','tests/enumeration/sc_1/enum.in_shape_10_p2')
        else:
            self._compare_files('test_enum.in','tests/enumeration/sc_1/enum.in_shape_10_p3')
        system("rm test_enum.in")

    def test_EnumIn13(self):
        from phenum.enumeration import _enum_in
        from os import system
        import sys
        args = {'profile': None, 'savedist': False, 'verbose': None, 'exec': 'enum.x',
                'outfile': 'test_enum.in', 'enum': False, 'cellsdir': 'tests/enumeration/sc_1/',
                'lattice': 'lattice.in', 'dataformat': 'cells.{}', 'acceptrate': None,
                'examples': False, 'sizes': [3,4], 'debug': False, 'input': 'enum.in',
                'polya': False, 'super': False, 'distribution': ['all','10'],'seed':1010,
                'filter':['conc','enum.in_conc']}
        _enum_in(args)
        if sys.version_info[0] < 3:
            self._compare_files('test_enum.in','tests/enumeration/sc_1/enum.in_conc_34_10_p2')
        else:
            self._compare_files('test_enum.in','tests/enumeration/sc_1/enum.in_conc_34_10_p3')
        system("rm test_enum.in")

    def test_EnumIn14(self):
        from phenum.enumeration import _enum_in
        from os import system
        import sys
        args = {'profile': None, 'savedist': False, 'verbose': None, 'exec': 'enum.x',
                'outfile': 'test_enum.in', 'enum': False, 'cellsdir': 'tests/enumeration/sc_1/',
                'lattice': 'lattice.in', 'dataformat': 'cells.{}', 'acceptrate': None,
                'examples': False, 'sizes': [3, 5], 'debug': False, 'input': 'enum.in',
                'polya': False, 'super': False, 'distribution': ['all','10'],'seed':1010,
                'filter':['shape','enum.in_shape2']}
        _enum_in(args)
        if sys.version_info[0] < 3:
            self._compare_files('test_enum.in','tests/enumeration/sc_1/enum.in_shape2_10_p2')
        else:
            self._compare_files('test_enum.in','tests/enumeration/sc_1/enum.in_shape2_10_p3')
        system("rm test_enum.in")

class TestPlotHNFs(ut.TestCase):
    """Tests of the _enum_in subroutine."""

    def test_visualize1(self):
        from phenum.enumeration import _script_enum
        from os import system
        args = {'profile': None, 'savedist': False, 'verbose': None, 'exec': 'enum.x',
                'outfile': None, 'enum': False, 'cellsdir': 'tests/enumeration/sc_1/',
                'lattice': 'input/fcc/lattice.in', 'dataformat': 'cells.{}', 'acceptrate': None,
                'examples': False, 'sizes': None, 'debug': False,
                'input': 'tests/enumeration/sc_1/enum.in_100',
                'polya': False, 'super': False, 'distribution': None,'seed':None,
                'filter':None,'visualize':True,'shapes':True,'show':False}
        _script_enum(args,testmode=True)
    
    def test_visualize2(self):
        from phenum.enumeration import _script_enum
        from os import system
        args = {'profile': None, 'savedist': False, 'verbose': None, 'exec': 'enum.x',
                'outfile': None, 'enum': False, 'cellsdir': 'tests/enumeration/sc_1/',
                'lattice': 'input/fcc/lattice.in', 'dataformat': 'cells.{}', 'acceptrate': None,
                'examples': False, 'sizes': None, 'debug': False,
                'input': 'tests/enumeration/sc_1/enum.in_100',
                'polya': False, 'super': False, 'distribution': None,'seed':None,
                'filter':None,'visualize':True,'shapes':False,'show':False}
        _script_enum(args,testmode=True)   
