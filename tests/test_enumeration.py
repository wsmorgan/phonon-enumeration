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
        from os import system, chdir, getcwd
        curdir = getcwd()
        chdir('tests/enumeration/fcc_1/')
        args = {'profile': None, 'savedist': False, 'verbose': None,
                'outfile': 'test_enum.in', 'enum': False,'lattice': 'lattice.in',
                'acceptrate': None, 'examples': False, 'sizes': [2,6], 'debug': False,
                'input': 'polya.out', 'polya': False, 'super': False,
                'distribution': ['all','all'],'seed':None,'filter':None}
        _enum_in(args)
        chdir(curdir)
        self._compare_files('tests/enumeration/fcc_1/test_enum.in',
                            'tests/enumeration/fcc_1/enum.in_2_6')
        system("rm tests/enumeration/fcc_1/test_enum.in")

    def test_EnumIn2(self):
        from phenum.enumeration import _enum_in
        from os import system, chdir, getcwd
        curdir = getcwd()
        chdir('tests/enumeration/fcc_2/')
        args = {'profile': None, 'savedist': False, 'verbose': None,
                'outfile': 'test_enum.in', 'enum': False,
                'lattice': 'lattice.in', 'acceptrate': None,
                'examples': False, 'sizes': [3,4], 'debug': False, 'input': 'polya.out',
                'polya': False, 'super': False, 'distribution': ['all','all'],'seed':None,
                'filter':None}
        _enum_in(args)
        chdir(curdir)
        self._compare_files('tests/enumeration/fcc_2/test_enum.in',
                            'tests/enumeration/fcc_2/enum.in_3_4')
        system("rm tests/enumeration/fcc_2/test_enum.in")

    def test_EnumIn3(self):
        from phenum.enumeration import _enum_in
        import sys
        from os import system, chdir, getcwd
        curdir = getcwd()
        chdir('tests/enumeration/fcc_1/')
        args = {'profile': None, 'savedist': False, 'verbose': None,
                'outfile': 'test_enum.in', 'enum': False,
                'lattice': 'lattice.in', 'acceptrate': None,
                'examples': False, 'sizes': None, 'debug': False, 'input': 'polya.out',
                'polya': False, 'super': False, 'distribution': ['all','100'],'seed':999,
                'filter':None}
        _enum_in(args)
        chdir(curdir)
        if sys.version_info[0] < 3:
            self._compare_files('tests/enumeration/fcc_1/test_enum.in',
                                'tests/enumeration/fcc_1/enum.in_100_p2')
        else:
            self._compare_files('tests/enumeration/fcc_1/test_enum.in',
                                'tests/enumeration/fcc_1/enum.in_100_p3')
        system("rm tests/enumeration/fcc_1/test_enum.in")

    def test_EnumIn4(self):
        from phenum.enumeration import _enum_in
        from os import system, getcwd, chdir
        import sys
        curdir = getcwd()
        chdir('tests/enumeration/fcc_2/')
        args = {'profile': None, 'savedist': False, 'verbose': None,
                'outfile': 'test_enum.in', 'enum': False,
                'lattice': 'lattice.in', 'acceptrate': None,
                'examples': False, 'sizes': None, 'debug': False, 'input': 'polya.out',
                'polya': False, 'super': False, 'distribution': ['all','50'],'seed':2461,
                'filter':None}
        _enum_in(args)
        chdir(curdir)
        if sys.version_info[0] < 3:
            self._compare_files('tests/enumeration/fcc_2/test_enum.in',
                                'tests/enumeration/fcc_2/enum.in_50_p2')
        else:
            self._compare_files('tests/enumeration/fcc_2/test_enum.in',
                                'tests/enumeration/fcc_2/enum.in_50_p3')
        system("rm tests/enumeration/fcc_2/test_enum.in")

    def test_EnumIn5(self):
        from phenum.enumeration import _enum_in
        from os import system, getcwd, chdir
        curdir = getcwd()
        chdir('tests/enumeration/sc_1/')
        args = {'profile': None, 'savedist': False, 'verbose': None,
                'outfile': 'test_enum.in', 'enum': False,
                'lattice': 'lattice.in', 'acceptrate': None,
                'examples': False, 'sizes': None, 'debug': False, 'input': 'polya.out',
                'polya': False, 'super': False, 'distribution': ['all','100'],'seed':None,
                'filter':None}
        _enum_in(args)
        chdir(curdir)
        self._compare_files('tests/enumeration/sc_1/test_enum.in',
                            'tests/enumeration/sc_1/enum.in_100')
        system("rm tests/enumeration/sc_1/test_enum.in")

    def test_EnumIn6(self):
        from phenum.enumeration import _enum_in
        from os import system, getcwd, chdir
        curdir = getcwd()
        chdir('tests/enumeration/sc_1/')
        args = {'profile': None, 'savedist': True, 'verbose': None,
                'outfile': 'test_enum.in', 'enum': False,
                'lattice': 'lattice.in', 'acceptrate': None,
                'examples': False, 'sizes': None, 'debug': False, 'input': 'polya.out',
                'polya': False, 'super': False, 'distribution': ['shape','all'],'seed':None,
                'filter':None}
        _enum_in(args)
        chdir(curdir)
        self._compare_files('tests/enumeration/sc_1/test_enum.in',
                            'tests/enumeration/sc_1/enum.in_shape')
        system("rm tests/enumeration/sc_1/test_enum.in")

    def test_EnumIn7(self):
        from phenum.enumeration import _enum_in
        from os import system, getcwd, chdir
        curdir = getcwd()
        chdir('tests/enumeration/sc_1/')
        args = {'profile': None, 'savedist': True, 'verbose': None,
                'outfile': 'test_enum.in', 'enum': False,
                'lattice': 'lattice.in', 'acceptrate': None,
                'examples': False, 'sizes': None, 'debug': False, 'input': 'polya.out',
                'polya': False, 'super': False, 'distribution': ['size','all'],'seed':None,
                'filter':None}
        _enum_in(args)
        chdir(curdir)
        self._compare_files('tests/enumeration/sc_1/test_enum.in',
                            'tests/enumeration/sc_1/enum.in_size')
        system("rm tests/enumeration/sc_1/test_enum.in")

    def test_EnumIn8(self):
        from phenum.enumeration import _enum_in
        from os import system, getcwd, chdir
        curdir = getcwd()
        chdir('tests/enumeration/sc_1/')
        args = {'profile': None, 'savedist': True, 'verbose': None,
                'outfile': 'test_enum.in', 'enum': False,
                'lattice': 'lattice.in', 'acceptrate': None,
                'examples': False, 'sizes': None, 'debug': False, 'input': 'polya.out',
                'polya': False, 'super': False, 'distribution': ['conc','all'],'seed':None,
                'filter':None}
        _enum_in(args)
        chdir(curdir)
        self._compare_files('tests/enumeration/sc_1/test_enum.in',
                            'tests/enumeration/sc_1/enum.in_conc')
        system("rm tests/enumeration/sc_1/test_enum.in")

    def test_EnumIn9(self):
        from phenum.enumeration import _enum_in
        from os import system, getcwd, chdir
        curdir = getcwd()
        chdir('tests/enumeration/sc_1/')
        args = {'profile': None, 'savedist': False, 'verbose': None,
                'outfile': 'test_enum.in', 'enum': False,
                'lattice': 'lattice.in', 'acceptrate': None,
                'examples': False, 'sizes': None, 'debug': False, 'input': 'polya.out',
                'polya': False, 'super': False, 'distribution': ['all','all'],'seed':None,
                'filter':['shape','enum.in_shape']}
        _enum_in(args)
        chdir(curdir)
        self._compare_files('tests/enumeration/sc_1/test_enum.in',
                            'tests/enumeration/sc_1/enum.in_100')
        system("rm tests/enumeration/sc_1/test_enum.in")

    def test_EnumIn10(self):
        from phenum.enumeration import _enum_in
        from os import system, getcwd, chdir
        curdir = getcwd()
        chdir('tests/enumeration/sc_1/')
        args = {'profile': None, 'savedist': False, 'verbose': None,
                'outfile': 'test_enum.in', 'enum': False,
                'lattice': 'lattice.in', 'acceptrate': None,
                'examples': False, 'sizes': None, 'debug': False, 'input': 'polya.out',
                'polya': False, 'super': False, 'distribution': ['all','all'],'seed':None,
                'filter':['conc','enum.in_conc']}
        _enum_in(args)
        chdir(curdir)
        self._compare_files('tests/enumeration/sc_1/test_enum.in',
                            'tests/enumeration/sc_1/enum.in_100')
        system("rm tests/enumeration/sc_1/test_enum.in")

    def test_EnumIn11(self):
        from phenum.enumeration import _enum_in
        import sys
        from os import system, getcwd, chdir
        curdir = getcwd()
        chdir('tests/enumeration/sc_1/')
        args = {'profile': None, 'savedist': False, 'verbose': None,
                'outfile': 'test_enum.in', 'enum': False,
                'lattice': 'lattice.in', 'acceptrate': None,
                'examples': False, 'sizes': None, 'debug': False, 'input': 'polya.out',
                'polya': False, 'super': False, 'distribution': ['all','10'],'seed':1010,
                'filter':['conc','enum.in_conc']}
        _enum_in(args)
        chdir(curdir)
        if sys.version_info[0] < 3:
            self._compare_files('tests/enumeration/sc_1/test_enum.in',
                                'tests/enumeration/sc_1/enum.in_conc_10_p2')
        else:
            self._compare_files('tests/enumeration/sc_1/test_enum.in',
                                'tests/enumeration/sc_1/enum.in_conc_10_p3')
        system("rm tests/enumeration/sc_1/test_enum.in")

    def test_EnumIn12(self):
        from phenum.enumeration import _enum_in
        import sys
        from os import system, getcwd, chdir
        curdir = getcwd()
        chdir('tests/enumeration/sc_1/')
        args = {'profile': None, 'savedist': False, 'verbose': None,
                'outfile': 'test_enum.in', 'enum': False,
                'lattice': 'lattice.in', 'acceptrate': None,
                'examples': False, 'sizes': None, 'debug': False, 'input': 'polya.out',
                'polya': False, 'super': False, 'distribution': ['all','10'],'seed':1010,
                'filter':['shape','enum.in_shape']}
        _enum_in(args)
        chdir(curdir)
        if sys.version_info[0] < 3:
            self._compare_files('tests/enumeration/sc_1/test_enum.in',
                                'tests/enumeration/sc_1/enum.in_shape_10_p2')
        else:
            self._compare_files('tests/enumeration/sc_1/test_enum.in',
                                'tests/enumeration/sc_1/enum.in_shape_10_p3')
        system("rm tests/enumeration/sc_1/test_enum.in")

    def test_EnumIn13(self):
        from phenum.enumeration import _enum_in
        import sys
        from os import system, getcwd, chdir
        curdir = getcwd()
        chdir('tests/enumeration/sc_1/')
        args = {'profile': None, 'savedist': False, 'verbose': None,
                'outfile': 'test_enum.in', 'enum': False,
                'lattice': 'lattice.in', 'acceptrate': None,
                'examples': False, 'sizes': [3,4], 'debug': False, 'input': 'polya.out',
                'polya': False, 'super': False, 'distribution': ['all','10'],'seed':1010,
                'filter':['conc','enum.in_conc']}
        _enum_in(args)
        chdir(curdir)
        if sys.version_info[0] < 3:
            self._compare_files('tests/enumeration/sc_1/test_enum.in',
                                'tests/enumeration/sc_1/enum.in_conc_34_10_p2')
        else:
            self._compare_files('tests/enumeration/sc_1/test_enum.in',
                                'tests/enumeration/sc_1/enum.in_conc_34_10_p3')
        system("rm tests/enumeration/sc_1/test_enum.in")

    def test_EnumIn14(self):
        from phenum.enumeration import _enum_in
        import sys
        from os import system, getcwd, chdir
        curdir = getcwd()
        chdir('tests/enumeration/sc_1/')
        args = {'profile': None, 'savedist': False, 'verbose': None,
                'outfile': 'test_enum.in', 'enum': False,
                'lattice': 'lattice.in', 'acceptrate': None,
                'examples': False, 'sizes': [3, 5], 'debug': False, 'input': 'polya.out',
                'polya': False, 'super': False, 'distribution': ['all','10'],'seed':1010,
                'filter':['shape','enum.in_shape2']}
        _enum_in(args)
        chdir(curdir)
        if sys.version_info[0] < 3:
            self._compare_files('tests/enumeration/sc_1/test_enum.in',
                                'tests/enumeration/sc_1/enum.in_shape2_10_p2')
        else:
            self._compare_files('tests/enumeration/sc_1/test_enum.in',
                                'tests/enumeration/sc_1/enum.in_shape2_10_p3')
        system("rm tests/enumeration/sc_1/test_enum.in")

    def test_EnumIn15(self):
        from phenum.enumeration import _script_enum

        args = {'profile': None, 'savedist': False, 'verbose': None,
                'outfile': None, 'enum': False,
                'lattice': 'lattice.in', 'acceptrate': None,
                'examples': False, 'sizes': None, 'debug': False, 'input': None,
                'polya': False, 'super': False, 'distribution': ['all'],'seed':1010,
                'filter':None}

        with pytest.raises(ValueError):
            _script_enum(args)

    def test_EnumIn16(self):
        from phenum.enumeration import _script_enum
        import sys
        from os import system

        args = {'profile': None, 'savedist': False, 'verbose': None,
                'outfile': "enum.in", 'enum': False,
                'lattice': 'tests/enumeration/enum_out/lattice.in_hcp', 'acceptrate': None,
                'examples': False, 'sizes': None, 'debug': False,
                'input': None,
                'polya': False, 'super': False, 'distribution': ['all','all'],'seed':None,
                'filter':None,'visualize':False,'shapes':False,'show':False}

        _script_enum(args)
        if sys.version_info[0] < 3:
            self._compare_files("enum.in","tests/enumeration/enum_out/enum.in_hcp_p2")
        else:
            self._compare_files("enum.in","tests/enumeration/enum_out/enum.in_hcp_p3")
        system("rm enum.in polya.out.2 polya.out.1")
            
class TestPlotHNFs(ut.TestCase):
    """Tests of the _enum_in subroutine."""

    def test_visualize1(self):
        from phenum.enumeration import _script_enum
        from os import system
        args = {'profile': None, 'savedist': False, 'verbose': None,
                'outfile': None, 'enum': False,
                'lattice': 'input/fcc/lattice.in', 'acceptrate': None,
                'examples': False, 'sizes': None, 'debug': False,
                'input': 'tests/enumeration/sc_1/enum.in_100',
                'polya': False, 'super': False, 'distribution': None,'seed':None,
                'filter':None,'visualize':True,'shapes':True,'show':False}
        _script_enum(args,testmode=True)
    
    def test_visualize2(self):
        from phenum.enumeration import _script_enum
        from os import system
        args = {'profile': None, 'savedist': False, 'verbose': None,
                'outfile': None, 'enum': False,
                'lattice': 'input/fcc/lattice.in', 'acceptrate': None,
                'examples': False, 'sizes': None, 'debug': False,
                'input': 'tests/enumeration/sc_1/enum.in_100',
                'polya': False, 'super': False, 'distribution': None,'seed':None,
                'filter':None,'visualize':True,'shapes':False,'show':False}
        _script_enum(args,testmode=True)   

class TestEnumOut(ut.TestCase):
    "Tests the _enum_out subroutine."""
    def _assertMultiLineEqual(self, first, second, msg=None):
        import difflib
        """Assert that two multi-line strings are equal.

        If they aren't, show a nice diff.
        code contributed by Ned Batchelder: http://stackoverflow.com/questions/3942820/how-to-do-unit-testing-of-functions-writing-files-using-python-unittest
        """
        self.assertTrue(isinstance(first, str),
                'First argument is not a string')
        self.assertTrue(isinstance(second, str),
                'Second argument is not a string')

        if first != second:
            message = ''.join(difflib.ndiff(first.splitlines(True),
                                                second.splitlines(True)))
            if msg:
                message += " : " + msg
            self.fail("Multi-line strings are unequal:\n" + message)

    def _compare_files(self,f1,f2):
        with open(f1,"r") as f:
            f1_content = f.readlines()
        with open(f2,"r") as f:
            f2_content = f.readlines()

        if len(f1_content) == len(f2_content):
            for i in range(1,len(f1_content)):
                self._assertMultiLineEqual(f1_content[i],f2_content[i])
        else:
            self.assertEqual(len(f1_content),len(f2_content))
                
    def test_enum_out1(self):
        from phenum.enumeration import _script_enum
        import sys
        from os import system

        args = {'profile': None, 'savedist': False, 'verbose': None,
                'outfile': "enum.out", 'enum': True,
                'lattice': 'tests/enumeration/enum_out/lattice.in_hcp', 'acceptrate': None,
                'examples': False, 'sizes': None, 'debug': False,
                'input': None,
                'polya': False, 'super': False, 'distribution': None,'seed':None,
                'filter':None,'visualize':False,'shapes':False,'show':False}

        _script_enum(args)
        if sys.version_info[0] < 3:
            self._compare_files("enum.out","tests/enumeration/enum_out/enum.out_hcp_p2")
        else:
            self._compare_files("enum.out","tests/enumeration/enum_out/enum.out_hcp_p3")
        system("rm enum.in polya.out.2 polya.out.1")

    def test_enum_out2(self):
        from phenum.enumeration import _script_enum
        from os import system
        import sys
        args = {'profile': None, 'savedist': False, 'verbose': None,
                'outfile': None, 'enum': True,
                'lattice': 'tests/enumeration/enum_out/lattice.in_sub1', 'acceptrate': None,
                'examples': False, 'sizes': None, 'debug': False,
                'input': 'tests/enumeration/enum_out/enum.in_sub1',
                'polya': False, 'super': False, 'distribution': None,'seed':0,
                'filter':None,'visualize':False,'shapes':False,'show':False}

        _script_enum(args)
        if sys.version_info[0] < 3:
            self._compare_files("enum.out","tests/enumeration/enum_out/enum.out_sub1_p2")
        else:
            self._compare_files("enum.out","tests/enumeration/enum_out/enum.out_sub1_p3")
        system("rm enum.out")

    def test_enum_out3(self):
        from phenum.enumeration import _script_enum
        from os import system
        import sys
        args = {'profile': None, 'savedist': False, 'verbose': None,
                'outfile': None, 'enum': True,
                'lattice': 'tests/enumeration/enum_out/lattice.in_ar', 'acceptrate': None,
                'examples': False, 'sizes': None, 'debug': False,
                'input': None,
                'polya': False, 'super': True, 'distribution': None,'seed': None,
                'filter':None,'visualize':False,'shapes':False,'show':False}

        _script_enum(args)
        if sys.version_info[0] < 3:
            self._compare_files("enum.out","tests/enumeration/enum_out/enum.out_ar_p2")
        else:
            self._compare_files("enum.out","tests/enumeration/enum_out/enum.out_ar_p3")
        system("rm enum.out enum.in polya.out.*")

    def test_enum_out4(self):
        from phenum.enumeration import _script_enum
        from os import system
        import sys
        args = {'profile': None, 'savedist': False, 'verbose': None,
                'outfile': None, 'enum': True,
                'lattice': 'tests/enumeration/enum_out/lattice.in_sub_ar', 'acceptrate': 0.5,
                'examples': False, 'sizes': None, 'debug': False,
                'input': "tests/enumeration/enum_out/enum.in_sub_ar",
                'polya': False, 'super': False, 'distribution': None,'seed': 0,
                'filter':None,'visualize':False,'shapes':False,'show':False}

        _script_enum(args)
        if sys.version_info[0] < 3:
            self._compare_files("enum.out","tests/enumeration/enum_out/enum.out_sub_ar_p2")
        else:
            self._compare_files("enum.out","tests/enumeration/enum_out/enum.out_sub_ar_p3")
        system("rm enum.out")

    def test_enum_out5(self):
        from phenum.enumeration import _script_enum
        from os import system
        import sys
        args = {'profile': None, 'savedist': False, 'verbose': None,
                'outfile': None, 'enum': True,
                'lattice': 'tests/enumeration/enum_out/lattice.in_guess_ar', 'acceptrate': None,
                'examples': False, 'sizes': None, 'debug': False,
                'input': "tests/enumeration/enum_out/enum.in_guess_ar",
                'polya': False, 'super': False, 'distribution': None,'seed': 0,
                'filter':None,'visualize':False,'shapes':False,'show':False}

        _script_enum(args)
        if sys.version_info[0] < 3:
            self._compare_files("enum.out","tests/enumeration/enum_out/enum.out_guess_ar_p2")
        else:
            self._compare_files("enum.out","tests/enumeration/enum_out/enum.out_guess_ar_p3")
        system("rm enum.out")
        
    def test_enum_out6(self):
        from phenum.enumeration import _script_enum
        from os import system
        import sys
        args = {'profile': None, 'savedist': False, 'verbose': None,
                'outfile': None, 'enum': True,
                'lattice': 'tests/enumeration/enum_out/lattice.in_guess_ar', 'acceptrate': 2.0,
                'examples': False, 'sizes': None, 'debug': False,
                'input': "tests/enumeration/enum_out/enum.in_guess_ar",
                'polya': False, 'super': False, 'distribution': None,'seed': 0,
                'filter':None,'visualize':False,'shapes':False,'show':False}

        with pytest.raises(ValueError):
            _script_enum(args)

    def test_enum_out7(self):
        from phenum.enumeration import _script_enum
        from os import system
        import sys
        args = {'profile': None, 'savedist': False, 'verbose': None,
                'outfile': None, 'enum': True,
                'lattice': 'tests/enumeration/enum_out/lattice.in_guess', 'acceptrate': None,
                'examples': False, 'sizes': None, 'debug': False,
                'input': "tests/enumeration/enum_out/enum.in_guess",
                'polya': False, 'super': False, 'distribution': None,'seed': 0,
                'filter':None,'visualize':False,'shapes':False,'show':False}

        _script_enum(args)
        if sys.version_info[0] < 3:
            self._compare_files("enum.out","tests/enumeration/enum_out/enum.out_guess_p2")
        else:
            self._compare_files("enum.out","tests/enumeration/enum_out/enum.out_guess_p3")
        system("rm enum.out")

    def test_enum_out8(self):
        from phenum.enumeration import _script_enum
        from os import system
        import sys
        args = {'profile': None, 'savedist': False, 'verbose': None,
                'outfile': None, 'enum': True,
                'lattice': 'tests/enumeration/enum_out/lattice.in_ar_sm', 'acceptrate': None,
                'examples': False, 'sizes': None, 'debug': False,
                'input': "tests/enumeration/enum_out/enum.in_ar_sm",
                'polya': False, 'super': False, 'distribution': None,'seed': 0,
                'filter':None,'visualize':False,'shapes':False,'show':False}

        _script_enum(args)
        if sys.version_info[0] < 3:
            self._compare_files("enum.out","tests/enumeration/enum_out/enum.out_ar_sm_p2")
        else:
            self._compare_files("enum.out","tests/enumeration/enum_out/enum.out_ar_sm_p3")
        system("rm enum.out")

    def test_enum_out9(self):
        from phenum.enumeration import _script_enum
        from os import system
        import sys
        args = {'profile': None, 'savedist': False, 'verbose': None,
                'outfile': None, 'enum': True,
                'lattice': 'tests/enumeration/enum_out/lattice.in_ar_sm2', 'acceptrate': None,
                'examples': False, 'sizes': None, 'debug': False,
                'input': "tests/enumeration/enum_out/enum.in_ar_sm2",
                'polya': False, 'super': False, 'distribution': None,'seed': 0,
                'filter':None,'visualize':False,'shapes':False,'show':False}

        _script_enum(args)
        if sys.version_info[0] < 3:
            self._compare_files("enum.out","tests/enumeration/enum_out/enum.out_ar_sm2_p2")
        else:
            self._compare_files("enum.out","tests/enumeration/enum_out/enum.out_ar_sm2_p3")
        system("rm enum.out")
        
class TestPolyaOut(ut.TestCase):
    """Tests the _polya_out subroutine."""

    def _assertMultiLineEqual(self, first, second, msg=None):
        import difflib
        """Assert that two multi-line strings are equal.

        If they aren't, show a nice diff.
        code contributed by Ned Batchelder: http://stackoverflow.com/questions/3942820/how-to-do-unit-testing-of-functions-writing-files-using-python-unittest
        """
        self.assertTrue(isinstance(first, str),
                'First argument is not a string')
        self.assertTrue(isinstance(second, str),
                'Second argument is not a string')

        if first != second:
            message = ''.join(difflib.ndiff(first.splitlines(True),
                                                second.splitlines(True)))
            if msg:
                message += " : " + msg
            self.fail("Multi-line strings are unequal:\n" + message)

    def _compare_files(self,f1,f2):
        with open(f1,"r") as f:
            f1_content = f.readlines()
        with open(f2,"r") as f:
            f2_content = f.readlines()

        if len(f1_content) == len(f2_content):
            for i in range(1,len(f1_content)):
                self._assertMultiLineEqual(f1_content[i],f2_content[i])
        else:
            self.assertEqual(len(f1_content),len(f2_content))
                
    def test_polya1(self):
        from phenum.enumeration import _script_enum

        args = {'profile': None, 'savedist': False, 'verbose': None,
                'outfile': None, 'enum': False,
                'lattice': "stuff", 'acceptrate': None,
                'examples': False, 'sizes': None, 'debug': False,
                'input': None,'polya': True, 'super': False, 'distribution': None,'seed':0,
                'filter':None,'visualize':False,'shapes':False,'show':False}
        with pytest.raises(IOError):
            _script_enum(args)

    def test_polya2(self):
        from phenum.enumeration import _script_enum
        from os import system

        args = {'profile': None, 'savedist': False, 'verbose': None,
                'outfile': None, 'enum': False,
                'lattice': "tests/enumeration/enum_out/lattice.in_hcp", 'acceptrate': None,
                'examples': False, 'sizes': None, 'debug': False,
                'input': 'tests/enumeration/enum_out/enum.in_sub1',
                'polya': True, 'super': False, 'distribution': None,'seed':0,
                'filter':None,'visualize':False,'shapes':False,'show':False}

        _script_enum(args)
        self._compare_files("polya.out.1","tests/enumeration/enum_out/hcp_polya.out.1")
        self._compare_files("polya.out.2","tests/enumeration/enum_out/hcp_polya.out.2")
        system("rm polya.out.1 polya.out.2")
