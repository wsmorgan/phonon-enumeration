"""Tests of the makeStr.py module."""
import unittest as ut
import os
import pytest

def get_sargs(args):
    """Returns the list of arguments parsed from sys.argv.
    """
    import sys
    sys.argv = args
    from phenum.makeStr import _parser_options
    return _parser_options()    

def test_examples():
    """Makes sure the script examples work properly.
    """
    argv = ["py.test", "-examples"]
    assert get_sargs(argv) is None


class TestMakeStructures(ut.TestCase):
    """Tests of the _make_structures subroutine."""

    def _compare_files(self,file1,file2):
        out1 = []
        out2 = []
        with open(file1,"r") as o1:
            for line in o1:
                out1.append(line.strip().split())
        with open(file2,"r") as o2:
            for line in o2:
                out2.append(line.strip().split())
        self.assertEqual(out1,out2)
                
    def test_str1(self):
        from phenum.makeStr import _make_structures
        from os import system
        args = {"structures":[10],
                "debug":False,
                "examples":False,
                "displace":0.0,
                "input":"tests/enumeration/sc_1/enum.out_100",
                "mink":True,
                "species":None,
                "verbose":None,
                "outfile":"vasp.{}",
                "rattle":0.0
                }
        _make_structures(args)
        self._compare_files("vasp.{}".format(args["structures"][0]),"tests/enumeration/sc_1/vasp.000010")
        system("rm vasp*")
                
    def test_str2(self):
        from phenum.makeStr import _make_structures
        from os import system
        args = {"structures":[20],
                "debug":False,
                "examples":False,
                "displace":0.0,
                "input":"tests/enumeration/sc_1/enum.out_100",
                "mink":True,
                "species":None,
                "verbose":None,
                "outfile":"vasp.{}",
                "rattle":0.0
                }
        _make_structures(args)
        self._compare_files("vasp.{}".format(args["structures"][0]),"tests/enumeration/sc_1/vasp.000020")
        system("rm vasp*")
                
    def test_str3(self):
        from phenum.makeStr import _make_structures
        from os import system
        args = {"structures":[33],
                "debug":False,
                "examples":False,
                "displace":0.0,
                "input":"tests/enumeration/sc_1/enum.out_100",
                "mink":True,
                "species":None,
                "verbose":None,
                "outfile":"vasp.{}",
                "rattle":0.0
                }
        _make_structures(args)
        self._compare_files("vasp.{}".format(args["structures"][0]),"tests/enumeration/sc_1/vasp.000033")
        system("rm vasp*")

    def test_str4(self):
        from phenum.makeStr import _make_structures
        from os import system
        args = {"structures":[1],
                "debug":False,
                "examples":False,
                "displace":0.0,
                "input":"tests/enumeration/fcc_1/enum.out_2_6",
                "mink":True,
                "species":None,
                "verbose":None,
                "outfile":"vasp.{}",
                "rattle":0.0
                }
        _make_structures(args)
        self._compare_files("vasp.{}".format(args["structures"][0]),"tests/enumeration/fcc_1/vasp.000001")
        system("rm vasp*")
                
    def test_str5(self):
        from phenum.makeStr import _make_structures
        from os import system
        args = {"structures":[55],
                "debug":False,
                "examples":False,
                "displace":0.0,
                "input":"tests/enumeration/fcc_1/enum.out_2_6",
                "mink":True,
                "species":None,
                "verbose":None,
                "outfile":"vasp.{}",
                "rattle":0.0
                }
        _make_structures(args)
        self._compare_files("vasp.{}".format(args["structures"][0]),"tests/enumeration/fcc_1/vasp.000055")
        system("rm vasp*")
                
    def test_str6(self):
        from phenum.makeStr import _make_structures
        from os import system
        args = {"structures":[50],
                "debug":False,
                "examples":False,
                "displace":0.0,
                "input":"tests/enumeration/fcc_1/enum.out_2_6",
                "mink":True,
                "species":None,
                "verbose":None,
                "outfile":"vasp.{}",
                "rattle":0.0
                }
        _make_structures(args)
        self._compare_files("vasp.{}".format(args["structures"][0]),"tests/enumeration/fcc_1/vasp.000050")
        system("rm vasp*")
                
    def test_str7(self):
        from phenum.makeStr import _make_structures
        from os import system
        args = {"structures":[88],
                "debug":False,
                "examples":False,
                "displace":0.0,
                "input":"tests/enumeration/fcc_1/enum.out_100_p2",
                "mink":True,
                "species":None,
                "verbose":None,
                "outfile":"vasp.{}",
                "rattle":0.0
                }
        _make_structures(args)
        self._compare_files("vasp.{}".format(args["structures"][0]),"tests/enumeration/fcc_1/vasp.000088")
        system("rm vasp*")
                
    def test_str8(self):
        from phenum.makeStr import _make_structures
        from os import system
        args = {"structures":[1],
                "debug":False,
                "examples":False,
                "displace":0.0,
                "input":"tests/enumeration/fcc_2/enum.out_3_4",
                "mink":True,
                "species":None,
                "verbose":None,
                "outfile":"vasp.{}",
                "rattle":0.0
                }
        _make_structures(args)
        self._compare_files("vasp.{}".format(args["structures"][0]),"tests/enumeration/fcc_2/vasp.000001")
        system("rm vasp*")
                
    def test_str9(self):
        from phenum.makeStr import _make_structures
        from os import system
        args = {"structures":[2],
                "debug":False,
                "examples":False,
                "displace":0.0,
                "input":"tests/enumeration/fcc_2/enum.out_3_4",
                "mink":True,
                "species":None,
                "verbose":None,
                "outfile":"vasp.{}",
                "rattle":0.0
                }
        _make_structures(args)
        self._compare_files("vasp.{}".format(args["structures"][0]),"tests/enumeration/fcc_2/vasp.000002")
        system("rm vasp*")
                
    def test_str10(self):
        from phenum.makeStr import _make_structures
        from os import system
        args = {"structures":[3],
                "debug":False,
                "examples":False,
                "displace":0.1,
                "input":"tests/enumeration/fcc_2/enum.out_3_4",
                "mink":True,
                "species":None,
                "verbose":None,
                "outfile":"vasp.{}",
                "rattle":0.0
                }
        _make_structures(args)
        self._compare_files("vasp.{}".format(args["structures"][0]),"tests/enumeration/fcc_2/vasp.000003")
        system("rm vasp*")

    def test_str11(self):
        from phenum.makeStr import _make_structures
        from os import system
        args = {"structures":[3],
                "debug":False,
                "examples":False,
                "displace":0.0,
                "input":"tests/enumeration/fcc_2/enum.out_3_4",
                "mink":True,
                "species":['Ni','Al','Cu'],
                "verbose":None,
                "outfile":"vasp.{}",
                "rattle":0.0
                }
        _make_structures(args)
        self._compare_files("vasp.{}".format(args["structures"][0]),"tests/enumeration/fcc_2/vasp.3.NiAlCu")
        system("rm vasp*")

    def test_str12(self):
        from phenum.makeStr import _make_structures
        from os import system
        args = {"structures":[3],
                "debug":False,
                "examples":False,
                "displace":0.0,
                "input":"tests/enumeration/fcc_2/enum.out_3_4",
                "mink":True,
                "species":['Co','W','V'],
                "verbose":None,
                "outfile":"vasp.{}",
                "rattle":0.0
                }
        _make_structures(args)
        self._compare_files("vasp.{}".format(args["structures"][0]),"tests/enumeration/fcc_2/vasp.3.CoWV")
        system("rm vasp*")

    def test_str13(self):
        from phenum.makeStr import _make_structures
        from os import system
        args = {"structures":[3],
                "debug":False,
                "examples":False,
                "displace":0.0,
                "input":"tests/enumeration/fcc_1/enum.out_2_6",
                "mink":True,
                "species":['Ti','S'],
                "verbose":None,
                "outfile":"vasp.{}",
                "rattle":0.0
                }
        _make_structures(args)
        self._compare_files("vasp.{}".format(args["structures"][0]),"tests/enumeration/fcc_1/vasp.3.TiS")
        system("rm vasp*")

    def test_str14(self):
        from phenum.makeStr import _make_structures
        from os import system
        args = {"structures":[3],
                "debug":False,
                "examples":False,
                "displace":0.0,
                "input":"tests/enumeration/fcc_1/enum.out_2_6",
                "mink":True,
                "species":['H','Pt'],
                "verbose":None,
                "outfile":"vasp.{}",
                "rattle":0.0
                }
        _make_structures(args)
        self._compare_files("vasp.{}".format(args["structures"][0]),"tests/enumeration/fcc_1/vasp.3.HPt")
        system("rm vasp*")

    def test_str15(self):
        from phenum.makeStr import run
        from os import system
        args = {"structures":[3],
                "debug":False,
                "examples":False,
                "displace":0.0,
                "input":"tests/enumeration/fcc_1/enum.out_2_6",
                "mink":True,
                "species":['H','Pt'],
                "verbose":None,
                "outfile":"vasp.{}",
                "rattle":0.0
                }
        run(args)
        self._compare_files("vasp.{}".format(args["structures"][0]),"tests/enumeration/fcc_1/vasp.3.HPt")
        system("rm vasp*")

    def test_str16(self):
        from phenum.makeStr import run
        from os import system
        args = {"structures": None,
                "debug":False,
                "examples":False,
                "displace":0.0,
                "input":"tests/enumeration/fcc_1/enum.out_2_6",
                "mink":True,
                "species":['H','Pt'],
                "verbose":None,
                "outfile":"vasp.{}",
                "rattle":0.0
                }
        with pytest.raises(ValueError):
            run(args)
        

    def test_str17(self):
        from phenum.makeStr import run
        from os import system
        args = {"structures": ['bite'],
                "debug":False,
                "examples":False,
                "displace":0.0,
                "input":"tests/enumeration/fcc_1/enum.out_2_6",
                "mink":True,
                "species":['H','Pt'],
                "verbose":None,
                "outfile":"vasp.{}",
                "rattle":0.0
                }
        with pytest.raises(ValueError):
            run(args)

    def test_str18(self):
        from phenum.makeStr import run
        from os import system
        args = {"structures": None,
                "debug":False,
                "examples":False,
                "displace":0.0,
                "input":"tests/enumeration/fcc_1/enum.out_2_6",
                "mink":True,
                "species":['H','Pt'],
                "verbose":None,
                "outfile":"vasp.{}",
                "rattle":0.0
                }
        with pytest.raises(ValueError):
            run(args)

    def test_str19(self):
        from phenum.makeStr import run
        from os import system
        args = {"structures":['all'],
                "debug":False,
                "examples":False,
                "displace":0.0,
                "input":"tests/enumeration/fcc_1/enum.out_2_6",
                "mink":True,
                "species":['H','Pt'],
                "verbose":None,
                "outfile":"vasp.{}",
                "rattle":0.0
                }
        run(args)
        system("rm vasp*")

    def test_str20(self):
        from phenum.makeStr import run
        from os import system
        args = {"structures":['1','3'],
                "debug":False,
                "examples":False,
                "displace":0.0,
                "input":"tests/enumeration/fcc_1/enum.out_2_6",
                "mink":True,
                "species":['H','Pt'],
                "verbose":None,
                "outfile":"vasp.{}",
                "rattle":0.0
                }
        run(args)
        system("rm vasp*")

    def test_str21(self):
        from phenum.makeStr import run
        from os import system
        args = {"structures":['1','175'],
                "debug":False,
                "examples":False,
                "displace":0.0,
                "input":"tests/enumeration/hcp_1/enum.out_1_4",
                "mink":True,
                "species":None,
                "verbose":None,
                "outfile":"vasp.{}",
                "rattle":0.0
                }
        run(args)
        for i in [1,45,90,175]:
            self._compare_files("vasp.{}".format(i),"tests/enumeration/hcp_1/vasp.{}.fin".format(i))
        system("rm vasp*")


    def test_str22(self):
        from phenum.makeStr import run
        from os import system
        args = None
        with pytest.raises(ValueError):
            run(args)
        
