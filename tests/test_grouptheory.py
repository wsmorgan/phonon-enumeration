"""Methods for testing the subroutines in the grouptheory module."""
import unittest as ut
from phenum.grouptheory import ArrowPerm, RotPermList, OpList
import pytest
import numpy as np

gpath =  "tests/grouptheory/"

def _read_fixOp_1D(fname):
    import os
    i = 1
    growing = True
    out = []
    while growing:
        if os.path.isfile(fname+"/_-"+str(i)+"-rot") or os.path.isfile(fname+"/_-"+str(i)+"-shift"):
            i += 1
        else:
            growing = False

    for j in range(1,i):
        if os.path.isfile(fname+"/_-"+str(j)+"-rot"):
            rot = [np.transpose(t) for t in _read_float_3D(fname+"/_-"+str(j)+"-rot")]
        else:
            rot = None
        if os.path.isfile(fname+"/_-"+str(j)+"-shift"):
            shift = list(map(list,zip(*_read_float_2D(fname+"/_-"+str(j)+"-shift"))))
        else:
            shift = None
        temp = OpList(rot=rot,shift=shift)

        out.append(temp)
    return out

def _read_RotPermList_1D(fname,arrowp = None):
    import os
    i = 1
    growing = True
    out = []
    
    while growing :
        if os.path.isfile(fname+"/_-"+str(i)+"-nL") or os.path.isfile(fname+"/_-"+str(i)+"-v") or os.path.isfile(fname+"/_-"+str(i)+"-RotIndx") or os.path.isfile(fname+"/_-"+str(i)+"-perm"):
            i += 1
        else:
            growing = False
            
    for j in range(1,i):
        if os.path.isfile(fname+"/_-"+str(j)+"-nL"):
            nL = _read_int(fname+"/_-"+str(j)+"-nL")
        else:
            nL = None
        if os.path.isfile(fname+"/_-"+str(j)+"-v"):
            v = _read_float_3D(fname+"/_-"+str(j)+"-v")
        else:
            v = None
        if os.path.isfile(fname+"/_-"+str(j)+"-perm"):
            perm = _read_int_2D(fname+"/_-"+str(j)+"-perm")
            perm = [[i-1 for i in t] for t in perm]
        else:
            perm = None
        if arrowp == None:
            a_perm = None
        if os.path.isfile(fname+"/_-"+str(j)+"-RotIndx"):
            RotIndx = _read_int_1D(fname+"/_-"+str(j)+"-RotIndx")
            RotIndx = [i-1 for i in RotIndx]
        else:
            RotIndx = None
        temp = RotPermList(nL = nL, v = v, perm = perm, arrows=a_perm, RotIndx= RotIndx)

        out.append(temp)
    return out
    
def _read_fixOp(fname):
    import os
    if os.path.isfile(fname+"/_-rot"):
        rot = _read_float_3D(fname+"/_-rot")
    else:
        rot = None
    if os.path.isfile(fname+"/_-shift"):
        shift = list(map(list,zip(*_read_float_2D(fname+"/_-shift"))))
    else:
        shift = None
    out = OpList(rot=rot,shift=shift)
    return out

def _read_RotPermList(fname,arrowp = None):
    import os
    if os.path.isfile(fname+"/_-nL"):
        nL = _read_int(fname+"/_-nL")
    else:
        nL = None
    if os.path.isfile(fname+"/_-v"):
        v = _read_float_3D(fname+"/_-v")
    else:
        v = None
    if os.path.isfile(fname+"/_-perm"):
        perm = _read_int_2D(fname+"/_-perm")
        perm = [[i-1 for i in j] for j in perm]

    else:
        perm = None
    if arrowp == None:
        a_perm = None
    if os.path.isfile(fname+"/_-RotIndx"):
        RotIndx = _read_int_1D(fname+"/_-RotIndx")
        RotIndx = [i-1 for i in RotIndx]
    else:
        RotIndx = None
    out = RotPermList(nL = nL, v = v, perm = perm, arrows=a_perm, RotIndx= RotIndx)
    return out

def _read_float_3D(fname):
        with open(fname,"r") as inf:
            temp = inf.readline()
            sizes = inf.readline()
            sizes = [int(x) for x in sizes.strip().split() if x !="##"]
            temp = inf.readline()
            in_data = []
            in_temp = []
            for line in inf:
                if "#" not in line:
                    in_temp.append([float(i) for i in line.strip().split()])
                else:
                    in_data.append(in_temp)
                    in_temp = []
            in_data.append(in_temp)
            
        out = []
        for i  in range(sizes[2]):
            out_t = []
            for j in range(sizes[1]):
                out_t.append([k[j][i] for k in in_data])
            out.append(out_t)

        return(out)

def _read_int_3D(fname):
        with open(fname,"r") as inf:
            temp = inf.readline()
            sizes = inf.readline()
            sizes = [int(x) for x in sizes.strip().split() if x !="##"]
            temp = inf.readline()
            in_data = []
            in_temp = []
            for line in inf:
                if "#" not in line:
                    in_temp.append([int(i) for i in line.strip().split()])
                else:
                    in_data.append(in_temp)
                    in_temp = []
            in_data.append(in_temp)
            
        out = []
        for i  in range(sizes[2]):
            out_t = []
            for j in range(sizes[1]):
                out_t.append([k[j][i] for k in in_data])
            out.append(np.transpose(out_t))

        return(out)

def _read_output(test):
    values = []
    with open("tests/grouptheory/"+test) as f:
        for line in f:
            values.append(eval(line))
    return values

def _read_float_2D(fname):
    array = []
    with open(fname,"r") as f1:
        for line in f1:
            if "#" not in line:
                array.append([float(i) for i in line.strip().split()])
    return array

def _read_float_1D(fname):
    array = []
    from os import getcwd
    with open(fname,"r") as f1:
        for line in f1:
            if "#" not in line:
                array = [float(i) for i in line.strip().split()]
    return array

def _read_int_2D(fname):
    array = []
    with open(fname,"r") as f1:
        for line in f1:
            if "#" not in line:
                array.append([int(i) for i in line.strip().split()])
    return array

def _read_int_1D(fname):
    array = []
    with open(fname,"r") as f1:
        for line in f1:
            if "#" not in line:
                array = [int(i) for i in line.strip().split()]
    return array

def _read_int(fname):
    with open(fname,"r") as f1:
        line = f1.readline()
        if "#" in line:
            line = f1.readline()
        val = int(line.strip())
    return val
    
def _read_float(fname):
    with open(fname,"r") as f1:
        line = f1.readline()
        if "#" in line:
            line = f1.readline()
        val = float(line.strip())
    return val    

def _read_logical(fname):
    with open(fname,"r") as f1:
        line = f1.readline()
        if "#" in line:
            line = f1.readline()
    if "t" in line.lower():
        val = True
    else:
        val = False
    return val

class TestGetFullHNF(ut.TestCase):
    """ Tests of the get_full_HNF subroutine."""

    def test_1(self):
        from phenum.grouptheory import get_full_HNF
        from numpy import array
        HNF = array([1,0,1,0,0,1])
        out = [[1,0,0],[0,1,0],[0,0,1]]
        self.assertEqual(get_full_HNF(HNF),out)

    def test_2(self):
        from phenum.grouptheory import get_full_HNF
        from numpy import array
        HNF = array([2,1,2,1,0,4])
        out = [[2,0,0],[1,2,0],[1,0,4]]
        self.assertEqual(get_full_HNF(HNF),out)

    def test_3(self):
        from phenum.grouptheory import get_full_HNF
        from numpy import array
        HNF = array([1,0,3,1,2,3])
        out = [[1,0,0],[0,3,0],[1,2,3]]
        self.assertEqual(get_full_HNF(HNF),out)

    def test_4(self):
        from phenum.grouptheory import get_full_HNF
        from numpy import array
        HNF = [0,0,0,0,0,0]
        out = [[0,0,0],[0,0,0],[0,0,0]]
        self.assertEqual(get_full_HNF(HNF),out)

    def test_5(self):
        from phenum.grouptheory import get_full_HNF
        from numpy import array
        HNF = array([3,0,3,0,0,3])
        out = [[3,0,0],[0,3,0],[0,0,3]]
        self.assertEqual(get_full_HNF(HNF),out)

    def test_1(self):
        from phenum.grouptheory import get_full_HNF
        from numpy import array
        HNF = array([1,1,2,0,2,2])
        out = [[1,0,0],[1,2,0],[0,2,2]]
        self.assertEqual(get_full_HNF(HNF),out)

    def test_1(self):
        from phenum.grouptheory import get_full_HNF
        from numpy import array
        HNF = array([2,0,2,0,2,4])
        out = [[2,0,0],[0,2,0],[0,2,4]]
        self.assertEqual(get_full_HNF(HNF),out)

class TestSmithNormalForm(ut.TestCase):
    """ Tests of the SmithNormalForm subroutine."""

    def test_1(self):
        from phenum.grouptheory import SmithNormalForm
        HNF =  [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
        out =  ([[1, 0, 0], [0, 1, 0], [0, 0, 1]], [[1, 0, 0], [0, 1, 0], [0, 0, 1]], [[1, 0, 0], [0, 1, 0], [0, 0, 1]])
        self.assertEqual(SmithNormalForm(HNF),out)

    def test_2(self):
        from phenum.grouptheory import SmithNormalForm
        HNF =  [[1, 0, 0], [0, 1, 0], [0, 1, 2]]
        out =  ([[1, 0, 0], [0, 1, 0], [0, 0, 2]], [[1, 0, 0], [0, 1, 0], [0, -1, 1]], [[1, 0, 0], [0, 1, 0], [0, 0, 1]])
        self.assertEqual(SmithNormalForm(HNF),out)

    def test_3(self):
        from phenum.grouptheory import SmithNormalForm
        HNF =  [[1, 0, 0], [0, 1, 0], [0, 0, 3]]
        out =  ([[1, 0, 0], [0, 1, 0], [0, 0, 3]], [[1, 0, 0], [0, 1, 0], [0, 0, 1]], [[1, 0, 0], [0, 1, 0], [0, 0, 1]])
        self.assertEqual(SmithNormalForm(HNF),out)

    def test_4(self):
        from phenum.grouptheory import SmithNormalForm
        HNF =  [[1, 0, 0], [0, 2, 0], [0, 0, 2]]
        out =  ([[1, 0, 0], [0, 2, 0], [0, 0, 2]], [[1, 0, 0], [0, 1, 0], [0, 0, 1]], [[1, 0, 0], [0, 1, 0], [0, 0, 1]])
        self.assertEqual(SmithNormalForm(HNF),out)

    def test_5(self):
        from phenum.grouptheory import SmithNormalForm
        HNF =  [[1, 0, 0], [0, 1, 0], [1, 2, 5]]
        out =  ([[1, 0, 0], [0, 1, 0], [0, 0, 5]], [[1, 0, 0], [0, 1, 0], [-1, -2, 1]], [[1, 0, 0], [0, 1, 0], [0, 0, 1]])
        self.assertEqual(SmithNormalForm(HNF),out)

    def test_6(self):
        from phenum.grouptheory import SmithNormalForm
        HNF =  [[1, 0, 0], [0, 1, 0], [2, 3, 6]]
        out =  ([[1, 0, 0], [0, 1, 0], [0, 0, 6]], [[1, 0, 0], [0, 1, 0], [-2, -3, 1]], [[1, 0, 0], [0, 1, 0], [0, 0, 1]])
        self.assertEqual(SmithNormalForm(HNF),out)

    def test_7(self):
        from phenum.grouptheory import SmithNormalForm
        HNF =  [[1, 0, 0], [0, 1, 0], [0, 6, 7]]
        out =  ([[1, 0, 0], [0, 1, 0], [0, 0, 7]], [[1, 0, 0], [0, 1, 0], [0, -6, 1]], [[1, 0, 0], [0, 1, 0], [0, 0, 1]])
        self.assertEqual(SmithNormalForm(HNF),out)

    def test_8(self):
        from phenum.grouptheory import SmithNormalForm
        HNF =  [[1, 0, 0], [1, 2, 0], [1, 0, 4]]
        out =  ([[1, 0, 0], [0, 2, 0], [0, 0, 4]], [[1, 0, 0], [-1, 1, 0], [-1, 0, 1]], [[1, 0, 0], [0, 1, 0], [0, 0, 1]])
        self.assertEqual(SmithNormalForm(HNF),out)

    def test_9(self):
        from phenum.grouptheory import SmithNormalForm
        HNF =  [[2, 0, 0], [0, 2, 0], [0, 0, 2]]
        out =  ([[2, 0, 0], [0, 2, 0], [0, 0, 2]], [[1, 0, 0], [0, 1, 0], [0, 0, 1]], [[1, 0, 0], [0, 1, 0], [0, 0, 1]])
        self.assertEqual(SmithNormalForm(HNF),out)

    def test_10(self):
        from phenum.grouptheory import SmithNormalForm
        HNF =  [[1, 0, 0], [0, 1, 0], [1, 5, 10]]
        out =  ([[1, 0, 0], [0, 1, 0], [0, 0, 10]], [[1, 0, 0], [0, 1, 0], [-1, -5, 1]], [[1, 0, 0], [0, 1, 0], [0, 0, 1]])
        self.assertEqual(SmithNormalForm(HNF),out)

    def test_11(self):
        from phenum.grouptheory import SmithNormalForm
        HNF =  [[-1, 0, 0], [0, 1, 0], [0, 0, 1]]
        out =  ()
        with pytest.raises(ValueError):
            self.assertEqual(SmithNormalForm(HNF),out)

    def test_12(self):
        from phenum.grouptheory import SmithNormalForm
        HNF =  [[0, 1, 0], [0, 0, 1], [1, 0, 0]]
        out =  ([[1, 0, 0], [0, 1, 0], [0, 0, 1]], [[0, 0, 1], [1, 0, 0], [0, 1, 0]], [[1, 0, 0], [0, 1, 0], [0, 0, 1]]) 
        self.assertEqual(SmithNormalForm(HNF),out)

    def test_13(self):
        from phenum.grouptheory import SmithNormalForm
        HNF =  [[1, -1, -2], [1, 2, -3], [1, 2, 4]]
        out =  ([[1, 0, 0], [0, 1, 0], [0, 0, 21]], [[1, 0, 0], [-1, 1, 0], [-7, 6, 1]], [[1, -2, 7], [0, 0, 1], [0, -1, 3]])
        self.assertEqual(SmithNormalForm(HNF),out)

    def test_14(self):
        from phenum.grouptheory import SmithNormalForm
        HNF =  [[-1, -2, -3], [-1, -1, -2], [-1, -2, -4]]
        out =  ([[1, 0, 0], [0, 1, 0], [0, 0, 1]], [[1, 0, 0], [-1, 1, 0], [-1, 0, 1]], [[-1, -2, 1], [0, 1, 1], [0, 0, -1]])
        self.assertEqual(SmithNormalForm(HNF),out)

    def test_15(self):
        from phenum.grouptheory import SmithNormalForm
        HNF =  [[1, 2.5, 0], [0, 1.5, 1.66], [1.5, 1.25, 1.3]]
        out =  ([[0.5, 0.0, 0.0], [0.0, 0.5, 0.0], [0.0, 0.0, 16.5]], [[-1.0, 0.0, 1.0], [3.0, -3.0, -2.0], [-9.0, 10.0, 6.0]], [[1, 2.5, 23.25], [0, 1.0, 10.5], [0, 0.0, 1.0]])
        with pytest.raises(RuntimeError):
            SmithNormalForm(HNF)

    def test_16(self):
        from phenum.grouptheory import SmithNormalForm
        HNF =  [[2, 0, 0], [0, 2, 0], [0, 0, 1]]
        out =  ([[1, 0, 0], [0, 2, 0], [0, 0, 2]], [[1, 0, 1], [0, 1, 0], [-1, 0, 0]], [[0, 0, -1], [0, 1, 0], [1, 0, 2]])
        self.assertEqual(SmithNormalForm(HNF),out)

    def test_17(self):
        """Test of the bug reported in issue #56."""
        from phenum.grouptheory import SmithNormalForm
        HNF = [[1,2,4],[3,3,4],[3,4,2]]
        S, L, R = SmithNormalForm(HNF)
        self.assertTrue(np.allclose(list(np.dot(np.dot(L,HNF),R)),S))

class TestAGroup(ut.TestCase):
    """ Tests of the a_group subroutine."""

    def test_1(self):
        from phenum.grouptheory import a_group
        trans = [[0,1],[1,0]]
        rots = [[[0,1],[0,1,2,3,4,5]],[[1,0],[2,3,0,1,5,4]],[[1,0],[2,1,0,3,5,4]],[[0,1],[0,3,2,1,5,4]]]
        out = _read_output("agroup.out.1")
        self.assertEqual(a_group(trans,rots),out)

    def test_2(self):
        from phenum.grouptheory import a_group
        trans = [[j-1 for j in i] for i in [[1, 2, 3, 4], [2, 1, 4, 3], [3, 4, 1, 2], [4, 3, 2, 1]]]
        rots = [[[j-1 for j in i] for i in t] for t in [[[1, 2, 3, 4], [1, 2, 3, 4, 5, 6]], [[1, 4, 3, 2], [1, 3, 2, 4, 6, 5]], [[1, 2, 3, 4], [4, 2, 3, 1, 5, 6]], [[1, 4, 3, 2], [4, 3, 2, 1, 6, 5]], [[1, 2, 3, 4], [1, 5, 3, 4, 2, 6]], [[1, 4, 3, 2], [1, 3, 5, 4, 6, 2]], [[1, 2, 3, 4], [4, 5, 3, 1, 2, 6]], [[1, 4, 3, 2], [4, 3, 5, 1, 6, 2]], [[1, 2, 3, 4], [1, 2, 6, 4, 5, 3]], [[1, 4, 3, 2], [1, 6, 2, 4, 3, 5]], [[1, 2, 3, 4], [4, 2, 6, 1, 5, 3]], [[1, 4, 3, 2], [4, 6, 2, 1, 3, 5]], [[1, 2, 3, 4], [1, 5, 6, 4, 2, 3]], [[1, 4, 3, 2], [1, 6, 5, 4, 3, 2]], [[1, 2, 3, 4], [4, 5, 6, 1, 2, 3]], [[1, 4, 3, 2], [4, 6, 5, 1, 3, 2]]]]
        out = _read_output("agroup.out.2")
        self.assertEqual(a_group(trans,rots),out)

    def test_3(self):
        from phenum.grouptheory import a_group
        trans = [[j-1 for j in i] for i in [[1, 2, 3, 4, 5, 6, 7, 8], [2, 1, 4, 3, 6, 5, 8, 7], [3, 4, 5, 6, 7, 8, 1, 2], [4, 3, 6, 5, 8, 7, 2, 1], [5, 6, 7, 8, 1, 2, 3, 4], [6, 5, 8, 7, 2, 1, 4, 3], [7, 8, 1, 2, 3, 4, 5, 6], [8, 7, 2, 1, 4, 3, 6, 5]]]
        rots = [[[j-1 for j in i] for i in t] for t in [[[1, 2, 3, 4, 5, 6, 7, 8], [1, 2, 3, 4, 5, 6]], [[1, 2, 3, 4, 5, 6, 7, 8], [4, 2, 3, 1, 5, 6]], [[1, 2, 3, 4, 5, 6, 7, 8], [1, 5, 3, 4, 2, 6]], [[1, 2, 3, 4, 5, 6, 7, 8], [4, 5, 3, 1, 2, 6]], [[1, 2, 7, 8, 5, 6, 3, 4], [1, 2, 6, 4, 5, 3]], [[1, 2, 7, 8, 5, 6, 3, 4], [4, 2, 6, 1, 5, 3]], [[1, 2, 7, 8, 5, 6, 3, 4], [1, 5, 6, 4, 2, 3]], [[1, 2, 7, 8, 5, 6, 3, 4], [4, 5, 6, 1, 2, 3]]]]
        out = _read_output("agroup.out.3")
        self.assertEqual(a_group(trans,rots),out)

    def test_4(self):
        from phenum.grouptheory import a_group
        trans =[[j - 1 for j in i] for i in[[1,2,3,4,5,6,7,8], [2,1,4,3,6,5,8,7], [3,4,5,6,7,8,1,2], [4,3,6,5,8,7,2,1], [5,6,7,8,1,2,3,4], [6,5,8,7,2,1,4,3], [7,8,1,2,3,4,5,6], [8,7,2,1,4,3,6,5]]]
        rots = [[[0,1,2,3,4,5,6,7],[0,1,2,3]],[[0,1,2,3,4,5,6,7],[2,1,0,3]],[[0,1,6,7,4,5,2,3],[0,3,2,1]],[[0,1,6,7,4,5,2,3],[2,3,0,1]]]
        out = _read_output("agroup.out.4")
        self.assertEqual(a_group(trans,rots),out)

    def test_5(self):
        from phenum.grouptheory import a_group
        trans =[[j - 1 for j in i] for i in[[1, 2, 3, 4, 5, 6, 7, 8], [2, 1, 4, 3, 6, 5, 8, 7], [3, 4, 5, 6, 7, 8, 1, 2], [4, 3, 6, 5, 8, 7, 2, 1], [5, 6, 7, 8, 1, 2, 3, 4], [6, 5, 8, 7, 2, 1, 4, 3], [7, 8, 1, 2, 3, 4, 5, 6], [8, 7, 2, 1, 4, 3, 6, 5]]]
        rots = [[[0,1,2,3,4,5,6,7],[0,1,2,3]],[[0,1,2,3,4,5,6,7],[2,1,0,3]],[[0,1,6,7,4,5,2,3],[0,3,2,1]],[[0,1,6,7,4,5,2,3],[2,3,0,1]]]
        out = _read_output("agroup.out.5")
        self.assertEqual(a_group(trans,rots),out)

class TestAGroupGen(ut.TestCase):
    """ Tests of the a_group subroutine."""

    def test_1(self):
        from phenum.grouptheory import a_group_gen
        trans = [[0,1],[1,0]]
        rots = [[[0,1],[0,1,2,3,4,5]],[[1,0],[2,3,0,1,5,4]],[[1,0],[2,1,0,3,5,4]],[[0,1],[0,3,2,1,5,4]]]
        out = _read_output("agroupgen.out.1")
        self.assertEqual(a_group_gen(trans,rots),out)

    def test_2(self):
        from phenum.grouptheory import a_group_gen
        trans = [[j-1 for j in i] for i in [[1, 2, 3, 4], [2, 1, 4, 3], [3, 4, 1, 2], [4, 3, 2, 1]]]
        rots = [[[j-1 for j in i] for i in t] for t in [[[1, 2, 3, 4], [1, 2, 3, 4, 5, 6]], [[1, 4, 3, 2], [1, 3, 2, 4, 6, 5]], [[1, 2, 3, 4], [4, 2, 3, 1, 5, 6]], [[1, 4, 3, 2], [4, 3, 2, 1, 6, 5]], [[1, 2, 3, 4], [1, 5, 3, 4, 2, 6]], [[1, 4, 3, 2], [1, 3, 5, 4, 6, 2]], [[1, 2, 3, 4], [4, 5, 3, 1, 2, 6]], [[1, 4, 3, 2], [4, 3, 5, 1, 6, 2]], [[1, 2, 3, 4], [1, 2, 6, 4, 5, 3]], [[1, 4, 3, 2], [1, 6, 2, 4, 3, 5]], [[1, 2, 3, 4], [4, 2, 6, 1, 5, 3]], [[1, 4, 3, 2], [4, 6, 2, 1, 3, 5]], [[1, 2, 3, 4], [1, 5, 6, 4, 2, 3]], [[1, 4, 3, 2], [1, 6, 5, 4, 3, 2]], [[1, 2, 3, 4], [4, 5, 6, 1, 2, 3]], [[1, 4, 3, 2], [4, 6, 5, 1, 3, 2]]]]
        out = _read_output("agroupgen.out.2")
        self.assertEqual(a_group_gen(trans,rots),out)

    def test_3(self):
        from phenum.grouptheory import a_group_gen
        trans = [[j-1 for j in i] for i in [[1, 2, 3, 4, 5, 6, 7, 8], [2, 1, 4, 3, 6, 5, 8, 7], [3, 4, 5, 6, 7, 8, 1, 2], [4, 3, 6, 5, 8, 7, 2, 1], [5, 6, 7, 8, 1, 2, 3, 4], [6, 5, 8, 7, 2, 1, 4, 3], [7, 8, 1, 2, 3, 4, 5, 6], [8, 7, 2, 1, 4, 3, 6, 5]]]
        rots = [[[j-1 for j in i] for i in t] for t in [[[1, 2, 3, 4, 5, 6, 7, 8], [1, 2, 3, 4, 5, 6]], [[1, 2, 3, 4, 5, 6, 7, 8], [4, 2, 3, 1, 5, 6]], [[1, 2, 3, 4, 5, 6, 7, 8], [1, 5, 3, 4, 2, 6]], [[1, 2, 3, 4, 5, 6, 7, 8], [4, 5, 3, 1, 2, 6]], [[1, 2, 7, 8, 5, 6, 3, 4], [1, 2, 6, 4, 5, 3]], [[1, 2, 7, 8, 5, 6, 3, 4], [4, 2, 6, 1, 5, 3]], [[1, 2, 7, 8, 5, 6, 3, 4], [1, 5, 6, 4, 2, 3]], [[1, 2, 7, 8, 5, 6, 3, 4], [4, 5, 6, 1, 2, 3]]]]
        out = _read_output("agroupgen.out.3")
        self.assertEqual(a_group_gen(trans,rots),out)

    def test_4(self):
        from phenum.grouptheory import a_group_gen
        trans =[[j - 1 for j in i] for i in[[1,2,3,4,5,6,7,8], [2,1,4,3,6,5,8,7], [3,4,5,6,7,8,1,2], [4,3,6,5,8,7,2,1], [5,6,7,8,1,2,3,4], [6,5,8,7,2,1,4,3], [7,8,1,2,3,4,5,6], [8,7,2,1,4,3,6,5]]]
        rots = [[[0,1,2,3,4,5,6,7],[0,1,2,3]],[[0,1,2,3,4,5,6,7],[2,1,0,3]],[[0,1,6,7,4,5,2,3],[0,3,2,1]],[[0,1,6,7,4,5,2,3],[2,3,0,1]]]
        out = _read_output("agroupgen.out.4")
        self.assertEqual(a_group_gen(trans,rots),out)

    def test_5(self):
        from phenum.grouptheory import a_group_gen
        trans =[[j - 1 for j in i] for i in[[1, 2, 3, 4, 5, 6, 7, 8], [2, 1, 4, 3, 6, 5, 8, 7], [3, 4, 5, 6, 7, 8, 1, 2], [4, 3, 6, 5, 8, 7, 2, 1], [5, 6, 7, 8, 1, 2, 3, 4], [6, 5, 8, 7, 2, 1, 4, 3], [7, 8, 1, 2, 3, 4, 5, 6], [8, 7, 2, 1, 4, 3, 6, 5]]]
        rots = [[[0,1,2,3,4,5,6,7],[0,1,2,3]],[[0,1,2,3,4,5,6,7],[2,1,0,3]],[[0,1,6,7,4,5,2,3],[0,3,2,1]],[[0,1,6,7,4,5,2,3],[2,3,0,1]]]
        out = _read_output("agroupgen.out.5")
        self.assertEqual(a_group_gen(trans,rots),out)

class TestMakeMemberList(ut.TestCase):
    """Tests of the _make_member_list subroutine."""

    def test_1(self):
        from phenum.grouptheory import _make_member_list
        case = 1
        n = _read_float_1D(gpath+"make_member_list_n.in."+str(case))
        out = list(map(list,zip(*_read_float_2D(gpath+"make_member_list_p.out."+str(case)))))
        self.assertTrue(np.allclose(_make_member_list(n),out))

    def test_2(self):
        from phenum.grouptheory import _make_member_list
        case = 2
        n = _read_float_1D(gpath+"make_member_list_n.in."+str(case))
        out = list(map(list,zip(*_read_float_2D(gpath+"make_member_list_p.out."+str(case)))))
        self.assertTrue(np.allclose(_make_member_list(n),out))

    def test_3(self):
        from phenum.grouptheory import _make_member_list
        case = 3
        n = _read_float_1D(gpath+"make_member_list_n.in."+str(case))
        out = list(map(list,zip(*_read_float_2D(gpath+"make_member_list_p.out."+str(case)))))
        self.assertTrue(np.allclose(_make_member_list(n),out))

    def test_4(self):
        from phenum.grouptheory import _make_member_list
        case = 4
        n = _read_float_1D(gpath+"make_member_list_n.in."+str(case))
        out = list(map(list,zip(*_read_float_2D(gpath+"make_member_list_p.out."+str(case)))))
        self.assertTrue(np.allclose(_make_member_list(n),out))

    def test_5(self):
        from phenum.grouptheory import _make_member_list
        case = 5
        n = _read_float_1D(gpath+"make_member_list_n.in."+str(case))
        out = list(map(list,zip(*_read_float_2D(gpath+"make_member_list_p.out."+str(case)))))
        self.assertTrue(np.allclose(_make_member_list(n),out))

    def test_6(self):
        from phenum.grouptheory import _make_member_list
        case = 6
        n = _read_float_1D(gpath+"make_member_list_n.in."+str(case))
        out = list(map(list,zip(*_read_float_2D(gpath+"make_member_list_p.out."+str(case)))))
        self.assertTrue(np.allclose(_make_member_list(n),out))

    def test_7(self):
        from phenum.grouptheory import _make_member_list
        case = 7
        n = _read_float_1D(gpath+"make_member_list_n.in."+str(case))
        out = list(map(list,zip(*_read_float_2D(gpath+"make_member_list_p.out."+str(case)))))
        self.assertTrue(np.allclose(_make_member_list(n),out))

    def test_8(self):
        from phenum.grouptheory import _make_member_list
        case = 8
        n = _read_float_1D(gpath+"make_member_list_n.in."+str(case))
        out = list(map(list,zip(*_read_float_2D(gpath+"make_member_list_p.out."+str(case)))))
        self.assertTrue(np.allclose(_make_member_list(n),out))

    def test_9(self):
        from phenum.grouptheory import _make_member_list
        case = 9
        n = _read_float_1D(gpath+"make_member_list_n.in."+str(case))
        out = list(map(list,zip(*_read_float_2D(gpath+"make_member_list_p.out."+str(case)))))
        self.assertTrue(np.allclose(_make_member_list(n),out))

    def test_10(self):
        from phenum.grouptheory import _make_member_list
        case = 10
        n = _read_float_1D(gpath+"make_member_list_n.in."+str(case))
        out = list(map(list,zip(*_read_float_2D(gpath+"make_member_list_p.out."+str(case)))))
        self.assertTrue(np.allclose(_make_member_list(n),out))

    def test_11(self):
        from phenum.grouptheory import _make_member_list
        case = 11
        n = _read_float_1D(gpath+"make_member_list_n.in."+str(case))
        out = list(map(list,zip(*_read_float_2D(gpath+"make_member_list_p.out."+str(case)))))
        self.assertTrue(np.allclose(_make_member_list(n),out))

    def test_12(self):
        from phenum.grouptheory import _make_member_list
        case = 12
        n = _read_float_1D(gpath+"make_member_list_n.in."+str(case))
        out = list(map(list,zip(*_read_float_2D(gpath+"make_member_list_p.out."+str(case)))))
        self.assertTrue(np.allclose(_make_member_list(n),out))

    def test_13(self):
        from phenum.grouptheory import _make_member_list
        case = 13
        n = _read_float_1D(gpath+"make_member_list_n.in."+str(case))
        out = list(map(list,zip(*_read_float_2D(gpath+"make_member_list_p.out."+str(case)))))
        self.assertTrue(np.allclose(_make_member_list(n),out))

    def test_14(self):
        from phenum.grouptheory import _make_member_list
        case = 14
        n = _read_float_1D(gpath+"make_member_list_n.in."+str(case))
        out = list(map(list,zip(*_read_float_2D(gpath+"make_member_list_p.out."+str(case)))))
        self.assertTrue(np.allclose(_make_member_list(n),out))

    def test_15(self):
        from phenum.grouptheory import _make_member_list
        case = 15
        n = _read_float_1D(gpath+"make_member_list_n.in."+str(case))
        out = list(map(list,zip(*_read_float_2D(gpath+"make_member_list_p.out."+str(case)))))
        self.assertTrue(np.allclose(_make_member_list(n),out))

    def test_16(self):
        from phenum.grouptheory import _make_member_list
        case = 16
        n = _read_float_1D(gpath+"make_member_list_n.in."+str(case))
        out = list(map(list,zip(*_read_float_2D(gpath+"make_member_list_p.out."+str(case)))))
        self.assertTrue(np.allclose(_make_member_list(n),out))

    def test_17(self):
        from phenum.grouptheory import _make_member_list
        case = 17
        n = _read_float_1D(gpath+"make_member_list_n.in."+str(case))
        out = list(map(list,zip(*_read_float_2D(gpath+"make_member_list_p.out."+str(case)))))
        self.assertTrue(np.allclose(_make_member_list(n),out))

    def test_18(self):
        from phenum.grouptheory import _make_member_list
        case = 18
        n = _read_float_1D(gpath+"make_member_list_n.in."+str(case))
        out = list(map(list,zip(*_read_float_2D(gpath+"make_member_list_p.out."+str(case)))))
        self.assertTrue(np.allclose(_make_member_list(n),out))

    def test_19(self):
        from phenum.grouptheory import _make_member_list
        case = 19
        n = _read_float_1D(gpath+"make_member_list_n.in."+str(case))
        out = list(map(list,zip(*_read_float_2D(gpath+"make_member_list_p.out."+str(case)))))
        self.assertTrue(np.allclose(_make_member_list(n),out))

    def test_20(self):
        from phenum.grouptheory import _make_member_list
        case = 20
        n = _read_float_1D(gpath+"make_member_list_n.in."+str(case))
        out = list(map(list,zip(*_read_float_2D(gpath+"make_member_list_p.out."+str(case)))))
        self.assertTrue(np.allclose(_make_member_list(n),out))

class TestFindPermutationOfGroup(ut.TestCase):
    """Tests of the _find_permutation_of_group subroutine."""

    def test_1(self):
        from phenum.grouptheory import _find_permutation_of_group
        case = 1
        g = list(map(list,zip(*_read_float_2D(gpath+"find_permutation_of_group_g.in."+str(case)))))
        gp = list(map(list,zip(*_read_float_2D(gpath+"find_permutation_of_group_gp.in."+str(case)))))
        out = [i-1 for i in _read_int_1D(gpath+"find_permutation_of_group_perm.out."+str(case))]
        self.assertEqual(_find_permutation_of_group(g,gp),out)

    def test_2(self):
        from phenum.grouptheory import _find_permutation_of_group
        case = 2
        g = list(map(list,zip(*_read_float_2D(gpath+"find_permutation_of_group_g.in."+str(case)))))
        gp = list(map(list,zip(*_read_float_2D(gpath+"find_permutation_of_group_gp.in."+str(case)))))
        out = [i-1 for i in _read_int_1D(gpath+"find_permutation_of_group_perm.out."+str(case))]
        self.assertEqual(_find_permutation_of_group(g,gp),out)

    def test_3(self):
        from phenum.grouptheory import _find_permutation_of_group
        case = 3
        g = list(map(list,zip(*_read_float_2D(gpath+"find_permutation_of_group_g.in."+str(case)))))
        gp = list(map(list,zip(*_read_float_2D(gpath+"find_permutation_of_group_gp.in."+str(case)))))
        out = [i-1 for i in _read_int_1D(gpath+"find_permutation_of_group_perm.out."+str(case))]
        self.assertEqual(_find_permutation_of_group(g,gp),out)

    def test_4(self):
        from phenum.grouptheory import _find_permutation_of_group
        case = 4
        g = list(map(list,zip(*_read_float_2D(gpath+"find_permutation_of_group_g.in."+str(case)))))
        gp = list(map(list,zip(*_read_float_2D(gpath+"find_permutation_of_group_gp.in."+str(case)))))
        out = [i-1 for i in _read_int_1D(gpath+"find_permutation_of_group_perm.out."+str(case))]
        self.assertEqual(_find_permutation_of_group(g,gp),out)

    def test_5(self):
        from phenum.grouptheory import _find_permutation_of_group
        case = 5
        g = list(map(list,zip(*_read_float_2D(gpath+"find_permutation_of_group_g.in."+str(case)))))
        gp = list(map(list,zip(*_read_float_2D(gpath+"find_permutation_of_group_gp.in."+str(case)))))
        out = [i-1 for i in _read_int_1D(gpath+"find_permutation_of_group_perm.out."+str(case))]
        self.assertEqual(_find_permutation_of_group(g,gp),out)

    def test_6(self):
        from phenum.grouptheory import _find_permutation_of_group
        case = 6
        g = list(map(list,zip(*_read_float_2D(gpath+"find_permutation_of_group_g.in."+str(case)))))
        gp = list(map(list,zip(*_read_float_2D(gpath+"find_permutation_of_group_gp.in."+str(case)))))
        out = [i-1 for i in _read_int_1D(gpath+"find_permutation_of_group_perm.out."+str(case))]
        self.assertEqual(_find_permutation_of_group(g,gp),out)

    def test_7(self):
        from phenum.grouptheory import _find_permutation_of_group
        case = 7
        g = list(map(list,zip(*_read_float_2D(gpath+"find_permutation_of_group_g.in."+str(case)))))
        gp = list(map(list,zip(*_read_float_2D(gpath+"find_permutation_of_group_gp.in."+str(case)))))
        out = [i-1 for i in _read_int_1D(gpath+"find_permutation_of_group_perm.out."+str(case))]
        self.assertEqual(_find_permutation_of_group(g,gp),out)

    def test_8(self):
        from phenum.grouptheory import _find_permutation_of_group
        case = 8
        g = list(map(list,zip(*_read_float_2D(gpath+"find_permutation_of_group_g.in."+str(case)))))
        gp = list(map(list,zip(*_read_float_2D(gpath+"find_permutation_of_group_gp.in."+str(case)))))
        out = [i-1 for i in _read_int_1D(gpath+"find_permutation_of_group_perm.out."+str(case))]
        self.assertEqual(_find_permutation_of_group(g,gp),out)

    def test_9(self):
        from phenum.grouptheory import _find_permutation_of_group
        case = 9
        g = list(map(list,zip(*_read_float_2D(gpath+"find_permutation_of_group_g.in."+str(case)))))
        gp = list(map(list,zip(*_read_float_2D(gpath+"find_permutation_of_group_gp.in."+str(case)))))
        out = [i-1 for i in _read_int_1D(gpath+"find_permutation_of_group_perm.out."+str(case))]
        self.assertEqual(_find_permutation_of_group(g,gp),out)

    def test_10(self):
        from phenum.grouptheory import _find_permutation_of_group
        case = 10
        g = list(map(list,zip(*_read_float_2D(gpath+"find_permutation_of_group_g.in."+str(case)))))
        gp = list(map(list,zip(*_read_float_2D(gpath+"find_permutation_of_group_gp.in."+str(case)))))
        out = [i-1 for i in _read_int_1D(gpath+"find_permutation_of_group_perm.out."+str(case))]
        self.assertEqual(_find_permutation_of_group(g,gp),out)

    def test_11(self):
        from phenum.grouptheory import _find_permutation_of_group
        case = 11
        g = list(map(list,zip(*_read_float_2D(gpath+"find_permutation_of_group_g.in."+str(case)))))
        gp = list(map(list,zip(*_read_float_2D(gpath+"find_permutation_of_group_gp.in."+str(case)))))
        out = [i-1 for i in _read_int_1D(gpath+"find_permutation_of_group_perm.out."+str(case))]
        self.assertEqual(_find_permutation_of_group(g,gp),out)

    def test_12(self):
        from phenum.grouptheory import _find_permutation_of_group
        case = 12
        g = list(map(list,zip(*_read_float_2D(gpath+"find_permutation_of_group_g.in."+str(case)))))
        gp = list(map(list,zip(*_read_float_2D(gpath+"find_permutation_of_group_gp.in."+str(case)))))
        out = [i-1 for i in _read_int_1D(gpath+"find_permutation_of_group_perm.out."+str(case))]
        self.assertEqual(_find_permutation_of_group(g,gp),out)

    def test_13(self):
        from phenum.grouptheory import _find_permutation_of_group
        case = 13
        g = list(map(list,zip(*_read_float_2D(gpath+"find_permutation_of_group_g.in."+str(case)))))
        gp = list(map(list,zip(*_read_float_2D(gpath+"find_permutation_of_group_gp.in."+str(case)))))
        out = [i-1 for i in _read_int_1D(gpath+"find_permutation_of_group_perm.out."+str(case))]
        self.assertEqual(_find_permutation_of_group(g,gp),out)

    def test_14(self):
        from phenum.grouptheory import _find_permutation_of_group
        case = 14
        g = list(map(list,zip(*_read_float_2D(gpath+"find_permutation_of_group_g.in."+str(case)))))
        gp = list(map(list,zip(*_read_float_2D(gpath+"find_permutation_of_group_gp.in."+str(case)))))
        out = [i-1 for i in _read_int_1D(gpath+"find_permutation_of_group_perm.out."+str(case))]
        self.assertEqual(_find_permutation_of_group(g,gp),out)

    def test_15(self):
        from phenum.grouptheory import _find_permutation_of_group
        case = 15
        g = list(map(list,zip(*_read_float_2D(gpath+"find_permutation_of_group_g.in."+str(case)))))
        gp = list(map(list,zip(*_read_float_2D(gpath+"find_permutation_of_group_gp.in."+str(case)))))
        out = [i-1 for i in _read_int_1D(gpath+"find_permutation_of_group_perm.out."+str(case))]
        self.assertEqual(_find_permutation_of_group(g,gp),out)

    def test_16(self):
        from phenum.grouptheory import _find_permutation_of_group
        case = 16
        g = list(map(list,zip(*_read_float_2D(gpath+"find_permutation_of_group_g.in."+str(case)))))
        gp = list(map(list,zip(*_read_float_2D(gpath+"find_permutation_of_group_gp.in."+str(case)))))
        out = [i-1 for i in _read_int_1D(gpath+"find_permutation_of_group_perm.out."+str(case))]
        self.assertEqual(_find_permutation_of_group(g,gp),out)

    def test_17(self):
        from phenum.grouptheory import _find_permutation_of_group
        case = 17
        g = list(map(list,zip(*_read_float_2D(gpath+"find_permutation_of_group_g.in."+str(case)))))
        gp = list(map(list,zip(*_read_float_2D(gpath+"find_permutation_of_group_gp.in."+str(case)))))
        out = [i-1 for i in _read_int_1D(gpath+"find_permutation_of_group_perm.out."+str(case))]
        self.assertEqual(_find_permutation_of_group(g,gp),out)

    def test_18(self):
        from phenum.grouptheory import _find_permutation_of_group
        case = 18
        g = list(map(list,zip(*_read_float_2D(gpath+"find_permutation_of_group_g.in."+str(case)))))
        gp = list(map(list,zip(*_read_float_2D(gpath+"find_permutation_of_group_gp.in."+str(case)))))
        out = [i-1 for i in _read_int_1D(gpath+"find_permutation_of_group_perm.out."+str(case))]
        self.assertEqual(_find_permutation_of_group(g,gp),out)

    def test_19(self):
        from phenum.grouptheory import _find_permutation_of_group
        case = 19
        g = list(map(list,zip(*_read_float_2D(gpath+"find_permutation_of_group_g.in."+str(case)))))
        gp = list(map(list,zip(*_read_float_2D(gpath+"find_permutation_of_group_gp.in."+str(case)))))
        out = [i-1 for i in _read_int_1D(gpath+"find_permutation_of_group_perm.out."+str(case))]
        self.assertEqual(_find_permutation_of_group(g,gp),out)

    def test_20(self):
        from phenum.grouptheory import _find_permutation_of_group
        case = 20
        g = list(map(list,zip(*_read_float_2D(gpath+"find_permutation_of_group_g.in."+str(case)))))
        gp = list(map(list,zip(*_read_float_2D(gpath+"find_permutation_of_group_gp.in."+str(case)))))
        out = [i-1 for i in _read_int_1D(gpath+"find_permutation_of_group_perm.out."+str(case))]
        self.assertEqual(_find_permutation_of_group(g,gp),out)

class TestIsEquivLattice(ut.TestCase):
    """Tests of the _is_equiv_lattice subroutine."""

    def test_1(self):
        from phenum.grouptheory import _is_equiv_lattice
        case = 1
        lat1 = _read_float_2D(gpath+"is_equiv_lattice_lat1.in."+str(case))
        lat2 = _read_float_2D(gpath+"is_equiv_lattice_lat2.in."+str(case))
        eps = _read_float(gpath+"is_equiv_lattice_eps.in."+str(case))
        out = _read_logical(gpath+"is_equiv_lattice.out."+str(case))
        self.assertEqual(_is_equiv_lattice(lat1,lat2,eps),out)

    def test_2(self):
        from phenum.grouptheory import _is_equiv_lattice
        case = 2
        lat1 = _read_float_2D(gpath+"is_equiv_lattice_lat1.in."+str(case))
        lat2 = _read_float_2D(gpath+"is_equiv_lattice_lat2.in."+str(case))
        eps = _read_float(gpath+"is_equiv_lattice_eps.in."+str(case))
        out = _read_logical(gpath+"is_equiv_lattice.out."+str(case))
        self.assertEqual(_is_equiv_lattice(lat1,lat2,eps),out)

    def test_3(self):
        from phenum.grouptheory import _is_equiv_lattice
        case = 3
        lat1 = _read_float_2D(gpath+"is_equiv_lattice_lat1.in."+str(case))
        lat2 = _read_float_2D(gpath+"is_equiv_lattice_lat2.in."+str(case))
        eps = _read_float(gpath+"is_equiv_lattice_eps.in."+str(case))
        out = _read_logical(gpath+"is_equiv_lattice.out."+str(case))
        self.assertEqual(_is_equiv_lattice(lat1,lat2,eps),out)

    def test_4(self):
        from phenum.grouptheory import _is_equiv_lattice
        case = 4
        lat1 = _read_float_2D(gpath+"is_equiv_lattice_lat1.in."+str(case))
        lat2 = _read_float_2D(gpath+"is_equiv_lattice_lat2.in."+str(case))
        eps = _read_float(gpath+"is_equiv_lattice_eps.in."+str(case))
        out = _read_logical(gpath+"is_equiv_lattice.out."+str(case))
        self.assertEqual(_is_equiv_lattice(lat1,lat2,eps),out)

    def test_5(self):
        from phenum.grouptheory import _is_equiv_lattice
        case = 5
        lat1 = _read_float_2D(gpath+"is_equiv_lattice_lat1.in."+str(case))
        lat2 = _read_float_2D(gpath+"is_equiv_lattice_lat2.in."+str(case))
        eps = _read_float(gpath+"is_equiv_lattice_eps.in."+str(case))
        out = _read_logical(gpath+"is_equiv_lattice.out."+str(case))
        self.assertEqual(_is_equiv_lattice(lat1,lat2,eps),out)

    def test_6(self):
        from phenum.grouptheory import _is_equiv_lattice
        case = 6
        lat1 = _read_float_2D(gpath+"is_equiv_lattice_lat1.in."+str(case))
        lat2 = _read_float_2D(gpath+"is_equiv_lattice_lat2.in."+str(case))
        eps = _read_float(gpath+"is_equiv_lattice_eps.in."+str(case))
        out = _read_logical(gpath+"is_equiv_lattice.out."+str(case))
        self.assertEqual(_is_equiv_lattice(lat1,lat2,eps),out)

    def test_7(self):
        from phenum.grouptheory import _is_equiv_lattice
        case = 7
        lat1 = _read_float_2D(gpath+"is_equiv_lattice_lat1.in."+str(case))
        lat2 = _read_float_2D(gpath+"is_equiv_lattice_lat2.in."+str(case))
        eps = _read_float(gpath+"is_equiv_lattice_eps.in."+str(case))
        out = _read_logical(gpath+"is_equiv_lattice.out."+str(case))
        self.assertEqual(_is_equiv_lattice(lat1,lat2,eps),out)

    def test_8(self):
        from phenum.grouptheory import _is_equiv_lattice
        case = 8
        lat1 = _read_float_2D(gpath+"is_equiv_lattice_lat1.in."+str(case))
        lat2 = _read_float_2D(gpath+"is_equiv_lattice_lat2.in."+str(case))
        eps = _read_float(gpath+"is_equiv_lattice_eps.in."+str(case))
        out = _read_logical(gpath+"is_equiv_lattice.out."+str(case))
        self.assertEqual(_is_equiv_lattice(lat1,lat2,eps),out)

    def test_9(self):
        from phenum.grouptheory import _is_equiv_lattice
        case = 9
        lat1 = _read_float_2D(gpath+"is_equiv_lattice_lat1.in."+str(case))
        lat2 = _read_float_2D(gpath+"is_equiv_lattice_lat2.in."+str(case))
        eps = _read_float(gpath+"is_equiv_lattice_eps.in."+str(case))
        out = _read_logical(gpath+"is_equiv_lattice.out."+str(case))
        self.assertEqual(_is_equiv_lattice(lat1,lat2,eps),out)

    def test_10(self):
        from phenum.grouptheory import _is_equiv_lattice
        case = 10
        lat1 = _read_float_2D(gpath+"is_equiv_lattice_lat1.in."+str(case))
        lat2 = _read_float_2D(gpath+"is_equiv_lattice_lat2.in."+str(case))
        eps = _read_float(gpath+"is_equiv_lattice_eps.in."+str(case))
        out = _read_logical(gpath+"is_equiv_lattice.out."+str(case))
        self.assertEqual(_is_equiv_lattice(lat1,lat2,eps),out)

    def test_11(self):
        from phenum.grouptheory import _is_equiv_lattice
        case = 11
        lat1 = _read_float_2D(gpath+"is_equiv_lattice_lat1.in."+str(case))
        lat2 = _read_float_2D(gpath+"is_equiv_lattice_lat2.in."+str(case))
        eps = _read_float(gpath+"is_equiv_lattice_eps.in."+str(case))
        out = _read_logical(gpath+"is_equiv_lattice.out."+str(case))
        self.assertEqual(_is_equiv_lattice(lat1,lat2,eps),out)

    def test_12(self):
        from phenum.grouptheory import _is_equiv_lattice
        case = 12
        lat1 = _read_float_2D(gpath+"is_equiv_lattice_lat1.in."+str(case))
        lat2 = _read_float_2D(gpath+"is_equiv_lattice_lat2.in."+str(case))
        eps = _read_float(gpath+"is_equiv_lattice_eps.in."+str(case))
        out = _read_logical(gpath+"is_equiv_lattice.out."+str(case))
        self.assertEqual(_is_equiv_lattice(lat1,lat2,eps),out)

    def test_13(self):
        from phenum.grouptheory import _is_equiv_lattice
        case = 13
        lat1 = _read_float_2D(gpath+"is_equiv_lattice_lat1.in."+str(case))
        lat2 = _read_float_2D(gpath+"is_equiv_lattice_lat2.in."+str(case))
        eps = _read_float(gpath+"is_equiv_lattice_eps.in."+str(case))
        out = _read_logical(gpath+"is_equiv_lattice.out."+str(case))
        self.assertEqual(_is_equiv_lattice(lat1,lat2,eps),out)

    def test_14(self):
        from phenum.grouptheory import _is_equiv_lattice
        case = 14
        lat1 = _read_float_2D(gpath+"is_equiv_lattice_lat1.in."+str(case))
        lat2 = _read_float_2D(gpath+"is_equiv_lattice_lat2.in."+str(case))
        eps = _read_float(gpath+"is_equiv_lattice_eps.in."+str(case))
        out = _read_logical(gpath+"is_equiv_lattice.out."+str(case))
        self.assertEqual(_is_equiv_lattice(lat1,lat2,eps),out)

    def test_15(self):
        from phenum.grouptheory import _is_equiv_lattice
        case = 15
        lat1 = _read_float_2D(gpath+"is_equiv_lattice_lat1.in."+str(case))
        lat2 = _read_float_2D(gpath+"is_equiv_lattice_lat2.in."+str(case))
        eps = _read_float(gpath+"is_equiv_lattice_eps.in."+str(case))
        out = _read_logical(gpath+"is_equiv_lattice.out."+str(case))
        self.assertEqual(_is_equiv_lattice(lat1,lat2,eps),out)

    def test_16(self):
        from phenum.grouptheory import _is_equiv_lattice
        case = 16
        lat1 = _read_float_2D(gpath+"is_equiv_lattice_lat1.in."+str(case))
        lat2 = _read_float_2D(gpath+"is_equiv_lattice_lat2.in."+str(case))
        eps = _read_float(gpath+"is_equiv_lattice_eps.in."+str(case))
        out = _read_logical(gpath+"is_equiv_lattice.out."+str(case))
        self.assertEqual(_is_equiv_lattice(lat1,lat2,eps),out)

    def test_17(self):
        from phenum.grouptheory import _is_equiv_lattice
        case = 17
        lat1 = _read_float_2D(gpath+"is_equiv_lattice_lat1.in."+str(case))
        lat2 = _read_float_2D(gpath+"is_equiv_lattice_lat2.in."+str(case))
        eps = _read_float(gpath+"is_equiv_lattice_eps.in."+str(case))
        out = _read_logical(gpath+"is_equiv_lattice.out."+str(case))
        self.assertEqual(_is_equiv_lattice(lat1,lat2,eps),out)

    def test_18(self):
        from phenum.grouptheory import _is_equiv_lattice
        case = 18
        lat1 = _read_float_2D(gpath+"is_equiv_lattice_lat1.in."+str(case))
        lat2 = _read_float_2D(gpath+"is_equiv_lattice_lat2.in."+str(case))
        eps = _read_float(gpath+"is_equiv_lattice_eps.in."+str(case))
        out = _read_logical(gpath+"is_equiv_lattice.out."+str(case))
        self.assertEqual(_is_equiv_lattice(lat1,lat2,eps),out)

    def test_19(self):
        from phenum.grouptheory import _is_equiv_lattice
        case = 19
        lat1 = _read_float_2D(gpath+"is_equiv_lattice_lat1.in."+str(case))
        lat2 = _read_float_2D(gpath+"is_equiv_lattice_lat2.in."+str(case))
        eps = _read_float(gpath+"is_equiv_lattice_eps.in."+str(case))
        out = _read_logical(gpath+"is_equiv_lattice.out."+str(case))
        self.assertEqual(_is_equiv_lattice(lat1,lat2,eps),out)

    def test_20(self):
        from phenum.grouptheory import _is_equiv_lattice
        case = 20
        lat1 = _read_float_2D(gpath+"is_equiv_lattice_lat1.in."+str(case))
        lat2 = _read_float_2D(gpath+"is_equiv_lattice_lat2.in."+str(case))
        eps = _read_float(gpath+"is_equiv_lattice_eps.in."+str(case))
        out = _read_logical(gpath+"is_equiv_lattice.out."+str(case))
        self.assertEqual(_is_equiv_lattice(lat1,lat2,eps),out)
        
class TestGetSLVFixingOperations(ut.TestCase):
    """Tests of the _get_sLV_fixing_operations subroutine."""

    def _compare_outputs(self,out1,out2):
        fix1 = out1[0]
        fix2 = out2[0]
        rot1 = out1[1]
        rot2 = out2[1]
        deg1  = out1[2]
        deg2 = out2[2]

        self.assertEqual(deg1,deg2)

        if len(fix1.rot) == len(fix2.rot):
            for i in range(len(fix1.rot)):
                for j in range(3):
                    for k in range(3):
                        self.assertAlmostEqual(fix1.rot[i][j][k],fix2.rot[i][j][k],places=12)
        else:
            self.assertEqual(len(fix1.rot),len(fix2.rot))

        if len(fix1.shift) == len(fix2.shift):
            for i in range(len(fix1.shift)):
                for j in range(3):
                    self.assertAlmostEqual(fix1.shift[i][j],fix2.shift[i][j],places=12)
        else:
            self.assertEqual(len(fix1.shift),len(fix2.shift))

        self.assertEqual(rot1.nL,rot2.nL)
        self.assertEqual(rot1.RotIndx, rot2.RotIndx)

        if (rot1.v != None) and (rot2.v != None):
            if len(rot1.v) == len(rot2.v):
                for i in range(len(rot1.v)):
                    for j in range(len(rot1.v[i])):
                        for k in range(len(rot1.v[i][j])):
                            self.assertAlmostEqual(rot1.v[i][j][k],rot2.v[i][j][k],places=12)
            else:
                self.assertEqual(len(rot1.v),len(rot2.v))
        else:
            self.assertEqual(rot1.v,rot2.v)

        # if (rot1.perm.site_perm != None) and (rot2.perm.site_perm != None):
        #     if len(rot1.perm.site_perm) == len(rot2.perm.site_perm):
        #         for i in range(len(rot1.perm.site_perm)):
        #             for j in range(len(rot1.perm.site_perm[i])):
        #                 self.assertEqual(rot1.perm.site_perm[i][j],rot2.perm.site_perm[i][j])
        #     else:
        #         self.assertEqual(len(rot1.perm.site_perm),len(rot2.perm.site_perm))
        # else:
        #     self.assertEqual(rot1.perm.site_perm,rot2.perm.site_perm)
            
        if (rot1.perm.arrow_perm != None) and (rot2.perm.arrow_perm != None):
            if len(rot1.perm.arrow_perm) == len(rot2.perm.arrow_perm):
                for i in range(len(rot1.perm.arrow_perm)):
                    for j in range(len(rot1.perm.arrow_perm[i])):
                        self.assertEqual(rot1.perm.arrow_perm[i][j],rot2.perm.arrow_perm[i][j])
            else:
                self.assertEqual(len(rot1.perm.arrow_perm),len(rot2.perm.arrow_perm))
        else:
            self.assertEqual(rot1.perm.arrow_perm,rot2.perm.arrow_perm)
        
    def test_1(self):
        from phenum.grouptheory import _get_sLV_fixing_operations
        case = 1
        HNF = _read_int_2D(gpath+"get_sLV_fixing_operations_HNF.in."+str(case))
        pLV = np.transpose(_read_float_2D(gpath+"get_sLV_fixing_operations_pLV.in."+str(case)))
        nD = _read_int(gpath+"get_sLV_fixing_operations_nD.in."+str(case))
        rot = _read_float_3D(gpath+"get_sLV_fixing_operations_rot.in."+str(case))
        shift = list(map(list,zip(*_read_float_2D(gpath+"get_sLV_fixing_operations_shift.in."+str(case)))))
        eps = _read_float(gpath+"get_sLV_fixing_operations_eps.in."+str(case))
        dPerm = _read_RotPermList(gpath+"get_sLV_fixing_operations_dPerm.in."+str(case))
        fixOp, rotPerm, degeneracy = _get_sLV_fixing_operations(HNF,pLV,rot,shift,dPerm,eps)
        rotPerm_out = _read_RotPermList(gpath+"get_sLV_fixing_operations_rotPerm.out."+str(case))
        degen_out = _read_int(gpath+"get_sLV_fixing_operations_degeneracy.out."+str(case))
        fixOp_out = _read_fixOp(gpath+"get_sLV_fixing_operations_fixOp.out."+str(case))

        self._compare_outputs([fixOp,rotPerm,degeneracy],[fixOp_out,rotPerm_out,degen_out])
        
    def test_2(self):
        from phenum.grouptheory import _get_sLV_fixing_operations
        case = 10
        HNF = _read_int_2D(gpath+"get_sLV_fixing_operations_HNF.in."+str(case))
        pLV = np.transpose(_read_float_2D(gpath+"get_sLV_fixing_operations_pLV.in."+str(case)))
        nD = _read_int(gpath+"get_sLV_fixing_operations_nD.in."+str(case))
        rot = _read_float_3D(gpath+"get_sLV_fixing_operations_rot.in."+str(case))
        shift = list(map(list,zip(*_read_float_2D(gpath+"get_sLV_fixing_operations_shift.in."+str(case)))))
        eps = _read_float(gpath+"get_sLV_fixing_operations_eps.in."+str(case))
        dPerm = _read_RotPermList(gpath+"get_sLV_fixing_operations_dPerm.in."+str(case))
        fixOp, rotPerm, degeneracy = _get_sLV_fixing_operations(HNF,pLV,rot,shift,dPerm,eps)
        rotPerm_out = _read_RotPermList(gpath+"get_sLV_fixing_operations_rotPerm.out."+str(case))
        degen_out = _read_int(gpath+"get_sLV_fixing_operations_degeneracy.out."+str(case))
        fixOp_out = _read_fixOp(gpath+"get_sLV_fixing_operations_fixOp.out."+str(case))

        self._compare_outputs([fixOp,rotPerm,degeneracy],[fixOp_out,rotPerm_out,degen_out])
        
    def test_3(self):
        from phenum.grouptheory import _get_sLV_fixing_operations
        case = 20
        HNF = _read_int_2D(gpath+"get_sLV_fixing_operations_HNF.in."+str(case))
        pLV = np.transpose(_read_float_2D(gpath+"get_sLV_fixing_operations_pLV.in."+str(case)))
        nD = _read_int(gpath+"get_sLV_fixing_operations_nD.in."+str(case))
        rot = _read_float_3D(gpath+"get_sLV_fixing_operations_rot.in."+str(case))
        shift = list(map(list,zip(*_read_float_2D(gpath+"get_sLV_fixing_operations_shift.in."+str(case)))))
        eps = _read_float(gpath+"get_sLV_fixing_operations_eps.in."+str(case))
        dPerm = _read_RotPermList(gpath+"get_sLV_fixing_operations_dPerm.in."+str(case))
        fixOp, rotPerm, degeneracy = _get_sLV_fixing_operations(HNF,pLV,rot,shift,dPerm,eps)
        rotPerm_out = _read_RotPermList(gpath+"get_sLV_fixing_operations_rotPerm.out."+str(case))
        degen_out = _read_int(gpath+"get_sLV_fixing_operations_degeneracy.out."+str(case))
        fixOp_out = _read_fixOp(gpath+"get_sLV_fixing_operations_fixOp.out."+str(case))

        self._compare_outputs([fixOp,rotPerm,degeneracy],[fixOp_out,rotPerm_out,degen_out])
        
    def test_4(self):
        from phenum.grouptheory import _get_sLV_fixing_operations
        case = 30
        HNF = _read_int_2D(gpath+"get_sLV_fixing_operations_HNF.in."+str(case))
        pLV = np.transpose(_read_float_2D(gpath+"get_sLV_fixing_operations_pLV.in."+str(case)))
        nD = _read_int(gpath+"get_sLV_fixing_operations_nD.in."+str(case))
        rot = _read_float_3D(gpath+"get_sLV_fixing_operations_rot.in."+str(case))
        shift = list(map(list,zip(*_read_float_2D(gpath+"get_sLV_fixing_operations_shift.in."+str(case)))))
        eps = _read_float(gpath+"get_sLV_fixing_operations_eps.in."+str(case))
        dPerm = _read_RotPermList(gpath+"get_sLV_fixing_operations_dPerm.in."+str(case))
        fixOp, rotPerm, degeneracy = _get_sLV_fixing_operations(HNF,pLV,rot,shift,dPerm,eps)
        rotPerm_out = _read_RotPermList(gpath+"get_sLV_fixing_operations_rotPerm.out."+str(case))
        degen_out = _read_int(gpath+"get_sLV_fixing_operations_degeneracy.out."+str(case))
        fixOp_out = _read_fixOp(gpath+"get_sLV_fixing_operations_fixOp.out."+str(case))

        self._compare_outputs([fixOp,rotPerm,degeneracy],[fixOp_out,rotPerm_out,degen_out])
        
    def test_5(self):
        from phenum.grouptheory import _get_sLV_fixing_operations
        case = 40
        HNF = _read_int_2D(gpath+"get_sLV_fixing_operations_HNF.in."+str(case))
        pLV = np.transpose(_read_float_2D(gpath+"get_sLV_fixing_operations_pLV.in."+str(case)))
        nD = _read_int(gpath+"get_sLV_fixing_operations_nD.in."+str(case))
        rot = _read_float_3D(gpath+"get_sLV_fixing_operations_rot.in."+str(case))
        shift = list(map(list,zip(*_read_float_2D(gpath+"get_sLV_fixing_operations_shift.in."+str(case)))))
        eps = _read_float(gpath+"get_sLV_fixing_operations_eps.in."+str(case))
        dPerm = _read_RotPermList(gpath+"get_sLV_fixing_operations_dPerm.in."+str(case))
        fixOp, rotPerm, degeneracy = _get_sLV_fixing_operations(HNF,pLV,rot,shift,dPerm,eps)
        rotPerm_out = _read_RotPermList(gpath+"get_sLV_fixing_operations_rotPerm.out."+str(case))
        degen_out = _read_int(gpath+"get_sLV_fixing_operations_degeneracy.out."+str(case))
        fixOp_out = _read_fixOp(gpath+"get_sLV_fixing_operations_fixOp.out."+str(case))

        self._compare_outputs([fixOp,rotPerm,degeneracy],[fixOp_out,rotPerm_out,degen_out])
        
    def test_6(self):
        from phenum.grouptheory import _get_sLV_fixing_operations
        case = 50
        HNF = _read_int_2D(gpath+"get_sLV_fixing_operations_HNF.in."+str(case))
        pLV = np.transpose(_read_float_2D(gpath+"get_sLV_fixing_operations_pLV.in."+str(case)))
        nD = _read_int(gpath+"get_sLV_fixing_operations_nD.in."+str(case))
        rot = _read_float_3D(gpath+"get_sLV_fixing_operations_rot.in."+str(case))
        shift = list(map(list,zip(*_read_float_2D(gpath+"get_sLV_fixing_operations_shift.in."+str(case)))))
        eps = _read_float(gpath+"get_sLV_fixing_operations_eps.in."+str(case))
        dPerm = _read_RotPermList(gpath+"get_sLV_fixing_operations_dPerm.in."+str(case))
        fixOp, rotPerm, degeneracy = _get_sLV_fixing_operations(HNF,pLV,rot,shift,dPerm,eps)
        rotPerm_out = _read_RotPermList(gpath+"get_sLV_fixing_operations_rotPerm.out."+str(case))
        degen_out = _read_int(gpath+"get_sLV_fixing_operations_degeneracy.out."+str(case))
        fixOp_out = _read_fixOp(gpath+"get_sLV_fixing_operations_fixOp.out."+str(case))

        self._compare_outputs([fixOp,rotPerm,degeneracy],[fixOp_out,rotPerm_out,degen_out])
        
    def test_7(self):
        from phenum.grouptheory import _get_sLV_fixing_operations
        case = 60
        HNF = _read_int_2D(gpath+"get_sLV_fixing_operations_HNF.in."+str(case))
        pLV = np.transpose(_read_float_2D(gpath+"get_sLV_fixing_operations_pLV.in."+str(case)))
        nD = _read_int(gpath+"get_sLV_fixing_operations_nD.in."+str(case))
        rot = _read_float_3D(gpath+"get_sLV_fixing_operations_rot.in."+str(case))
        shift = list(map(list,zip(*_read_float_2D(gpath+"get_sLV_fixing_operations_shift.in."+str(case)))))
        eps = _read_float(gpath+"get_sLV_fixing_operations_eps.in."+str(case))
        dPerm = _read_RotPermList(gpath+"get_sLV_fixing_operations_dPerm.in."+str(case))
        fixOp, rotPerm, degeneracy = _get_sLV_fixing_operations(HNF,pLV,rot,shift,dPerm,eps)
        rotPerm_out = _read_RotPermList(gpath+"get_sLV_fixing_operations_rotPerm.out."+str(case))
        degen_out = _read_int(gpath+"get_sLV_fixing_operations_degeneracy.out."+str(case))
        fixOp_out = _read_fixOp(gpath+"get_sLV_fixing_operations_fixOp.out."+str(case))

        self._compare_outputs([fixOp,rotPerm,degeneracy],[fixOp_out,rotPerm_out,degen_out])
        
    def test_8(self):
        from phenum.grouptheory import _get_sLV_fixing_operations
        case = 70
        HNF = _read_int_2D(gpath+"get_sLV_fixing_operations_HNF.in."+str(case))
        pLV = np.transpose(_read_float_2D(gpath+"get_sLV_fixing_operations_pLV.in."+str(case)))
        nD = _read_int(gpath+"get_sLV_fixing_operations_nD.in."+str(case))
        rot = _read_float_3D(gpath+"get_sLV_fixing_operations_rot.in."+str(case))
        shift = list(map(list,zip(*_read_float_2D(gpath+"get_sLV_fixing_operations_shift.in."+str(case)))))
        eps = _read_float(gpath+"get_sLV_fixing_operations_eps.in."+str(case))
        dPerm = _read_RotPermList(gpath+"get_sLV_fixing_operations_dPerm.in."+str(case))
        fixOp, rotPerm, degeneracy = _get_sLV_fixing_operations(HNF,pLV,rot,shift,dPerm,eps)
        rotPerm_out = _read_RotPermList(gpath+"get_sLV_fixing_operations_rotPerm.out."+str(case))
        degen_out = _read_int(gpath+"get_sLV_fixing_operations_degeneracy.out."+str(case))
        fixOp_out = _read_fixOp(gpath+"get_sLV_fixing_operations_fixOp.out."+str(case))

        self._compare_outputs([fixOp,rotPerm,degeneracy],[fixOp_out,rotPerm_out,degen_out])
        
    def test_9(self):
        from phenum.grouptheory import _get_sLV_fixing_operations
        case = 80
        HNF = _read_int_2D(gpath+"get_sLV_fixing_operations_HNF.in."+str(case))
        pLV = np.transpose(_read_float_2D(gpath+"get_sLV_fixing_operations_pLV.in."+str(case)))
        nD = _read_int(gpath+"get_sLV_fixing_operations_nD.in."+str(case))
        rot = _read_float_3D(gpath+"get_sLV_fixing_operations_rot.in."+str(case))
        shift = list(map(list,zip(*_read_float_2D(gpath+"get_sLV_fixing_operations_shift.in."+str(case)))))
        eps = _read_float(gpath+"get_sLV_fixing_operations_eps.in."+str(case))
        dPerm = _read_RotPermList(gpath+"get_sLV_fixing_operations_dPerm.in."+str(case))
        fixOp, rotPerm, degeneracy = _get_sLV_fixing_operations(HNF,pLV,rot,shift,dPerm,eps)
        rotPerm_out = _read_RotPermList(gpath+"get_sLV_fixing_operations_rotPerm.out."+str(case))
        degen_out = _read_int(gpath+"get_sLV_fixing_operations_degeneracy.out."+str(case))
        fixOp_out = _read_fixOp(gpath+"get_sLV_fixing_operations_fixOp.out."+str(case))

        self._compare_outputs([fixOp,rotPerm,degeneracy],[fixOp_out,rotPerm_out,degen_out])
        
    def test_10(self):
        from phenum.grouptheory import _get_sLV_fixing_operations
        case = 90
        HNF = _read_int_2D(gpath+"get_sLV_fixing_operations_HNF.in."+str(case))
        pLV = np.transpose(_read_float_2D(gpath+"get_sLV_fixing_operations_pLV.in."+str(case)))
        nD = _read_int(gpath+"get_sLV_fixing_operations_nD.in."+str(case))
        rot = _read_float_3D(gpath+"get_sLV_fixing_operations_rot.in."+str(case))
        shift = list(map(list,zip(*_read_float_2D(gpath+"get_sLV_fixing_operations_shift.in."+str(case)))))
        eps = _read_float(gpath+"get_sLV_fixing_operations_eps.in."+str(case))
        dPerm = _read_RotPermList(gpath+"get_sLV_fixing_operations_dPerm.in."+str(case))
        fixOp, rotPerm, degeneracy = _get_sLV_fixing_operations(HNF,pLV,rot,shift,dPerm,eps)
        rotPerm_out = _read_RotPermList(gpath+"get_sLV_fixing_operations_rotPerm.out."+str(case))
        degen_out = _read_int(gpath+"get_sLV_fixing_operations_degeneracy.out."+str(case))
        fixOp_out = _read_fixOp(gpath+"get_sLV_fixing_operations_fixOp.out."+str(case))
        self._compare_outputs([fixOp,rotPerm,degeneracy],[fixOp_out,rotPerm_out,degen_out])
        
    def test_11(self):
        from phenum.grouptheory import _get_sLV_fixing_operations
        case = 100
        HNF = _read_int_2D(gpath+"get_sLV_fixing_operations_HNF.in."+str(case))
        pLV = np.transpose(_read_float_2D(gpath+"get_sLV_fixing_operations_pLV.in."+str(case)))
        nD = _read_int(gpath+"get_sLV_fixing_operations_nD.in."+str(case))
        rot = _read_float_3D(gpath+"get_sLV_fixing_operations_rot.in."+str(case))
        shift = list(map(list,zip(*_read_float_2D(gpath+"get_sLV_fixing_operations_shift.in."+str(case)))))
        eps = _read_float(gpath+"get_sLV_fixing_operations_eps.in."+str(case))
        dPerm = _read_RotPermList(gpath+"get_sLV_fixing_operations_dPerm.in."+str(case))
        fixOp, rotPerm, degeneracy = _get_sLV_fixing_operations(HNF,pLV,rot,shift,dPerm,eps)
        rotPerm_out = _read_RotPermList(gpath+"get_sLV_fixing_operations_rotPerm.out."+str(case))
        degen_out = _read_int(gpath+"get_sLV_fixing_operations_degeneracy.out."+str(case))
        fixOp_out = _read_fixOp(gpath+"get_sLV_fixing_operations_fixOp.out."+str(case))

        self._compare_outputs([fixOp,rotPerm,degeneracy],[fixOp_out,rotPerm_out,degen_out])
        
    def test_12(self):
        from phenum.grouptheory import _get_sLV_fixing_operations
        case = 110
        HNF = _read_int_2D(gpath+"get_sLV_fixing_operations_HNF.in."+str(case))
        pLV = np.transpose(_read_float_2D(gpath+"get_sLV_fixing_operations_pLV.in."+str(case)))
        nD = _read_int(gpath+"get_sLV_fixing_operations_nD.in."+str(case))
        rot = _read_float_3D(gpath+"get_sLV_fixing_operations_rot.in."+str(case))
        shift = list(map(list,zip(*_read_float_2D(gpath+"get_sLV_fixing_operations_shift.in."+str(case)))))
        eps = _read_float(gpath+"get_sLV_fixing_operations_eps.in."+str(case))
        dPerm = _read_RotPermList(gpath+"get_sLV_fixing_operations_dPerm.in."+str(case))
        fixOp, rotPerm, degeneracy = _get_sLV_fixing_operations(HNF,pLV,rot,shift,dPerm,eps)
        rotPerm_out = _read_RotPermList(gpath+"get_sLV_fixing_operations_rotPerm.out."+str(case))
        degen_out = _read_int(gpath+"get_sLV_fixing_operations_degeneracy.out."+str(case))
        fixOp_out = _read_fixOp(gpath+"get_sLV_fixing_operations_fixOp.out."+str(case))

        self._compare_outputs([fixOp,rotPerm,degeneracy],[fixOp_out,rotPerm_out,degen_out])
        
    def test_13(self):
        from phenum.grouptheory import _get_sLV_fixing_operations
        case = 120
        HNF = _read_int_2D(gpath+"get_sLV_fixing_operations_HNF.in."+str(case))
        pLV = np.transpose(_read_float_2D(gpath+"get_sLV_fixing_operations_pLV.in."+str(case)))
        nD = _read_int(gpath+"get_sLV_fixing_operations_nD.in."+str(case))
        rot = _read_float_3D(gpath+"get_sLV_fixing_operations_rot.in."+str(case))
        shift = list(map(list,zip(*_read_float_2D(gpath+"get_sLV_fixing_operations_shift.in."+str(case)))))
        eps = _read_float(gpath+"get_sLV_fixing_operations_eps.in."+str(case))
        dPerm = _read_RotPermList(gpath+"get_sLV_fixing_operations_dPerm.in."+str(case))
        fixOp, rotPerm, degeneracy = _get_sLV_fixing_operations(HNF,pLV,rot,shift,dPerm,eps)
        rotPerm_out = _read_RotPermList(gpath+"get_sLV_fixing_operations_rotPerm.out."+str(case))
        degen_out = _read_int(gpath+"get_sLV_fixing_operations_degeneracy.out."+str(case))
        fixOp_out = _read_fixOp(gpath+"get_sLV_fixing_operations_fixOp.out."+str(case))

        self._compare_outputs([fixOp,rotPerm,degeneracy],[fixOp_out,rotPerm_out,degen_out])
        
    def test_14(self):
        from phenum.grouptheory import _get_sLV_fixing_operations
        case = 130
        HNF = _read_int_2D(gpath+"get_sLV_fixing_operations_HNF.in."+str(case))
        pLV = np.transpose(_read_float_2D(gpath+"get_sLV_fixing_operations_pLV.in."+str(case)))
        nD = _read_int(gpath+"get_sLV_fixing_operations_nD.in."+str(case))
        rot = _read_float_3D(gpath+"get_sLV_fixing_operations_rot.in."+str(case))
        shift = list(map(list,zip(*_read_float_2D(gpath+"get_sLV_fixing_operations_shift.in."+str(case)))))
        eps = _read_float(gpath+"get_sLV_fixing_operations_eps.in."+str(case))
        dPerm = _read_RotPermList(gpath+"get_sLV_fixing_operations_dPerm.in."+str(case))
        fixOp, rotPerm, degeneracy = _get_sLV_fixing_operations(HNF,pLV,rot,shift,dPerm,eps)
        rotPerm_out = _read_RotPermList(gpath+"get_sLV_fixing_operations_rotPerm.out."+str(case))
        degen_out = _read_int(gpath+"get_sLV_fixing_operations_degeneracy.out."+str(case))
        fixOp_out = _read_fixOp(gpath+"get_sLV_fixing_operations_fixOp.out."+str(case))

        self._compare_outputs([fixOp,rotPerm,degeneracy],[fixOp_out,rotPerm_out,degen_out])
        
    def test_15(self):
        from phenum.grouptheory import _get_sLV_fixing_operations
        case = 140
        HNF = _read_int_2D(gpath+"get_sLV_fixing_operations_HNF.in."+str(case))
        pLV = np.transpose(_read_float_2D(gpath+"get_sLV_fixing_operations_pLV.in."+str(case)))
        nD = _read_int(gpath+"get_sLV_fixing_operations_nD.in."+str(case))
        rot = _read_float_3D(gpath+"get_sLV_fixing_operations_rot.in."+str(case))
        shift = list(map(list,zip(*_read_float_2D(gpath+"get_sLV_fixing_operations_shift.in."+str(case)))))
        eps = _read_float(gpath+"get_sLV_fixing_operations_eps.in."+str(case))
        dPerm = _read_RotPermList(gpath+"get_sLV_fixing_operations_dPerm.in."+str(case))
        fixOp, rotPerm, degeneracy = _get_sLV_fixing_operations(HNF,pLV,rot,shift,dPerm,eps)
        rotPerm_out = _read_RotPermList(gpath+"get_sLV_fixing_operations_rotPerm.out."+str(case))
        degen_out = _read_int(gpath+"get_sLV_fixing_operations_degeneracy.out."+str(case))
        fixOp_out = _read_fixOp(gpath+"get_sLV_fixing_operations_fixOp.out."+str(case))

        self._compare_outputs([fixOp,rotPerm,degeneracy],[fixOp_out,rotPerm_out,degen_out])
        
    def test_16(self):
        from phenum.grouptheory import _get_sLV_fixing_operations
        case = 150
        HNF = _read_int_2D(gpath+"get_sLV_fixing_operations_HNF.in."+str(case))
        pLV = np.transpose(_read_float_2D(gpath+"get_sLV_fixing_operations_pLV.in."+str(case)))
        nD = _read_int(gpath+"get_sLV_fixing_operations_nD.in."+str(case))
        rot = _read_float_3D(gpath+"get_sLV_fixing_operations_rot.in."+str(case))
        shift = list(map(list,zip(*_read_float_2D(gpath+"get_sLV_fixing_operations_shift.in."+str(case)))))
        eps = _read_float(gpath+"get_sLV_fixing_operations_eps.in."+str(case))
        dPerm = _read_RotPermList(gpath+"get_sLV_fixing_operations_dPerm.in."+str(case))
        fixOp, rotPerm, degeneracy = _get_sLV_fixing_operations(HNF,pLV,rot,shift,dPerm,eps)
        rotPerm_out = _read_RotPermList(gpath+"get_sLV_fixing_operations_rotPerm.out."+str(case))
        degen_out = _read_int(gpath+"get_sLV_fixing_operations_degeneracy.out."+str(case))
        fixOp_out = _read_fixOp(gpath+"get_sLV_fixing_operations_fixOp.out."+str(case))

        self._compare_outputs([fixOp,rotPerm,degeneracy],[fixOp_out,rotPerm_out,degen_out])
            
class TestMapDvectorPermutation(ut.TestCase):
    """Tests of the _map_dvector_permutation subroutine."""

    def _compare_outputs(self,out1,out2):
        if len(out1) == len(out2):
            for i in range(len(out1)):
                self.assertEqual(out1[i],out2[i])
        else:
            self.assertEqual(len(out1),len(out2))
                
    def test_1(self):
        from phenum.grouptheory import _map_dvector_permutation
        case = 1
        rd = list(map(list,zip(*_read_float_2D(gpath+"map_dvector_permutation_rd.in."+str(case)))))
        d = list(map(list,zip(*_read_float_2D(gpath+"map_dvector_permutation_d.in."+str(case)))))
        eps = _read_float(gpath+"map_dvector_permutation_eps.in."+str(case))
        n = len(rd)
        out = _read_int_1D(gpath+"map_dvector_permutation_RP.out."+str(case))
        out = [i-1 for i in out]

        self._compare_outputs(_map_dvector_permutation(rd,d,eps),out)
        
    def test_2(self):
        from phenum.grouptheory import _map_dvector_permutation
        case = 10
        rd = list(map(list,zip(*_read_float_2D(gpath+"map_dvector_permutation_rd.in."+str(case)))))
        d = list(map(list,zip(*_read_float_2D(gpath+"map_dvector_permutation_d.in."+str(case)))))
        eps = _read_float(gpath+"map_dvector_permutation_eps.in."+str(case))
        n = len(rd)
        out = _read_int_1D(gpath+"map_dvector_permutation_RP.out."+str(case))
        out = [i-1 for i in out]

        self._compare_outputs(_map_dvector_permutation(rd,d,eps),out)
        
    def test_3(self):
        from phenum.grouptheory import _map_dvector_permutation
        case = 20
        rd = list(map(list,zip(*_read_float_2D(gpath+"map_dvector_permutation_rd.in."+str(case)))))
        d = list(map(list,zip(*_read_float_2D(gpath+"map_dvector_permutation_d.in."+str(case)))))
        eps = _read_float(gpath+"map_dvector_permutation_eps.in."+str(case))
        n = len(rd)
        out = _read_int_1D(gpath+"map_dvector_permutation_RP.out."+str(case))
        out = [i-1 for i in out]

        self._compare_outputs(_map_dvector_permutation(rd,d,eps),out)
        
    def test_4(self):
        from phenum.grouptheory import _map_dvector_permutation
        case = 30
        rd = list(map(list,zip(*_read_float_2D(gpath+"map_dvector_permutation_rd.in."+str(case)))))
        d = list(map(list,zip(*_read_float_2D(gpath+"map_dvector_permutation_d.in."+str(case)))))
        eps = _read_float(gpath+"map_dvector_permutation_eps.in."+str(case))
        n = len(rd)
        out = _read_int_1D(gpath+"map_dvector_permutation_RP.out."+str(case))
        out = [i-1 for i in out]

        self._compare_outputs(_map_dvector_permutation(rd,d,eps),out)
        
    def test_5(self):
        from phenum.grouptheory import _map_dvector_permutation
        case = 40
        rd = list(map(list,zip(*_read_float_2D(gpath+"map_dvector_permutation_rd.in."+str(case)))))
        d = list(map(list,zip(*_read_float_2D(gpath+"map_dvector_permutation_d.in."+str(case)))))
        eps = _read_float(gpath+"map_dvector_permutation_eps.in."+str(case))
        n = len(rd)
        out = _read_int_1D(gpath+"map_dvector_permutation_RP.out."+str(case))
        out = [i-1 for i in out]

        self._compare_outputs(_map_dvector_permutation(rd,d,eps),out)
        
    def test_6(self):
        from phenum.grouptheory import _map_dvector_permutation
        case = 50
        rd = list(map(list,zip(*_read_float_2D(gpath+"map_dvector_permutation_rd.in."+str(case)))))
        d = list(map(list,zip(*_read_float_2D(gpath+"map_dvector_permutation_d.in."+str(case)))))
        eps = _read_float(gpath+"map_dvector_permutation_eps.in."+str(case))
        n = len(rd)
        out = _read_int_1D(gpath+"map_dvector_permutation_RP.out."+str(case))
        out = [i-1 for i in out]

        self._compare_outputs(_map_dvector_permutation(rd,d,eps),out)
        
    def test_7(self):
        from phenum.grouptheory import _map_dvector_permutation
        case = 60
        rd = list(map(list,zip(*_read_float_2D(gpath+"map_dvector_permutation_rd.in."+str(case)))))
        d = list(map(list,zip(*_read_float_2D(gpath+"map_dvector_permutation_d.in."+str(case)))))
        eps = _read_float(gpath+"map_dvector_permutation_eps.in."+str(case))
        n = len(rd)
        out = _read_int_1D(gpath+"map_dvector_permutation_RP.out."+str(case))
        out = [i-1 for i in out]

        self._compare_outputs(_map_dvector_permutation(rd,d,eps),out)
        
    def test_8(self):
        from phenum.grouptheory import _map_dvector_permutation
        case = 70
        rd = list(map(list,zip(*_read_float_2D(gpath+"map_dvector_permutation_rd.in."+str(case)))))
        d = list(map(list,zip(*_read_float_2D(gpath+"map_dvector_permutation_d.in."+str(case)))))
        eps = _read_float(gpath+"map_dvector_permutation_eps.in."+str(case))
        n = len(rd)
        out = _read_int_1D(gpath+"map_dvector_permutation_RP.out."+str(case))
        out = [i-1 for i in out]

        self._compare_outputs(_map_dvector_permutation(rd,d,eps),out)
        
    def test_9(self):
        from phenum.grouptheory import _map_dvector_permutation
        case = 80
        rd = list(map(list,zip(*_read_float_2D(gpath+"map_dvector_permutation_rd.in."+str(case)))))
        d = list(map(list,zip(*_read_float_2D(gpath+"map_dvector_permutation_d.in."+str(case)))))
        eps = _read_float(gpath+"map_dvector_permutation_eps.in."+str(case))
        n = len(rd)
        out = _read_int_1D(gpath+"map_dvector_permutation_RP.out."+str(case))
        out = [i-1 for i in out]

        self._compare_outputs(_map_dvector_permutation(rd,d,eps),out)
        
    def test_10(self):
        from phenum.grouptheory import _map_dvector_permutation
        case = 90
        rd = list(map(list,zip(*_read_float_2D(gpath+"map_dvector_permutation_rd.in."+str(case)))))
        d = list(map(list,zip(*_read_float_2D(gpath+"map_dvector_permutation_d.in."+str(case)))))
        eps = _read_float(gpath+"map_dvector_permutation_eps.in."+str(case))
        n = len(rd)
        out = _read_int_1D(gpath+"map_dvector_permutation_RP.out."+str(case))
        out = [i-1 for i in out]

        self._compare_outputs(_map_dvector_permutation(rd,d,eps),out)
        
    def test_11(self):
        from phenum.grouptheory import _map_dvector_permutation
        case = 100
        rd = list(map(list,zip(*_read_float_2D(gpath+"map_dvector_permutation_rd.in."+str(case)))))
        d = list(map(list,zip(*_read_float_2D(gpath+"map_dvector_permutation_d.in."+str(case)))))
        eps = _read_float(gpath+"map_dvector_permutation_eps.in."+str(case))
        n = len(rd)
        out = _read_int_1D(gpath+"map_dvector_permutation_RP.out."+str(case))
        out = [i-1 for i in out]

        self._compare_outputs(_map_dvector_permutation(rd,d,eps),out)
        
    def test_12(self):
        from phenum.grouptheory import _map_dvector_permutation
        case = 110
        rd = list(map(list,zip(*_read_float_2D(gpath+"map_dvector_permutation_rd.in."+str(case)))))
        d = list(map(list,zip(*_read_float_2D(gpath+"map_dvector_permutation_d.in."+str(case)))))
        eps = _read_float(gpath+"map_dvector_permutation_eps.in."+str(case))
        n = len(rd)
        out = _read_int_1D(gpath+"map_dvector_permutation_RP.out."+str(case))
        out = [i-1 for i in out]

        self._compare_outputs(_map_dvector_permutation(rd,d,eps),out)
        
    def test_13(self):
        from phenum.grouptheory import _map_dvector_permutation
        case = 120
        rd = list(map(list,zip(*_read_float_2D(gpath+"map_dvector_permutation_rd.in."+str(case)))))
        d = list(map(list,zip(*_read_float_2D(gpath+"map_dvector_permutation_d.in."+str(case)))))
        eps = _read_float(gpath+"map_dvector_permutation_eps.in."+str(case))
        n = len(rd)
        out = _read_int_1D(gpath+"map_dvector_permutation_RP.out."+str(case))
        out = [i-1 for i in out]

        self._compare_outputs(_map_dvector_permutation(rd,d,eps),out)
        
    def test_14(self):
        from phenum.grouptheory import _map_dvector_permutation
        case = 130
        rd = list(map(list,zip(*_read_float_2D(gpath+"map_dvector_permutation_rd.in."+str(case)))))
        d = list(map(list,zip(*_read_float_2D(gpath+"map_dvector_permutation_d.in."+str(case)))))
        eps = _read_float(gpath+"map_dvector_permutation_eps.in."+str(case))
        n = len(rd)
        out = _read_int_1D(gpath+"map_dvector_permutation_RP.out."+str(case))
        out = [i-1 for i in out]

        self._compare_outputs(_map_dvector_permutation(rd,d,eps),out)
        
    def test_15(self):
        from phenum.grouptheory import _map_dvector_permutation
        case = 140
        rd = list(map(list,zip(*_read_float_2D(gpath+"map_dvector_permutation_rd.in."+str(case)))))
        d = list(map(list,zip(*_read_float_2D(gpath+"map_dvector_permutation_d.in."+str(case)))))
        eps = _read_float(gpath+"map_dvector_permutation_eps.in."+str(case))
        n = len(rd)
        out = _read_int_1D(gpath+"map_dvector_permutation_RP.out."+str(case))
        out = [i-1 for i in out]

        self._compare_outputs(_map_dvector_permutation(rd,d,eps),out)
        
    def test_16(self):
        from phenum.grouptheory import _map_dvector_permutation
        case = 150
        rd = list(map(list,zip(*_read_float_2D(gpath+"map_dvector_permutation_rd.in."+str(case)))))
        d = list(map(list,zip(*_read_float_2D(gpath+"map_dvector_permutation_d.in."+str(case)))))
        eps = _read_float(gpath+"map_dvector_permutation_eps.in."+str(case))
        n = len(rd)
        out = _read_int_1D(gpath+"map_dvector_permutation_RP.out."+str(case))
        out = [i-1 for i in out]

        self._compare_outputs(_map_dvector_permutation(rd,d,eps),out)
        
    def test_17(self):
        from phenum.grouptheory import _map_dvector_permutation
        case = 160
        rd = list(map(list,zip(*_read_float_2D(gpath+"map_dvector_permutation_rd.in."+str(case)))))
        d = list(map(list,zip(*_read_float_2D(gpath+"map_dvector_permutation_d.in."+str(case)))))
        eps = _read_float(gpath+"map_dvector_permutation_eps.in."+str(case))
        n = len(rd)
        out = _read_int_1D(gpath+"map_dvector_permutation_RP.out."+str(case))
        out = [i-1 for i in out]

        self._compare_outputs(_map_dvector_permutation(rd,d,eps),out)
        
    def test_18(self):
        from phenum.grouptheory import _map_dvector_permutation
        case = 170
        rd = list(map(list,zip(*_read_float_2D(gpath+"map_dvector_permutation_rd.in."+str(case)))))
        d = list(map(list,zip(*_read_float_2D(gpath+"map_dvector_permutation_d.in."+str(case)))))
        eps = _read_float(gpath+"map_dvector_permutation_eps.in."+str(case))
        n = len(rd)
        out = _read_int_1D(gpath+"map_dvector_permutation_RP.out."+str(case))
        out = [i-1 for i in out]

        self._compare_outputs(_map_dvector_permutation(rd,d,eps),out)
        
    def test_19(self):
        from phenum.grouptheory import _map_dvector_permutation
        case = 180
        rd = list(map(list,zip(*_read_float_2D(gpath+"map_dvector_permutation_rd.in."+str(case)))))
        d = list(map(list,zip(*_read_float_2D(gpath+"map_dvector_permutation_d.in."+str(case)))))
        eps = _read_float(gpath+"map_dvector_permutation_eps.in."+str(case))
        n = len(rd)
        out = _read_int_1D(gpath+"map_dvector_permutation_RP.out."+str(case))
        out = [i-1 for i in out]

        self._compare_outputs(_map_dvector_permutation(rd,d,eps),out)
        
    def test_20(self):
        from phenum.grouptheory import _map_dvector_permutation
        case = 190
        rd = list(map(list,zip(*_read_float_2D(gpath+"map_dvector_permutation_rd.in."+str(case)))))
        d = list(map(list,zip(*_read_float_2D(gpath+"map_dvector_permutation_d.in."+str(case)))))
        eps = _read_float(gpath+"map_dvector_permutation_eps.in."+str(case))
        n = len(rd)
        out = _read_int_1D(gpath+"map_dvector_permutation_RP.out."+str(case))
        out = [i-1 for i in out]

        self._compare_outputs(_map_dvector_permutation(rd,d,eps),out)
        
    def test_21(self):
        from phenum.grouptheory import _map_dvector_permutation
        case = 200
        rd = list(map(list,zip(*_read_float_2D(gpath+"map_dvector_permutation_rd.in."+str(case)))))
        d = list(map(list,zip(*_read_float_2D(gpath+"map_dvector_permutation_d.in."+str(case)))))
        eps = _read_float(gpath+"map_dvector_permutation_eps.in."+str(case))
        n = len(rd)
        out = _read_int_1D(gpath+"map_dvector_permutation_RP.out."+str(case))
        out = [i-1 for i in out]

        self._compare_outputs(_map_dvector_permutation(rd,d,eps),out)
                
class TestFindMinmaxIndices(ut.TestCase):
    """Tests of the _find_minmax_indices subroutine."""

    def test_1(self):
        from phenum.grouptheory import _find_minmax_indices
        case = 1
        invec = _read_int_1D(gpath+"get_minmax_indices_invec.in."+str(case))
        min,max = _find_minmax_indices(invec)
        min_out = _read_int(gpath+"get_minmax_indices_min.out."+str(case))-1
        max_out = _read_int(gpath+"get_minmax_indices_max.out."+str(case))-1
        self.assertEqual(min,min_out)
        self.assertEqual(max,max_out)

    def test_2(self):
        from phenum.grouptheory import _find_minmax_indices
        case = 5
        invec = _read_int_1D(gpath+"get_minmax_indices_invec.in."+str(case))
        min,max = _find_minmax_indices(invec)
        min_out = _read_int(gpath+"get_minmax_indices_min.out."+str(case))-1
        max_out = _read_int(gpath+"get_minmax_indices_max.out."+str(case))-1
        self.assertEqual(min,min_out)
        self.assertEqual(max,max_out)

    def test_3(self):
        from phenum.grouptheory import _find_minmax_indices
        case = 10
        invec = _read_int_1D(gpath+"get_minmax_indices_invec.in."+str(case))
        min,max = _find_minmax_indices(invec)
        min_out = _read_int(gpath+"get_minmax_indices_min.out."+str(case))-1
        max_out = _read_int(gpath+"get_minmax_indices_max.out."+str(case))-1
        self.assertEqual(min,min_out)
        self.assertEqual(max,max_out)

    def test_4(self):
        from phenum.grouptheory import _find_minmax_indices
        case = 15
        invec = _read_int_1D(gpath+"get_minmax_indices_invec.in."+str(case))
        min,max = _find_minmax_indices(invec)
        min_out = _read_int(gpath+"get_minmax_indices_min.out."+str(case))-1
        max_out = _read_int(gpath+"get_minmax_indices_max.out."+str(case))-1
        self.assertEqual(min,min_out)
        self.assertEqual(max,max_out)

    def test_5(self):
        from phenum.grouptheory import _find_minmax_indices
        case = 20
        invec = _read_int_1D(gpath+"get_minmax_indices_invec.in."+str(case))
        min,max = _find_minmax_indices(invec)
        min_out = _read_int(gpath+"get_minmax_indices_min.out."+str(case))-1
        max_out = _read_int(gpath+"get_minmax_indices_max.out."+str(case))-1
        self.assertEqual(min,min_out)
        self.assertEqual(max,max_out)

    def test_6(self):
        from phenum.grouptheory import _find_minmax_indices
        case = 25
        invec = _read_int_1D(gpath+"get_minmax_indices_invec.in."+str(case))
        min,max = _find_minmax_indices(invec)
        min_out = _read_int(gpath+"get_minmax_indices_min.out."+str(case))-1
        max_out = _read_int(gpath+"get_minmax_indices_max.out."+str(case))-1
        self.assertEqual(min,min_out)
        self.assertEqual(max,max_out)

    def test_7(self):
        from phenum.grouptheory import _find_minmax_indices
        case = 30
        invec = _read_int_1D(gpath+"get_minmax_indices_invec.in."+str(case))
        min,max = _find_minmax_indices(invec)
        min_out = _read_int(gpath+"get_minmax_indices_min.out."+str(case))-1
        max_out = _read_int(gpath+"get_minmax_indices_max.out."+str(case))-1
        self.assertEqual(min,min_out)
        self.assertEqual(max,max_out)

    def test_8(self):
        from phenum.grouptheory import _find_minmax_indices
        case = 35
        invec = _read_int_1D(gpath+"get_minmax_indices_invec.in."+str(case))
        min,max = _find_minmax_indices(invec)
        min_out = _read_int(gpath+"get_minmax_indices_min.out."+str(case))-1
        max_out = _read_int(gpath+"get_minmax_indices_max.out."+str(case))-1
        self.assertEqual(min,min_out)
        self.assertEqual(max,max_out)

    def test_9(self):
        from phenum.grouptheory import _find_minmax_indices
        case = 40
        invec = _read_int_1D(gpath+"get_minmax_indices_invec.in."+str(case))
        min,max = _find_minmax_indices(invec)
        min_out = _read_int(gpath+"get_minmax_indices_min.out."+str(case))-1
        max_out = _read_int(gpath+"get_minmax_indices_max.out."+str(case))-1
        self.assertEqual(min,min_out)
        self.assertEqual(max,max_out)

    def test_10(self):
        from phenum.grouptheory import _find_minmax_indices
        case = 45
        invec = _read_int_1D(gpath+"get_minmax_indices_invec.in."+str(case))
        min,max = _find_minmax_indices(invec)
        min_out = _read_int(gpath+"get_minmax_indices_min.out."+str(case))-1
        max_out = _read_int(gpath+"get_minmax_indices_max.out."+str(case))-1
        self.assertEqual(min,min_out)
        self.assertEqual(max,max_out)

    def test_11(self):
        from phenum.grouptheory import _find_minmax_indices
        case = 50
        invec = _read_int_1D(gpath+"get_minmax_indices_invec.in."+str(case))
        min,max = _find_minmax_indices(invec)
        min_out = _read_int(gpath+"get_minmax_indices_min.out."+str(case))-1
        max_out = _read_int(gpath+"get_minmax_indices_max.out."+str(case))-1
        self.assertEqual(min,min_out)
        self.assertEqual(max,max_out)

class TestGetDvectorPermutations(ut.TestCase):
    """Tests of the _get_dvector_permutations subroutine."""

    def _compare_outputs(self,rot1,rot2):
        # self.assertEqual(rot1.nL,rot2.nL)
        self.assertEqual(rot1.RotIndx, rot2.RotIndx)

        if (rot1.v != None) and (rot2.v != None):
            if len(rot1.v) == len(rot2.v):
                for i in range(len(rot1.v)):
                    for j in range(len(rot1.v[i])):
                        for k in range(len(rot1.v[i][j])):
                            self.assertAlmostEqual(rot1.v[i][j][k],rot2.v[i][j][k],places=12)
            else:
                self.assertEqual(len(rot1.v),len(rot2.v))
        else:
            self.assertEqual(rot1.v,rot2.v)

        if (rot1.perm.site_perm != None) and (rot2.perm.site_perm != None):
            if len(rot1.perm.site_perm) == len(rot2.perm.site_perm):
                for i in range(len(rot1.perm.site_perm)):
                    for j in range(len(rot1.perm.site_perm[i])):
                        self.assertEqual(rot1.perm.site_perm[i][j],rot2.perm.site_perm[i][j])
            else:
                self.assertEqual(len(rot1.perm.site_perm),len(rot2.perm.site_perm))
        else:
            self.assertEqual(rot1.perm.site_perm,rot2.perm.site_perm)
            
        if (rot1.perm.arrow_perm != None) and (rot2.perm.arrow_perm != None):
            if len(rot1.perm.arrow_perm) == len(rot2.perm.arrow_perm):
                for i in range(len(rot1.perm.arrow_perm)):
                    for j in range(len(rot1.perm.arrow_perm[i])):
                        self.assertEqual(rot1.perm.arrow_perm[i][j],rot2.perm.arrow_perm[i][j])
            else:
                self.assertEqual(len(rot1.perm.arrow_perm),len(rot2.perm.arrow_perm))
        else:
            self.assertEqual(rot1.perm.arrow_perm,rot2.perm.arrow_perm)
        
    def test_1(self):
        from phenum.grouptheory import _get_dvector_permutations
        case = 1
        par_lat = list(map(list,zip(*_read_float_2D(gpath+"get_dvector_permutations_pLV.in."+str(case)))))
        bas_vecs = list(map(list,zip(*_read_float_2D(gpath+"get_dvector_permutations_pd.in."+str(case)))))
        LatDim = _read_int(gpath+"get_dvector_permutations_LatDim.in."+str(case))
        eps = _read_float(gpath+"get_dvector_permutations_eps.in."+str(case))
        dRPList_out = _read_RotPermList(gpath+"get_dvector_permutations_dRPList.out."+str(case))
        dRPList = _get_dvector_permutations(par_lat,bas_vecs,LatDim,eps)
        self._compare_outputs(dRPList,dRPList_out)
        
    def test_2(self):
        from phenum.grouptheory import _get_dvector_permutations
        case = 2
        par_lat = list(map(list,zip(*_read_float_2D(gpath+"get_dvector_permutations_pLV.in."+str(case)))))
        bas_vecs = list(map(list,zip(*_read_float_2D(gpath+"get_dvector_permutations_pd.in."+str(case)))))
        LatDim = _read_int(gpath+"get_dvector_permutations_LatDim.in."+str(case))
        eps = _read_float(gpath+"get_dvector_permutations_eps.in."+str(case))
        dRPList_out = _read_RotPermList(gpath+"get_dvector_permutations_dRPList.out."+str(case))
        dRPList = _get_dvector_permutations(par_lat,bas_vecs,LatDim,eps)
        self._compare_outputs(dRPList,dRPList_out)

    def test_3(self):
        from phenum.grouptheory import _get_dvector_permutations
        case = 3
        par_lat = list(map(list,zip(*_read_float_2D(gpath+"get_dvector_permutations_pLV.in."+str(case)))))
        bas_vecs = list(map(list,zip(*_read_float_2D(gpath+"get_dvector_permutations_pd.in."+str(case)))))
        LatDim = _read_int(gpath+"get_dvector_permutations_LatDim.in."+str(case))
        eps = _read_float(gpath+"get_dvector_permutations_eps.in."+str(case))
        dRPList_out = _read_RotPermList(gpath+"get_dvector_permutations_dRPList.out."+str(case))
        dRPList = _get_dvector_permutations(par_lat,bas_vecs,LatDim,eps)
        self._compare_outputs(dRPList,dRPList_out)

    def test_4(self):
        from phenum.grouptheory import _get_dvector_permutations
        case = 4
        par_lat = list(map(list,zip(*_read_float_2D(gpath+"get_dvector_permutations_pLV.in."+str(case)))))
        bas_vecs = list(map(list,zip(*_read_float_2D(gpath+"get_dvector_permutations_pd.in."+str(case)))))
        LatDim = _read_int(gpath+"get_dvector_permutations_LatDim.in."+str(case))
        eps = _read_float(gpath+"get_dvector_permutations_eps.in."+str(case))
        dRPList_out = _read_RotPermList(gpath+"get_dvector_permutations_dRPList.out."+str(case))
        dRPList = _get_dvector_permutations(par_lat,bas_vecs,LatDim,eps)
        self._compare_outputs(dRPList,dRPList_out)

    def test_5(self):
        from phenum.grouptheory import _get_dvector_permutations
        case = 5
        par_lat = list(map(list,zip(*_read_float_2D(gpath+"get_dvector_permutations_pLV.in."+str(case)))))
        bas_vecs = list(map(list,zip(*_read_float_2D(gpath+"get_dvector_permutations_pd.in."+str(case)))))
        LatDim = _read_int(gpath+"get_dvector_permutations_LatDim.in."+str(case))
        eps = _read_float(gpath+"get_dvector_permutations_eps.in."+str(case))
        dRPList_out = _read_RotPermList(gpath+"get_dvector_permutations_dRPList.out."+str(case))
        dRPList = _get_dvector_permutations(par_lat,bas_vecs,LatDim,eps)
        self._compare_outputs(dRPList,dRPList_out)

    def test_6(self):
        from phenum.grouptheory import _get_dvector_permutations
        case = 6
        par_lat = list(map(list,zip(*_read_float_2D(gpath+"get_dvector_permutations_pLV.in."+str(case)))))
        bas_vecs = list(map(list,zip(*_read_float_2D(gpath+"get_dvector_permutations_pd.in."+str(case)))))
        LatDim = _read_int(gpath+"get_dvector_permutations_LatDim.in."+str(case))
        eps = _read_float(gpath+"get_dvector_permutations_eps.in."+str(case))
        dRPList_out = _read_RotPermList(gpath+"get_dvector_permutations_dRPList.out."+str(case))
        dRPList = _get_dvector_permutations(par_lat,bas_vecs,LatDim,eps)
        self._compare_outputs(dRPList,dRPList_out)

    def test_7(self):
        from phenum.grouptheory import _get_dvector_permutations
        case = 7
        par_lat = list(map(list,zip(*_read_float_2D(gpath+"get_dvector_permutations_pLV.in."+str(case)))))
        bas_vecs = list(map(list,zip(*_read_float_2D(gpath+"get_dvector_permutations_pd.in."+str(case)))))
        LatDim = _read_int(gpath+"get_dvector_permutations_LatDim.in."+str(case))
        eps = _read_float(gpath+"get_dvector_permutations_eps.in."+str(case))
        dRPList_out = _read_RotPermList(gpath+"get_dvector_permutations_dRPList.out."+str(case))
        dRPList = _get_dvector_permutations(par_lat,bas_vecs,LatDim,eps)
        self._compare_outputs(dRPList,dRPList_out)

    def test_8(self):
        from phenum.grouptheory import _get_dvector_permutations
        case = 8
        par_lat = list(map(list,zip(*_read_float_2D(gpath+"get_dvector_permutations_pLV.in."+str(case)))))
        bas_vecs = list(map(list,zip(*_read_float_2D(gpath+"get_dvector_permutations_pd.in."+str(case)))))
        LatDim = _read_int(gpath+"get_dvector_permutations_LatDim.in."+str(case))
        eps = _read_float(gpath+"get_dvector_permutations_eps.in."+str(case))
        dRPList_out = _read_RotPermList(gpath+"get_dvector_permutations_dRPList.out."+str(case))
        dRPList = _get_dvector_permutations(par_lat,bas_vecs,LatDim,eps)
        self._compare_outputs(dRPList,dRPList_out)
    
    def test_9(self):
        from phenum.grouptheory import _get_dvector_permutations
        case = 9
        par_lat = list(map(list,zip(*_read_float_2D(gpath+"get_dvector_permutations_pLV.in."+str(case)))))
        bas_vecs = list(map(list,zip(*_read_float_2D(gpath+"get_dvector_permutations_pd.in."+str(case)))))
        LatDim = _read_int(gpath+"get_dvector_permutations_LatDim.in."+str(case))
        eps = _read_float(gpath+"get_dvector_permutations_eps.in."+str(case))
        dRPList = _get_dvector_permutations(par_lat,bas_vecs,LatDim,eps)
        dRPList_out = _read_RotPermList(gpath+"get_dvector_permutations_dRPList.out."+str(case))
        self._compare_outputs(dRPList,dRPList_out)

    def test_10(self):
        from phenum.grouptheory import _get_dvector_permutations
        case = 10
        par_lat = list(map(list,zip(*_read_float_2D(gpath+"get_dvector_permutations_pLV.in."+str(case)))))
        bas_vecs = list(map(list,zip(*_read_float_2D(gpath+"get_dvector_permutations_pd.in."+str(case)))))
        LatDim = _read_int(gpath+"get_dvector_permutations_LatDim.in."+str(case))
        eps = _read_float(gpath+"get_dvector_permutations_eps.in."+str(case))
        dRPList_out = _read_RotPermList(gpath+"get_dvector_permutations_dRPList.out."+str(case))
        dRPList = _get_dvector_permutations(par_lat,bas_vecs,LatDim,eps)
        self._compare_outputs(dRPList,dRPList_out)

    def test_11(self):
        from phenum.grouptheory import _get_dvector_permutations
        case = 11
        par_lat = list(map(list,zip(*_read_float_2D(gpath+"get_dvector_permutations_pLV.in."+str(case)))))
        bas_vecs = list(map(list,zip(*_read_float_2D(gpath+"get_dvector_permutations_pd.in."+str(case)))))
        LatDim = _read_int(gpath+"get_dvector_permutations_LatDim.in."+str(case))
        eps = _read_float(gpath+"get_dvector_permutations_eps.in."+str(case))
        dRPList_out = _read_RotPermList(gpath+"get_dvector_permutations_dRPList.out."+str(case))
        dRPList = _get_dvector_permutations(par_lat,bas_vecs,LatDim,eps)
        self._compare_outputs(dRPList,dRPList_out)

    def test_12(self):
        from phenum.grouptheory import _get_dvector_permutations
        case = 12
        par_lat = list(map(list,zip(*_read_float_2D(gpath+"get_dvector_permutations_pLV.in."+str(case)))))
        bas_vecs = list(map(list,zip(*_read_float_2D(gpath+"get_dvector_permutations_pd.in."+str(case)))))
        LatDim = _read_int(gpath+"get_dvector_permutations_LatDim.in."+str(case))
        eps = _read_float(gpath+"get_dvector_permutations_eps.in."+str(case))
        dRPList_out = _read_RotPermList(gpath+"get_dvector_permutations_dRPList.out."+str(case))
        dRPList = _get_dvector_permutations(par_lat,bas_vecs,LatDim,eps)
        self._compare_outputs(dRPList,dRPList_out)

    def test_13(self):
        from phenum.grouptheory import _get_dvector_permutations
        case = 13
        par_lat = list(map(list,zip(*_read_float_2D(gpath+"get_dvector_permutations_pLV.in."+str(case)))))
        bas_vecs = list(map(list,zip(*_read_float_2D(gpath+"get_dvector_permutations_pd.in."+str(case)))))
        LatDim = _read_int(gpath+"get_dvector_permutations_LatDim.in."+str(case))
        eps = _read_float(gpath+"get_dvector_permutations_eps.in."+str(case))
        dRPList_out = _read_RotPermList(gpath+"get_dvector_permutations_dRPList.out."+str(case))
        dRPList = _get_dvector_permutations(par_lat,bas_vecs,LatDim,eps)
        self._compare_outputs(dRPList,dRPList_out)

    def test_14(self):
        from phenum.grouptheory import _get_dvector_permutations
        case = 14
        par_lat = list(map(list,zip(*_read_float_2D(gpath+"get_dvector_permutations_pLV.in."+str(case)))))
        bas_vecs = list(map(list,zip(*_read_float_2D(gpath+"get_dvector_permutations_pd.in."+str(case)))))
        LatDim = _read_int(gpath+"get_dvector_permutations_LatDim.in."+str(case))
        eps = _read_float(gpath+"get_dvector_permutations_eps.in."+str(case))
        dRPList_out = _read_RotPermList(gpath+"get_dvector_permutations_dRPList.out."+str(case))
        dRPList = _get_dvector_permutations(par_lat,bas_vecs,LatDim,eps)
        self._compare_outputs(dRPList,dRPList_out)

    def test_15(self):
        from phenum.grouptheory import _get_dvector_permutations
        case = 15
        par_lat = list(map(list,zip(*_read_float_2D(gpath+"get_dvector_permutations_pLV.in."+str(case)))))
        bas_vecs = list(map(list,zip(*_read_float_2D(gpath+"get_dvector_permutations_pd.in."+str(case)))))
        LatDim = _read_int(gpath+"get_dvector_permutations_LatDim.in."+str(case))
        eps = _read_float(gpath+"get_dvector_permutations_eps.in."+str(case))
        dRPList_out = _read_RotPermList(gpath+"get_dvector_permutations_dRPList.out."+str(case))
        dRPList = _get_dvector_permutations(par_lat,bas_vecs,LatDim,eps)
        self._compare_outputs(dRPList,dRPList_out)

    def test_16(self):
        from phenum.grouptheory import _get_dvector_permutations
        case = 16
        par_lat = list(map(list,zip(*_read_float_2D(gpath+"get_dvector_permutations_pLV.in."+str(case)))))
        bas_vecs = list(map(list,zip(*_read_float_2D(gpath+"get_dvector_permutations_pd.in."+str(case)))))
        LatDim = _read_int(gpath+"get_dvector_permutations_LatDim.in."+str(case))
        eps = _read_float(gpath+"get_dvector_permutations_eps.in."+str(case))
        dRPList_out = _read_RotPermList(gpath+"get_dvector_permutations_dRPList.out."+str(case))
        dRPList = _get_dvector_permutations(par_lat,bas_vecs,LatDim,eps)
        self._compare_outputs(dRPList,dRPList_out)

    def test_17(self):
        from phenum.grouptheory import _get_dvector_permutations
        case = 17
        par_lat = list(map(list,zip(*_read_float_2D(gpath+"get_dvector_permutations_pLV.in."+str(case)))))
        bas_vecs = list(map(list,zip(*_read_float_2D(gpath+"get_dvector_permutations_pd.in."+str(case)))))
        LatDim = _read_int(gpath+"get_dvector_permutations_LatDim.in."+str(case))
        eps = _read_float(gpath+"get_dvector_permutations_eps.in."+str(case))
        dRPList_out = _read_RotPermList(gpath+"get_dvector_permutations_dRPList.out."+str(case))
        dRPList = _get_dvector_permutations(par_lat,bas_vecs,LatDim,eps)
        self._compare_outputs(dRPList,dRPList_out)

    def test_18(self):
        from phenum.grouptheory import _get_dvector_permutations
        case = 18
        par_lat = list(map(list,zip(*_read_float_2D(gpath+"get_dvector_permutations_pLV.in."+str(case)))))
        bas_vecs = list(map(list,zip(*_read_float_2D(gpath+"get_dvector_permutations_pd.in."+str(case)))))
        LatDim = _read_int(gpath+"get_dvector_permutations_LatDim.in."+str(case))
        eps = _read_float(gpath+"get_dvector_permutations_eps.in."+str(case))
        dRPList_out = _read_RotPermList(gpath+"get_dvector_permutations_dRPList.out."+str(case))
        dRPList = _get_dvector_permutations(par_lat,bas_vecs,LatDim,eps)
        self._compare_outputs(dRPList,dRPList_out)

    def test_19(self):
        from phenum.grouptheory import _get_dvector_permutations
        case = 19
        par_lat = list(map(list,zip(*_read_float_2D(gpath+"get_dvector_permutations_pLV.in."+str(case)))))
        bas_vecs = list(map(list,zip(*_read_float_2D(gpath+"get_dvector_permutations_pd.in."+str(case)))))
        LatDim = _read_int(gpath+"get_dvector_permutations_LatDim.in."+str(case))
        eps = _read_float(gpath+"get_dvector_permutations_eps.in."+str(case))
        dRPList_out = _read_RotPermList(gpath+"get_dvector_permutations_dRPList.out."+str(case))
        dRPList = _get_dvector_permutations(par_lat,bas_vecs,LatDim,eps)
        self._compare_outputs(dRPList,dRPList_out)

    def test_20(self):
        from phenum.grouptheory import _get_dvector_permutations
        case = 20
        par_lat = list(map(list,zip(*_read_float_2D(gpath+"get_dvector_permutations_pLV.in."+str(case)))))
        bas_vecs = list(map(list,zip(*_read_float_2D(gpath+"get_dvector_permutations_pd.in."+str(case)))))
        LatDim = _read_int(gpath+"get_dvector_permutations_LatDim.in."+str(case))
        eps = _read_float(gpath+"get_dvector_permutations_eps.in."+str(case))
        dRPList_out = _read_RotPermList(gpath+"get_dvector_permutations_dRPList.out."+str(case))
        dRPList = _get_dvector_permutations(par_lat,bas_vecs,LatDim,eps)
        self._compare_outputs(dRPList,dRPList_out)

class TestGetRotationPermsLists(ut.TestCase):
    """Tests of the _get_rotation_perms_lists subroutine."""

    def _compare_outputs(self,out1,out2):
        if len(out1) == len(out2):
            for t in range(len(out1)):
                rot1 = out1[t]
                rot2 = out2[t]
                if rot1.nL == 0 or rot1.nL == None:
                    if rot2.nL == 0 or rot2.nL == None:
                        self.assertEqual(True,True)
                    else:
                        self.assertEqual(rot1.nL,rot2.nL)
                else:
                    self.assertEqual(rot1.nL,rot2.nL)

                self.assertEqual(rot1.RotIndx, rot2.RotIndx)

                if (rot1.v != None) and (rot2.v != None):
                    if len(rot1.v) == len(rot2.v):
                        for i in range(len(rot1.v)):
                            for j in range(len(rot1.v[i])):
                                for k in range(len(rot1.v[i][j])):
                                    self.assertAlmostEqual(rot1.v[i][j][k],rot2.v[i][j][k],places=12)
                    else:
                        self.assertEqual(len(rot1.v),len(rot2.v))
                else:
                    self.assertEqual(rot1.v,rot2.v)

                if (rot1.perm.site_perm != None) and (rot2.perm.site_perm != None):
                    if len(rot1.perm.site_perm) == len(rot2.perm.site_perm):
                        rot1.perm.site_perm = sorted(rot1.perm.site_perm)
                        rot2.perm.site_perm = sorted(rot2.perm.site_perm)
                        for i in range(len(rot1.perm.site_perm)):
                            for j in range(len(rot1.perm.site_perm[i])):
                                self.assertEqual(rot1.perm.site_perm[i][j],rot2.perm.site_perm[i][j])
                    else:
                        self.assertEqual(len(rot1.perm.site_perm),len(rot2.perm.site_perm))
                else:
                    self.assertEqual(rot1.perm.site_perm,rot2.perm.site_perm)
            
                if (rot1.perm.arrow_perm != None) and (rot2.perm.arrow_perm != None):
                    if len(rot1.perm.arrow_perm) == len(rot2.perm.arrow_perm):
                        for i in range(len(rot1.perm.arrow_perm)):
                            for j in range(len(rot1.perm.arrow_perm[i])):
                                self.assertEqual(rot1.perm.arrow_perm[i][j],rot2.perm.arrow_perm[i][j])
                    else:
                        self.assertEqual(len(rot1.perm.arrow_perm),len(rot2.perm.arrow_perm))
                else:
                    self.assertEqual(rot1.perm.arrow_perm,rot2.perm.arrow_perm)
        else:
            self.assertEqual(len(out1),len(out2))

    def test_1(self):
        from phenum.grouptheory import _get_rotation_perms_lists
        case = 1
        A = [[10,0,0],[0,10,0],[0,0,10]]
        HNF = [[[1,0,0],[0,1,0],[0,0,1]]]
        L = [[[1,0,0],[0,1,0],[0,0,1]]]
        SNF = [[[1,0,0],[0,1,0],[0,0,1]]]
        Op = _read_fixOp_1D(gpath+"get_rotation_perms_lists_Op.in."+str(case))
        dperms = _read_RotPermList(gpath+"get_rotation_perms_lists_dperms.in."+str(case))
        RPlist = _read_RotPermList_1D(gpath+"get_rotation_perms_lists_RPlist.in."+str(case))
        eps = _read_float(gpath+"get_rotation_perms_lists_eps.in."+str(case))
        out1 = _get_rotation_perms_lists(A,HNF,L,SNF,Op,RPlist,dperms,eps)
        for out in out1:
            out.perm.arrow_perm = None
        out2 = _read_RotPermList_1D(gpath+"get_rotation_perms_lists_RPlist.out."+str(case))
        self._compare_outputs(out1,out2)

    def test_2(self):
        from phenum.grouptheory import _get_rotation_perms_lists
        case = 2
        A = list(map(list,zip(*_read_float_2D(gpath+"get_rotation_perms_lists_A.in."+str(case)))))
        HNF = _read_int_3D(gpath+"get_rotation_perms_lists_HNF.in."+str(case))
        L = _read_int_3D(gpath+"get_rotation_perms_lists_L.in."+str(case))
        SNF = _read_int_3D(gpath+"get_rotation_perms_lists_SNF.in."+str(case))
        Op = _read_fixOp_1D(gpath+"get_rotation_perms_lists_Op.in."+str(case))
        dperms = _read_RotPermList(gpath+"get_rotation_perms_lists_dperms.in."+str(case))
        RPlist = _read_RotPermList_1D(gpath+"get_rotation_perms_lists_RPlist.in."+str(case))
        eps = _read_float(gpath+"get_rotation_perms_lists_eps.in."+str(case))
        out1 = _get_rotation_perms_lists(A,HNF,L,SNF,Op,RPlist,dperms,eps)
        for out in out1:
            out.perm.arrow_perm = None
        out2 = _read_RotPermList_1D(gpath+"get_rotation_perms_lists_RPlist.out."+str(case))
        self._compare_outputs(out1,out2)

    def test_3(self):
        from phenum.grouptheory import _get_rotation_perms_lists
        case = 3
        A = list(map(list,zip(*_read_float_2D(gpath+"get_rotation_perms_lists_A.in."+str(case)))))
        HNF = _read_int_3D(gpath+"get_rotation_perms_lists_HNF.in."+str(case))
        L = _read_int_3D(gpath+"get_rotation_perms_lists_L.in."+str(case))
        SNF = _read_int_3D(gpath+"get_rotation_perms_lists_SNF.in."+str(case))
        Op = _read_fixOp_1D(gpath+"get_rotation_perms_lists_Op.in."+str(case))
        dperms = _read_RotPermList(gpath+"get_rotation_perms_lists_dperms.in."+str(case))
        RPlist = _read_RotPermList_1D(gpath+"get_rotation_perms_lists_RPlist.in."+str(case))
        eps = _read_float(gpath+"get_rotation_perms_lists_eps.in."+str(case))
        out1 = _get_rotation_perms_lists(A,HNF,L,SNF,Op,RPlist,dperms,eps)
        out2 = _read_RotPermList_1D(gpath+"get_rotation_perms_lists_RPlist.out."+str(case))
        for out in out1:
            out.perm.arrow_perm = None
        self._compare_outputs(out1,out2)

    def test_4(self):
        from phenum.grouptheory import _get_rotation_perms_lists
        case = 4
        A = list(map(list,zip(*_read_float_2D(gpath+"get_rotation_perms_lists_A.in."+str(case)))))
        HNF = _read_int_3D(gpath+"get_rotation_perms_lists_HNF.in."+str(case))
        L = _read_int_3D(gpath+"get_rotation_perms_lists_L.in."+str(case))
        SNF = _read_int_3D(gpath+"get_rotation_perms_lists_SNF.in."+str(case))
        Op = _read_fixOp_1D(gpath+"get_rotation_perms_lists_Op.in."+str(case))
        dperms = _read_RotPermList(gpath+"get_rotation_perms_lists_dperms.in."+str(case))
        RPlist = _read_RotPermList_1D(gpath+"get_rotation_perms_lists_RPlist.in."+str(case))
        eps = _read_float(gpath+"get_rotation_perms_lists_eps.in."+str(case))
        out1 = _get_rotation_perms_lists(A,HNF,L,SNF,Op,RPlist,dperms,eps)
        out2 = _read_RotPermList_1D(gpath+"get_rotation_perms_lists_RPlist.out."+str(case))
        for out in out1:
            out.perm.arrow_perm = None
        self._compare_outputs(out1,out2)

    def test_5(self):
        from phenum.grouptheory import _get_rotation_perms_lists
        case = 5
        A = list(map(list,zip(*_read_float_2D(gpath+"get_rotation_perms_lists_A.in."+str(case)))))
        HNF = _read_int_3D(gpath+"get_rotation_perms_lists_HNF.in."+str(case))
        L = _read_int_3D(gpath+"get_rotation_perms_lists_L.in."+str(case))
        SNF = _read_int_3D(gpath+"get_rotation_perms_lists_SNF.in."+str(case))
        Op = _read_fixOp_1D(gpath+"get_rotation_perms_lists_Op.in."+str(case))
        dperms = _read_RotPermList(gpath+"get_rotation_perms_lists_dperms.in."+str(case))
        RPlist = _read_RotPermList_1D(gpath+"get_rotation_perms_lists_RPlist.in."+str(case))
        eps = _read_float(gpath+"get_rotation_perms_lists_eps.in."+str(case))
        out1 = _get_rotation_perms_lists(A,HNF,L,SNF,Op,RPlist,dperms,eps)
        out2 = _read_RotPermList_1D(gpath+"get_rotation_perms_lists_RPlist.out."+str(case))
        for out in out1:
            out.perm.arrow_perm = None
        self._compare_outputs(out1,out2)

    def test_6(self):
        from phenum.grouptheory import _get_rotation_perms_lists
        case = 6
        A = list(map(list,zip(*_read_float_2D(gpath+"get_rotation_perms_lists_A.in."+str(case)))))
        HNF = _read_int_3D(gpath+"get_rotation_perms_lists_HNF.in."+str(case))
        L = _read_int_3D(gpath+"get_rotation_perms_lists_L.in."+str(case))
        SNF = _read_int_3D(gpath+"get_rotation_perms_lists_SNF.in."+str(case))
        Op = _read_fixOp_1D(gpath+"get_rotation_perms_lists_Op.in."+str(case))
        dperms = _read_RotPermList(gpath+"get_rotation_perms_lists_dperms.in."+str(case))
        RPlist = _read_RotPermList_1D(gpath+"get_rotation_perms_lists_RPlist.in."+str(case))
        eps = _read_float(gpath+"get_rotation_perms_lists_eps.in."+str(case))
        out1 = _get_rotation_perms_lists(A,HNF,L,SNF,Op,RPlist,dperms,eps)
        out2 = _read_RotPermList_1D(gpath+"get_rotation_perms_lists_RPlist.out."+str(case))
        for out in out1:
            out.perm.arrow_perm = None
        self._compare_outputs(out1,out2)
        
    def test_7(self):
        from phenum.grouptheory import _get_rotation_perms_lists
        case = 7
        A = list(map(list,zip(*_read_float_2D(gpath+"get_rotation_perms_lists_A.in."+str(case)))))
        HNF = _read_int_3D(gpath+"get_rotation_perms_lists_HNF.in."+str(case))
        L = _read_int_3D(gpath+"get_rotation_perms_lists_L.in."+str(case))
        SNF = _read_int_3D(gpath+"get_rotation_perms_lists_SNF.in."+str(case))
        Op = _read_fixOp_1D(gpath+"get_rotation_perms_lists_Op.in."+str(case))
        dperms = _read_RotPermList(gpath+"get_rotation_perms_lists_dperms.in."+str(case))
        RPlist = _read_RotPermList_1D(gpath+"get_rotation_perms_lists_RPlist.in."+str(case))
        eps = _read_float(gpath+"get_rotation_perms_lists_eps.in."+str(case))
        out1 = _get_rotation_perms_lists(A,HNF,L,SNF,Op,RPlist,dperms,eps)
        out2 = _read_RotPermList_1D(gpath+"get_rotation_perms_lists_RPlist.out."+str(case))
        for out in out1:
            out.perm.arrow_perm = None
        self._compare_outputs(out1,out2)

    def test_8(self):
        from phenum.grouptheory import _get_rotation_perms_lists
        case = 8
        A = list(map(list,zip(*_read_float_2D(gpath+"get_rotation_perms_lists_A.in."+str(case)))))
        HNF = _read_int_3D(gpath+"get_rotation_perms_lists_HNF.in."+str(case))
        L = _read_int_3D(gpath+"get_rotation_perms_lists_L.in."+str(case))
        SNF = _read_int_3D(gpath+"get_rotation_perms_lists_SNF.in."+str(case))
        Op = _read_fixOp_1D(gpath+"get_rotation_perms_lists_Op.in."+str(case))
        dperms = _read_RotPermList(gpath+"get_rotation_perms_lists_dperms.in."+str(case))
        RPlist = _read_RotPermList_1D(gpath+"get_rotation_perms_lists_RPlist.in."+str(case))
        eps = _read_float(gpath+"get_rotation_perms_lists_eps.in."+str(case))
        out1 = _get_rotation_perms_lists(A,HNF,L,SNF,Op,RPlist,dperms,eps)
        out2 = _read_RotPermList_1D(gpath+"get_rotation_perms_lists_RPlist.out."+str(case))
        for out in out1:
            out.perm.arrow_perm = None
        self._compare_outputs(out1,out2)

    def test_9(self):
        from phenum.grouptheory import _get_rotation_perms_lists
        case = 9
        A = list(map(list,zip(*_read_float_2D(gpath+"get_rotation_perms_lists_A.in."+str(case)))))
        HNF = _read_int_3D(gpath+"get_rotation_perms_lists_HNF.in."+str(case))
        L = _read_int_3D(gpath+"get_rotation_perms_lists_L.in."+str(case))
        SNF = _read_int_3D(gpath+"get_rotation_perms_lists_SNF.in."+str(case))
        Op = _read_fixOp_1D(gpath+"get_rotation_perms_lists_Op.in."+str(case))
        dperms = _read_RotPermList(gpath+"get_rotation_perms_lists_dperms.in."+str(case))
        RPlist = _read_RotPermList_1D(gpath+"get_rotation_perms_lists_RPlist.in."+str(case))
        eps = _read_float(gpath+"get_rotation_perms_lists_eps.in."+str(case))
        out1 = _get_rotation_perms_lists(A,HNF,L,SNF,Op,RPlist,dperms,eps)
        out2 = _read_RotPermList_1D(gpath+"get_rotation_perms_lists_RPlist.out."+str(case))
        for out in out1:
            out.perm.arrow_perm = None
        self._compare_outputs(out1,out2)

    def test_10(self):
        from phenum.grouptheory import _get_rotation_perms_lists
        case = 10
        A = list(map(list,zip(*_read_float_2D(gpath+"get_rotation_perms_lists_A.in."+str(case)))))
        HNF = _read_int_3D(gpath+"get_rotation_perms_lists_HNF.in."+str(case))
        L = _read_int_3D(gpath+"get_rotation_perms_lists_L.in."+str(case))
        SNF = _read_int_3D(gpath+"get_rotation_perms_lists_SNF.in."+str(case))
        Op = _read_fixOp_1D(gpath+"get_rotation_perms_lists_Op.in."+str(case))
        dperms = _read_RotPermList(gpath+"get_rotation_perms_lists_dperms.in."+str(case))
        RPlist = _read_RotPermList_1D(gpath+"get_rotation_perms_lists_RPlist.in."+str(case))
        eps = _read_float(gpath+"get_rotation_perms_lists_eps.in."+str(case))
        out1 = _get_rotation_perms_lists(A,HNF,L,SNF,Op,RPlist,dperms,eps)
        out2 = _read_RotPermList_1D(gpath+"get_rotation_perms_lists_RPlist.out."+str(case))
        for out in out1:
            out.perm.arrow_perm = None
        self._compare_outputs(out1,out2)

    def test_11(self):
        from phenum.grouptheory import _get_rotation_perms_lists
        case = 11
        A = list(map(list,zip(*_read_float_2D(gpath+"get_rotation_perms_lists_A.in."+str(case)))))
        HNF = _read_int_3D(gpath+"get_rotation_perms_lists_HNF.in."+str(case))
        L = _read_int_3D(gpath+"get_rotation_perms_lists_L.in."+str(case))
        SNF = _read_int_3D(gpath+"get_rotation_perms_lists_SNF.in."+str(case))
        Op = _read_fixOp_1D(gpath+"get_rotation_perms_lists_Op.in."+str(case))
        dperms = _read_RotPermList(gpath+"get_rotation_perms_lists_dperms.in."+str(case))
        RPlist = _read_RotPermList_1D(gpath+"get_rotation_perms_lists_RPlist.in."+str(case))
        eps = _read_float(gpath+"get_rotation_perms_lists_eps.in."+str(case))
        out1 = _get_rotation_perms_lists(A,HNF,L,SNF,Op,RPlist,dperms,eps)
        out2 = _read_RotPermList_1D(gpath+"get_rotation_perms_lists_RPlist.out."+str(case))
        for out in out1:
            out.perm.arrow_perm = None
        self._compare_outputs(out1,out2)

    def test_12(self):
        from phenum.grouptheory import _get_rotation_perms_lists
        case = 12
        A = list(map(list,zip(*_read_float_2D(gpath+"get_rotation_perms_lists_A.in."+str(case)))))
        HNF = _read_int_3D(gpath+"get_rotation_perms_lists_HNF.in."+str(case))
        L = _read_int_3D(gpath+"get_rotation_perms_lists_L.in."+str(case))
        SNF = _read_int_3D(gpath+"get_rotation_perms_lists_SNF.in."+str(case))
        Op = _read_fixOp_1D(gpath+"get_rotation_perms_lists_Op.in."+str(case))
        dperms = _read_RotPermList(gpath+"get_rotation_perms_lists_dperms.in."+str(case))
        RPlist = _read_RotPermList_1D(gpath+"get_rotation_perms_lists_RPlist.in."+str(case))
        eps = _read_float(gpath+"get_rotation_perms_lists_eps.in."+str(case))
        out1 = _get_rotation_perms_lists(A,HNF,L,SNF,Op,RPlist,dperms,eps)
        out2 = _read_RotPermList_1D(gpath+"get_rotation_perms_lists_RPlist.out."+str(case))
        for out in out1:
            out.perm.arrow_perm = None
        self._compare_outputs(out1,out2)

    def test_13(self):
        from phenum.grouptheory import _get_rotation_perms_lists
        case = 13
        A = list(map(list,zip(*_read_float_2D(gpath+"get_rotation_perms_lists_A.in."+str(case)))))
        HNF = _read_int_3D(gpath+"get_rotation_perms_lists_HNF.in."+str(case))
        L = _read_int_3D(gpath+"get_rotation_perms_lists_L.in."+str(case))
        SNF = _read_int_3D(gpath+"get_rotation_perms_lists_SNF.in."+str(case))
        Op = _read_fixOp_1D(gpath+"get_rotation_perms_lists_Op.in."+str(case))
        dperms = _read_RotPermList(gpath+"get_rotation_perms_lists_dperms.in."+str(case))
        RPlist = _read_RotPermList_1D(gpath+"get_rotation_perms_lists_RPlist.in."+str(case))
        eps = _read_float(gpath+"get_rotation_perms_lists_eps.in."+str(case))
        out1 = _get_rotation_perms_lists(A,HNF,L,SNF,Op,RPlist,dperms,eps)
        out2 = _read_RotPermList_1D(gpath+"get_rotation_perms_lists_RPlist.out."+str(case))
        for out in out1:
            out.perm.arrow_perm = None
        self._compare_outputs(out1,out2)

    def test_14(self):
        from phenum.grouptheory import _get_rotation_perms_lists
        case = 14
        A = list(map(list,zip(*_read_float_2D(gpath+"get_rotation_perms_lists_A.in."+str(case)))))
        HNF = _read_int_3D(gpath+"get_rotation_perms_lists_HNF.in."+str(case))
        L = _read_int_3D(gpath+"get_rotation_perms_lists_L.in."+str(case))
        SNF = _read_int_3D(gpath+"get_rotation_perms_lists_SNF.in."+str(case))
        Op = _read_fixOp_1D(gpath+"get_rotation_perms_lists_Op.in."+str(case))
        dperms = _read_RotPermList(gpath+"get_rotation_perms_lists_dperms.in."+str(case))
        RPlist = _read_RotPermList_1D(gpath+"get_rotation_perms_lists_RPlist.in."+str(case))
        eps = _read_float(gpath+"get_rotation_perms_lists_eps.in."+str(case))
        out1 = _get_rotation_perms_lists(A,HNF,L,SNF,Op,RPlist,dperms,eps)
        out2 = _read_RotPermList_1D(gpath+"get_rotation_perms_lists_RPlist.out."+str(case))
        for out in out1:
            out.perm.arrow_perm = None
        self._compare_outputs(out1,out2)

    def test_15(self):
        from phenum.grouptheory import _get_rotation_perms_lists
        case = 15
        A = list(map(list,zip(*_read_float_2D(gpath+"get_rotation_perms_lists_A.in."+str(case)))))
        HNF = _read_int_3D(gpath+"get_rotation_perms_lists_HNF.in."+str(case))
        L = _read_int_3D(gpath+"get_rotation_perms_lists_L.in."+str(case))
        SNF = _read_int_3D(gpath+"get_rotation_perms_lists_SNF.in."+str(case))
        Op = _read_fixOp_1D(gpath+"get_rotation_perms_lists_Op.in."+str(case))
        dperms = _read_RotPermList(gpath+"get_rotation_perms_lists_dperms.in."+str(case))
        RPlist = _read_RotPermList_1D(gpath+"get_rotation_perms_lists_RPlist.in."+str(case))
        eps = _read_float(gpath+"get_rotation_perms_lists_eps.in."+str(case))
        out1 = _get_rotation_perms_lists(A,HNF,L,SNF,Op,RPlist,dperms,eps)
        out2 = _read_RotPermList_1D(gpath+"get_rotation_perms_lists_RPlist.out."+str(case))
        for out in out1:
            out.perm.arrow_perm = None
        self._compare_outputs(out1,out2)

    def test_16(self):
        from phenum.grouptheory import _get_rotation_perms_lists
        case = 16
        A = list(map(list,zip(*_read_float_2D(gpath+"get_rotation_perms_lists_A.in."+str(case)))))
        HNF = _read_int_3D(gpath+"get_rotation_perms_lists_HNF.in."+str(case))
        L = _read_int_3D(gpath+"get_rotation_perms_lists_L.in."+str(case))
        SNF = _read_int_3D(gpath+"get_rotation_perms_lists_SNF.in."+str(case))
        Op = _read_fixOp_1D(gpath+"get_rotation_perms_lists_Op.in."+str(case))
        dperms = _read_RotPermList(gpath+"get_rotation_perms_lists_dperms.in."+str(case))
        RPlist = _read_RotPermList_1D(gpath+"get_rotation_perms_lists_RPlist.in."+str(case))
        eps = _read_float(gpath+"get_rotation_perms_lists_eps.in."+str(case))
        out1 = _get_rotation_perms_lists(A,HNF,L,SNF,Op,RPlist,dperms,eps)
        out2 = _read_RotPermList_1D(gpath+"get_rotation_perms_lists_RPlist.out."+str(case))
        for out in out1:
            out.perm.arrow_perm = None
        self._compare_outputs(out1,out2)

    def test_17(self):
        from phenum.grouptheory import _get_rotation_perms_lists
        case = 17
        A = list(map(list,zip(*_read_float_2D(gpath+"get_rotation_perms_lists_A.in."+str(case)))))
        HNF = _read_int_3D(gpath+"get_rotation_perms_lists_HNF.in."+str(case))
        L = _read_int_3D(gpath+"get_rotation_perms_lists_L.in."+str(case))
        SNF = _read_int_3D(gpath+"get_rotation_perms_lists_SNF.in."+str(case))
        Op = _read_fixOp_1D(gpath+"get_rotation_perms_lists_Op.in."+str(case))
        dperms = _read_RotPermList(gpath+"get_rotation_perms_lists_dperms.in."+str(case))
        RPlist = _read_RotPermList_1D(gpath+"get_rotation_perms_lists_RPlist.in."+str(case))
        eps = _read_float(gpath+"get_rotation_perms_lists_eps.in."+str(case))
        out1 = _get_rotation_perms_lists(A,HNF,L,SNF,Op,RPlist,dperms,eps)
        out2 = _read_RotPermList_1D(gpath+"get_rotation_perms_lists_RPlist.out."+str(case))
        for out in out1:
            out.perm.arrow_perm = None
        self._compare_outputs(out1,out2)

    def test_18(self):
        from phenum.grouptheory import _get_rotation_perms_lists
        case = 18
        A = list(map(list,zip(*_read_float_2D(gpath+"get_rotation_perms_lists_A.in."+str(case)))))
        HNF = _read_int_3D(gpath+"get_rotation_perms_lists_HNF.in."+str(case))
        L = _read_int_3D(gpath+"get_rotation_perms_lists_L.in."+str(case))
        SNF = _read_int_3D(gpath+"get_rotation_perms_lists_SNF.in."+str(case))
        Op = _read_fixOp_1D(gpath+"get_rotation_perms_lists_Op.in."+str(case))
        dperms = _read_RotPermList(gpath+"get_rotation_perms_lists_dperms.in."+str(case))
        RPlist = _read_RotPermList_1D(gpath+"get_rotation_perms_lists_RPlist.in."+str(case))
        eps = _read_float(gpath+"get_rotation_perms_lists_eps.in."+str(case))
        out1 = _get_rotation_perms_lists(A,HNF,L,SNF,Op,RPlist,dperms,eps)
        out2 = _read_RotPermList_1D(gpath+"get_rotation_perms_lists_RPlist.out."+str(case))
        for out in out1:
            out.perm.arrow_perm = None
        self._compare_outputs(out1,out2)

    def test_19(self):
        from phenum.grouptheory import _get_rotation_perms_lists
        case = 19
        A = list(map(list,zip(*_read_float_2D(gpath+"get_rotation_perms_lists_A.in."+str(case)))))
        HNF = _read_int_3D(gpath+"get_rotation_perms_lists_HNF.in."+str(case))
        L = _read_int_3D(gpath+"get_rotation_perms_lists_L.in."+str(case))
        SNF = _read_int_3D(gpath+"get_rotation_perms_lists_SNF.in."+str(case))
        Op = _read_fixOp_1D(gpath+"get_rotation_perms_lists_Op.in."+str(case))
        dperms = _read_RotPermList(gpath+"get_rotation_perms_lists_dperms.in."+str(case))
        RPlist = _read_RotPermList_1D(gpath+"get_rotation_perms_lists_RPlist.in."+str(case))
        eps = _read_float(gpath+"get_rotation_perms_lists_eps.in."+str(case))
        out1 = _get_rotation_perms_lists(A,HNF,L,SNF,Op,RPlist,dperms,eps)
        out2 = _read_RotPermList_1D(gpath+"get_rotation_perms_lists_RPlist.out."+str(case))
        for out in out1:
            out.perm.arrow_perm = None
        self._compare_outputs(out1,out2)

    def test_20(self):
        from phenum.grouptheory import _get_rotation_perms_lists
        case = 20
        A = list(map(list,zip(*_read_float_2D(gpath+"get_rotation_perms_lists_A.in."+str(case)))))
        HNF = _read_int_3D(gpath+"get_rotation_perms_lists_HNF.in."+str(case))
        L = _read_int_3D(gpath+"get_rotation_perms_lists_L.in."+str(case))
        SNF = _read_int_3D(gpath+"get_rotation_perms_lists_SNF.in."+str(case))
        Op = _read_fixOp_1D(gpath+"get_rotation_perms_lists_Op.in."+str(case))
        dperms = _read_RotPermList(gpath+"get_rotation_perms_lists_dperms.in."+str(case))
        RPlist = _read_RotPermList_1D(gpath+"get_rotation_perms_lists_RPlist.in."+str(case))
        eps = _read_float(gpath+"get_rotation_perms_lists_eps.in."+str(case))
        out1 = _get_rotation_perms_lists(A,HNF,L,SNF,Op,RPlist,dperms,eps)
        out2 = _read_RotPermList_1D(gpath+"get_rotation_perms_lists_RPlist.out."+str(case))
        for out in out1:
            out.perm.arrow_perm = None
        self._compare_outputs(out1,out2)

    def test_21(self):
        from phenum.grouptheory import _get_rotation_perms_lists
        case = 21
        A = list(map(list,zip(*_read_float_2D(gpath+"get_rotation_perms_lists_A.in."+str(case)))))
        HNF = _read_int_3D(gpath+"get_rotation_perms_lists_HNF.in."+str(case))
        L = _read_int_3D(gpath+"get_rotation_perms_lists_L.in."+str(case))
        SNF = _read_int_3D(gpath+"get_rotation_perms_lists_SNF.in."+str(case))
        Op = _read_fixOp_1D(gpath+"get_rotation_perms_lists_Op.in."+str(case))
        dperms = _read_RotPermList(gpath+"get_rotation_perms_lists_dperms.in."+str(case))
        RPlist = _read_RotPermList_1D(gpath+"get_rotation_perms_lists_RPlist.in."+str(case))
        eps = _read_float(gpath+"get_rotation_perms_lists_eps.in."+str(case))
        out1 = _get_rotation_perms_lists(A,HNF,L,SNF,Op,RPlist,dperms,eps)
        out2 = _read_RotPermList_1D(gpath+"get_rotation_perms_lists_RPlist.out."+str(case))
        for out in out1:
            out.perm.arrow_perm = None
        self._compare_outputs(out1,out2)

    def test_22(self):
        from phenum.grouptheory import _get_rotation_perms_lists
        case = 22
        A = list(map(list,zip(*_read_float_2D(gpath+"get_rotation_perms_lists_A.in."+str(case)))))
        HNF = _read_int_3D(gpath+"get_rotation_perms_lists_HNF.in."+str(case))
        L = _read_int_3D(gpath+"get_rotation_perms_lists_L.in."+str(case))
        SNF = _read_int_3D(gpath+"get_rotation_perms_lists_SNF.in."+str(case))
        Op = _read_fixOp_1D(gpath+"get_rotation_perms_lists_Op.in."+str(case))
        dperms = _read_RotPermList(gpath+"get_rotation_perms_lists_dperms.in."+str(case))
        RPlist = _read_RotPermList_1D(gpath+"get_rotation_perms_lists_RPlist.in."+str(case))
        eps = _read_float(gpath+"get_rotation_perms_lists_eps.in."+str(case))
        out1 = _get_rotation_perms_lists(A,HNF,L,SNF,Op,RPlist,dperms,eps)
        out2 = _read_RotPermList_1D(gpath+"get_rotation_perms_lists_RPlist.out."+str(case))
        for out in out1:
            out.perm.arrow_perm = None
        self._compare_outputs(out1,out2)

    def test_23(self):
        from phenum.grouptheory import _get_rotation_perms_lists
        case = 23
        A = list(map(list,zip(*_read_float_2D(gpath+"get_rotation_perms_lists_A.in."+str(case)))))
        HNF = _read_int_3D(gpath+"get_rotation_perms_lists_HNF.in."+str(case))
        L = _read_int_3D(gpath+"get_rotation_perms_lists_L.in."+str(case))
        SNF = _read_int_3D(gpath+"get_rotation_perms_lists_SNF.in."+str(case))
        Op = _read_fixOp_1D(gpath+"get_rotation_perms_lists_Op.in."+str(case))
        dperms = _read_RotPermList(gpath+"get_rotation_perms_lists_dperms.in."+str(case))
        RPlist = _read_RotPermList_1D(gpath+"get_rotation_perms_lists_RPlist.in."+str(case))
        eps = _read_float(gpath+"get_rotation_perms_lists_eps.in."+str(case))
        out1 = _get_rotation_perms_lists(A,HNF,L,SNF,Op,RPlist,dperms,eps)
        out2 = _read_RotPermList_1D(gpath+"get_rotation_perms_lists_RPlist.out."+str(case))
        for out in out1:
            out.perm.arrow_perm = None
        self._compare_outputs(out1,out2)

    def test_24(self):
        from phenum.grouptheory import _get_rotation_perms_lists
        case = 24
        A = list(map(list,zip(*_read_float_2D(gpath+"get_rotation_perms_lists_A.in."+str(case)))))
        HNF = _read_int_3D(gpath+"get_rotation_perms_lists_HNF.in."+str(case))
        L = _read_int_3D(gpath+"get_rotation_perms_lists_L.in."+str(case))
        SNF = _read_int_3D(gpath+"get_rotation_perms_lists_SNF.in."+str(case))
        Op = _read_fixOp_1D(gpath+"get_rotation_perms_lists_Op.in."+str(case))
        dperms = _read_RotPermList(gpath+"get_rotation_perms_lists_dperms.in."+str(case))
        RPlist = _read_RotPermList_1D(gpath+"get_rotation_perms_lists_RPlist.in."+str(case))
        eps = _read_float(gpath+"get_rotation_perms_lists_eps.in."+str(case))
        out1 = _get_rotation_perms_lists(A,HNF,L,SNF,Op,RPlist,dperms,eps)
        out2 = _read_RotPermList_1D(gpath+"get_rotation_perms_lists_RPlist.out."+str(case))
        for out in out1:
            out.perm.arrow_perm = None
        self._compare_outputs(out1,out2)

    def test_25(self):
        from phenum.grouptheory import _get_rotation_perms_lists
        case = 25
        A = list(map(list,zip(*_read_float_2D(gpath+"get_rotation_perms_lists_A.in."+str(case)))))
        HNF = _read_int_3D(gpath+"get_rotation_perms_lists_HNF.in."+str(case))
        L = _read_int_3D(gpath+"get_rotation_perms_lists_L.in."+str(case))
        SNF = _read_int_3D(gpath+"get_rotation_perms_lists_SNF.in."+str(case))
        Op = _read_fixOp_1D(gpath+"get_rotation_perms_lists_Op.in."+str(case))
        dperms = _read_RotPermList(gpath+"get_rotation_perms_lists_dperms.in."+str(case))
        RPlist = _read_RotPermList_1D(gpath+"get_rotation_perms_lists_RPlist.in."+str(case))
        eps = _read_float(gpath+"get_rotation_perms_lists_eps.in."+str(case))
        out1 = _get_rotation_perms_lists(A,HNF,L,SNF,Op,RPlist,dperms,eps)
        out2 = _read_RotPermList_1D(gpath+"get_rotation_perms_lists_RPlist.out."+str(case))
        for out in out1:
            out.perm.arrow_perm = None
        self._compare_outputs(out1,out2)

    def test_26(self):
        from phenum.grouptheory import _get_rotation_perms_lists
        case = 26
        A = list(map(list,zip(*_read_float_2D(gpath+"get_rotation_perms_lists_A.in."+str(case)))))
        HNF = _read_int_3D(gpath+"get_rotation_perms_lists_HNF.in."+str(case))
        L = _read_int_3D(gpath+"get_rotation_perms_lists_L.in."+str(case))
        SNF = _read_int_3D(gpath+"get_rotation_perms_lists_SNF.in."+str(case))
        Op = _read_fixOp_1D(gpath+"get_rotation_perms_lists_Op.in."+str(case))
        dperms = _read_RotPermList(gpath+"get_rotation_perms_lists_dperms.in."+str(case))
        RPlist = _read_RotPermList_1D(gpath+"get_rotation_perms_lists_RPlist.in."+str(case))
        eps = _read_float(gpath+"get_rotation_perms_lists_eps.in."+str(case))
        out1 = _get_rotation_perms_lists(A,HNF,L,SNF,Op,RPlist,dperms,eps)
        out2 = _read_RotPermList_1D(gpath+"get_rotation_perms_lists_RPlist.out."+str(case))
        for out in out1:
            out.perm.arrow_perm = None
        self._compare_outputs(out1,out2)

    def test_27(self):
        from phenum.grouptheory import _get_rotation_perms_lists
        case = 27
        A = list(map(list,zip(*_read_float_2D(gpath+"get_rotation_perms_lists_A.in."+str(case)))))
        HNF = _read_int_3D(gpath+"get_rotation_perms_lists_HNF.in."+str(case))
        L = _read_int_3D(gpath+"get_rotation_perms_lists_L.in."+str(case))
        SNF = _read_int_3D(gpath+"get_rotation_perms_lists_SNF.in."+str(case))
        Op = _read_fixOp_1D(gpath+"get_rotation_perms_lists_Op.in."+str(case))
        dperms = _read_RotPermList(gpath+"get_rotation_perms_lists_dperms.in."+str(case))
        RPlist = _read_RotPermList_1D(gpath+"get_rotation_perms_lists_RPlist.in."+str(case))
        eps = _read_float(gpath+"get_rotation_perms_lists_eps.in."+str(case))
        out1 = _get_rotation_perms_lists(A,HNF,L,SNF,Op,RPlist,dperms,eps)
        out2 = _read_RotPermList_1D(gpath+"get_rotation_perms_lists_RPlist.out."+str(case))
        for out in out1:
            out.perm.arrow_perm = None
        self._compare_outputs(out1,out2)

    def test_28(self):
        from phenum.grouptheory import _get_rotation_perms_lists
        case = 28
        A = list(map(list,zip(*_read_float_2D(gpath+"get_rotation_perms_lists_A.in."+str(case)))))
        HNF = _read_int_3D(gpath+"get_rotation_perms_lists_HNF.in."+str(case))
        L = _read_int_3D(gpath+"get_rotation_perms_lists_L.in."+str(case))
        SNF = _read_int_3D(gpath+"get_rotation_perms_lists_SNF.in."+str(case))
        Op = _read_fixOp_1D(gpath+"get_rotation_perms_lists_Op.in."+str(case))
        dperms = _read_RotPermList(gpath+"get_rotation_perms_lists_dperms.in."+str(case))
        RPlist = _read_RotPermList_1D(gpath+"get_rotation_perms_lists_RPlist.in."+str(case))
        eps = _read_float(gpath+"get_rotation_perms_lists_eps.in."+str(case))
        out1 = _get_rotation_perms_lists(A,HNF,L,SNF,Op,RPlist,dperms,eps)
        out2 = _read_RotPermList_1D(gpath+"get_rotation_perms_lists_RPlist.out."+str(case))
        for out in out1:
            out.perm.arrow_perm = None
        self._compare_outputs(out1,out2)

    def test_29(self):
        from phenum.grouptheory import _get_rotation_perms_lists
        case = 29
        A = list(map(list,zip(*_read_float_2D(gpath+"get_rotation_perms_lists_A.in."+str(case)))))
        HNF = _read_int_3D(gpath+"get_rotation_perms_lists_HNF.in."+str(case))
        L = _read_int_3D(gpath+"get_rotation_perms_lists_L.in."+str(case))
        SNF = _read_int_3D(gpath+"get_rotation_perms_lists_SNF.in."+str(case))
        Op = _read_fixOp_1D(gpath+"get_rotation_perms_lists_Op.in."+str(case))
        dperms = _read_RotPermList(gpath+"get_rotation_perms_lists_dperms.in."+str(case))
        RPlist = _read_RotPermList_1D(gpath+"get_rotation_perms_lists_RPlist.in."+str(case))
        eps = _read_float(gpath+"get_rotation_perms_lists_eps.in."+str(case))
        out1 = _get_rotation_perms_lists(A,HNF,L,SNF,Op,RPlist,dperms,eps)
        out2 = _read_RotPermList_1D(gpath+"get_rotation_perms_lists_RPlist.out."+str(case))
        for out in out1:
            out.perm.arrow_perm = None
        self._compare_outputs(out1,out2)

    def test_30(self):
        from phenum.grouptheory import _get_rotation_perms_lists
        case = 30
        A = list(map(list,zip(*_read_float_2D(gpath+"get_rotation_perms_lists_A.in."+str(case)))))
        HNF = _read_int_3D(gpath+"get_rotation_perms_lists_HNF.in."+str(case))
        L = _read_int_3D(gpath+"get_rotation_perms_lists_L.in."+str(case))
        SNF = _read_int_3D(gpath+"get_rotation_perms_lists_SNF.in."+str(case))
        Op = _read_fixOp_1D(gpath+"get_rotation_perms_lists_Op.in."+str(case))
        dperms = _read_RotPermList(gpath+"get_rotation_perms_lists_dperms.in."+str(case))
        RPlist = _read_RotPermList_1D(gpath+"get_rotation_perms_lists_RPlist.in."+str(case))
        eps = _read_float(gpath+"get_rotation_perms_lists_eps.in."+str(case))
        out1 = _get_rotation_perms_lists(A,HNF,L,SNF,Op,RPlist,dperms,eps)
        out2 = _read_RotPermList_1D(gpath+"get_rotation_perms_lists_RPlist.out."+str(case))
        for out in out1:
            out.perm.arrow_perm = None
        self._compare_outputs(out1,out2)

    def test_31(self):
        from phenum.grouptheory import _get_rotation_perms_lists
        case = 31
        A = list(map(list,zip(*_read_float_2D(gpath+"get_rotation_perms_lists_A.in."+str(case)))))
        HNF = _read_int_3D(gpath+"get_rotation_perms_lists_HNF.in."+str(case))
        L = _read_int_3D(gpath+"get_rotation_perms_lists_L.in."+str(case))
        SNF = _read_int_3D(gpath+"get_rotation_perms_lists_SNF.in."+str(case))
        Op = _read_fixOp_1D(gpath+"get_rotation_perms_lists_Op.in."+str(case))
        dperms = _read_RotPermList(gpath+"get_rotation_perms_lists_dperms.in."+str(case))
        RPlist = _read_RotPermList_1D(gpath+"get_rotation_perms_lists_RPlist.in."+str(case))
        eps = _read_float(gpath+"get_rotation_perms_lists_eps.in."+str(case))
        out1 = _get_rotation_perms_lists(A,HNF,L,SNF,Op,RPlist,dperms,eps)
        out2 = _read_RotPermList_1D(gpath+"get_rotation_perms_lists_RPlist.out."+str(case))
        for out in out1:
            out.perm.arrow_perm = None
        self._compare_outputs(out1,out2)

    def test_32(self):
        from phenum.grouptheory import _get_rotation_perms_lists
        case = 32
        A = list(map(list,zip(*_read_float_2D(gpath+"get_rotation_perms_lists_A.in."+str(case)))))
        HNF = _read_int_3D(gpath+"get_rotation_perms_lists_HNF.in."+str(case))
        L = _read_int_3D(gpath+"get_rotation_perms_lists_L.in."+str(case))
        SNF = _read_int_3D(gpath+"get_rotation_perms_lists_SNF.in."+str(case))
        Op = _read_fixOp_1D(gpath+"get_rotation_perms_lists_Op.in."+str(case))
        dperms = _read_RotPermList(gpath+"get_rotation_perms_lists_dperms.in."+str(case))
        RPlist = _read_RotPermList_1D(gpath+"get_rotation_perms_lists_RPlist.in."+str(case))
        eps = _read_float(gpath+"get_rotation_perms_lists_eps.in."+str(case))
        out1 = _get_rotation_perms_lists(A,HNF,L,SNF,Op,RPlist,dperms,eps)
        out2 = _read_RotPermList_1D(gpath+"get_rotation_perms_lists_RPlist.out."+str(case))
        for out in out1:
            out.perm.arrow_perm = None
        self._compare_outputs(out1,out2)

    def test_33(self):
        from phenum.grouptheory import _get_rotation_perms_lists
        case = 33
        A = list(map(list,zip(*_read_float_2D(gpath+"get_rotation_perms_lists_A.in."+str(case)))))
        HNF = _read_int_3D(gpath+"get_rotation_perms_lists_HNF.in."+str(case))
        L = _read_int_3D(gpath+"get_rotation_perms_lists_L.in."+str(case))
        SNF = _read_int_3D(gpath+"get_rotation_perms_lists_SNF.in."+str(case))
        Op = _read_fixOp_1D(gpath+"get_rotation_perms_lists_Op.in."+str(case))
        dperms = _read_RotPermList(gpath+"get_rotation_perms_lists_dperms.in."+str(case))
        RPlist = _read_RotPermList_1D(gpath+"get_rotation_perms_lists_RPlist.in."+str(case))
        eps = _read_float(gpath+"get_rotation_perms_lists_eps.in."+str(case))
        out1 = _get_rotation_perms_lists(A,HNF,L,SNF,Op,RPlist,dperms,eps)
        out2 = _read_RotPermList_1D(gpath+"get_rotation_perms_lists_RPlist.out."+str(case))
        for out in out1:
            out.perm.arrow_perm = None
        self._compare_outputs(out1,out2)

    def test_34(self):
        from phenum.grouptheory import _get_rotation_perms_lists
        case = 34
        A = list(map(list,zip(*_read_float_2D(gpath+"get_rotation_perms_lists_A.in."+str(case)))))
        HNF = _read_int_3D(gpath+"get_rotation_perms_lists_HNF.in."+str(case))
        L = _read_int_3D(gpath+"get_rotation_perms_lists_L.in."+str(case))
        SNF = _read_int_3D(gpath+"get_rotation_perms_lists_SNF.in."+str(case))
        Op = _read_fixOp_1D(gpath+"get_rotation_perms_lists_Op.in."+str(case))
        dperms = _read_RotPermList(gpath+"get_rotation_perms_lists_dperms.in."+str(case))
        RPlist = _read_RotPermList_1D(gpath+"get_rotation_perms_lists_RPlist.in."+str(case))
        eps = _read_float(gpath+"get_rotation_perms_lists_eps.in."+str(case))
        out1 = _get_rotation_perms_lists(A,HNF,L,SNF,Op,RPlist,dperms,eps)
        out2 = _read_RotPermList_1D(gpath+"get_rotation_perms_lists_RPlist.out."+str(case))
        for out in out1:
            out.perm.arrow_perm = None
        self._compare_outputs(out1,out2)

    def test_35(self):
        from phenum.grouptheory import _get_rotation_perms_lists
        case = 35
        A = list(map(list,zip(*_read_float_2D(gpath+"get_rotation_perms_lists_A.in."+str(case)))))
        HNF = _read_int_3D(gpath+"get_rotation_perms_lists_HNF.in."+str(case))
        L = _read_int_3D(gpath+"get_rotation_perms_lists_L.in."+str(case))
        SNF = _read_int_3D(gpath+"get_rotation_perms_lists_SNF.in."+str(case))
        Op = _read_fixOp_1D(gpath+"get_rotation_perms_lists_Op.in."+str(case))
        dperms = _read_RotPermList(gpath+"get_rotation_perms_lists_dperms.in."+str(case))
        RPlist = _read_RotPermList_1D(gpath+"get_rotation_perms_lists_RPlist.in."+str(case))
        eps = _read_float(gpath+"get_rotation_perms_lists_eps.in."+str(case))
        out1 = _get_rotation_perms_lists(A,HNF,L,SNF,Op,RPlist,dperms,eps)
        out2 = _read_RotPermList_1D(gpath+"get_rotation_perms_lists_RPlist.out."+str(case))
        for out in out1:
            out.perm.arrow_perm = None
        self._compare_outputs(out1,out2)

    def test_36(self):
        from phenum.grouptheory import _get_rotation_perms_lists
        case = 36
        A = list(map(list,zip(*_read_float_2D(gpath+"get_rotation_perms_lists_A.in."+str(case)))))
        HNF = _read_int_3D(gpath+"get_rotation_perms_lists_HNF.in."+str(case))
        L = _read_int_3D(gpath+"get_rotation_perms_lists_L.in."+str(case))
        SNF = _read_int_3D(gpath+"get_rotation_perms_lists_SNF.in."+str(case))
        Op = _read_fixOp_1D(gpath+"get_rotation_perms_lists_Op.in."+str(case))
        dperms = _read_RotPermList(gpath+"get_rotation_perms_lists_dperms.in."+str(case))
        RPlist = _read_RotPermList_1D(gpath+"get_rotation_perms_lists_RPlist.in."+str(case))
        eps = _read_float(gpath+"get_rotation_perms_lists_eps.in."+str(case))
        out1 = _get_rotation_perms_lists(A,HNF,L,SNF,Op,RPlist,dperms,eps)
        out2 = _read_RotPermList_1D(gpath+"get_rotation_perms_lists_RPlist.out."+str(case))
        for out in out1:
            out.perm.arrow_perm = None
        self._compare_outputs(out1,out2)

    def test_37(self):
        from phenum.grouptheory import _get_rotation_perms_lists
        case = 37
        A = list(map(list,zip(*_read_float_2D(gpath+"get_rotation_perms_lists_A.in."+str(case)))))
        HNF = _read_int_3D(gpath+"get_rotation_perms_lists_HNF.in."+str(case))
        L = _read_int_3D(gpath+"get_rotation_perms_lists_L.in."+str(case))
        SNF = _read_int_3D(gpath+"get_rotation_perms_lists_SNF.in."+str(case))
        Op = _read_fixOp_1D(gpath+"get_rotation_perms_lists_Op.in."+str(case))
        dperms = _read_RotPermList(gpath+"get_rotation_perms_lists_dperms.in."+str(case))
        RPlist = _read_RotPermList_1D(gpath+"get_rotation_perms_lists_RPlist.in."+str(case))
        eps = _read_float(gpath+"get_rotation_perms_lists_eps.in."+str(case))
        out1 = _get_rotation_perms_lists(A,HNF,L,SNF,Op,RPlist,dperms,eps)
        out2 = _read_RotPermList_1D(gpath+"get_rotation_perms_lists_RPlist.out."+str(case))
        for out in out1:
            out.perm.arrow_perm = None
        self._compare_outputs(out1,out2)

    def test_38(self):
        from phenum.grouptheory import _get_rotation_perms_lists
        case = 38
        A = list(map(list,zip(*_read_float_2D(gpath+"get_rotation_perms_lists_A.in."+str(case)))))
        HNF = _read_int_3D(gpath+"get_rotation_perms_lists_HNF.in."+str(case))
        L = _read_int_3D(gpath+"get_rotation_perms_lists_L.in."+str(case))
        SNF = _read_int_3D(gpath+"get_rotation_perms_lists_SNF.in."+str(case))
        Op = _read_fixOp_1D(gpath+"get_rotation_perms_lists_Op.in."+str(case))
        dperms = _read_RotPermList(gpath+"get_rotation_perms_lists_dperms.in."+str(case))
        RPlist = _read_RotPermList_1D(gpath+"get_rotation_perms_lists_RPlist.in."+str(case))
        eps = _read_float(gpath+"get_rotation_perms_lists_eps.in."+str(case))
        out1 = _get_rotation_perms_lists(A,HNF,L,SNF,Op,RPlist,dperms,eps)
        out2 = _read_RotPermList_1D(gpath+"get_rotation_perms_lists_RPlist.out."+str(case))
        for out in out1:
            out.perm.arrow_perm = None
        self._compare_outputs(out1,out2)

    def test_39(self):
        from phenum.grouptheory import _get_rotation_perms_lists
        case = 39
        A = list(map(list,zip(*_read_float_2D(gpath+"get_rotation_perms_lists_A.in."+str(case)))))
        HNF = _read_int_3D(gpath+"get_rotation_perms_lists_HNF.in."+str(case))
        L = _read_int_3D(gpath+"get_rotation_perms_lists_L.in."+str(case))
        SNF = _read_int_3D(gpath+"get_rotation_perms_lists_SNF.in."+str(case))
        Op = _read_fixOp_1D(gpath+"get_rotation_perms_lists_Op.in."+str(case))
        dperms = _read_RotPermList(gpath+"get_rotation_perms_lists_dperms.in."+str(case))
        RPlist = _read_RotPermList_1D(gpath+"get_rotation_perms_lists_RPlist.in."+str(case))
        eps = _read_float(gpath+"get_rotation_perms_lists_eps.in."+str(case))
        out1 = _get_rotation_perms_lists(A,HNF,L,SNF,Op,RPlist,dperms,eps)
        out2 = _read_RotPermList_1D(gpath+"get_rotation_perms_lists_RPlist.out."+str(case))
        for out in out1:
            out.perm.arrow_perm = None
        self._compare_outputs(out1,out2)

    def test_50(self):
        from phenum.grouptheory import _get_rotation_perms_lists
        case = 50
        A = list(map(list,zip(*_read_float_2D(gpath+"get_rotation_perms_lists_A.in."+str(case)))))
        HNF = _read_int_3D(gpath+"get_rotation_perms_lists_HNF.in."+str(case))
        L = _read_int_3D(gpath+"get_rotation_perms_lists_L.in."+str(case))
        SNF = _read_int_3D(gpath+"get_rotation_perms_lists_SNF.in."+str(case))
        Op = _read_fixOp_1D(gpath+"get_rotation_perms_lists_Op.in."+str(case))
        dperms = _read_RotPermList(gpath+"get_rotation_perms_lists_dperms.in."+str(case))
        RPlist = _read_RotPermList_1D(gpath+"get_rotation_perms_lists_RPlist.in."+str(case))
        eps = _read_float(gpath+"get_rotation_perms_lists_eps.in."+str(case))
        out1 = _get_rotation_perms_lists(A,HNF,L,SNF,Op,RPlist,dperms,eps)
        out2 = _read_RotPermList_1D(gpath+"get_rotation_perms_lists_RPlist.out."+str(case))
        for out in out1:
            out.perm.arrow_perm = None
        self._compare_outputs(out1,out2)

    def test_51(self):
        from phenum.grouptheory import _get_rotation_perms_lists
        case = 51
        A = list(map(list,zip(*_read_float_2D(gpath+"get_rotation_perms_lists_A.in."+str(case)))))
        HNF = _read_int_3D(gpath+"get_rotation_perms_lists_HNF.in."+str(case))
        L = _read_int_3D(gpath+"get_rotation_perms_lists_L.in."+str(case))
        SNF = _read_int_3D(gpath+"get_rotation_perms_lists_SNF.in."+str(case))
        Op = _read_fixOp_1D(gpath+"get_rotation_perms_lists_Op.in."+str(case))
        dperms = _read_RotPermList(gpath+"get_rotation_perms_lists_dperms.in."+str(case))
        RPlist = _read_RotPermList_1D(gpath+"get_rotation_perms_lists_RPlist.in."+str(case))
        eps = _read_float(gpath+"get_rotation_perms_lists_eps.in."+str(case))
        out1 = _get_rotation_perms_lists(A,HNF,L,SNF,Op,RPlist,dperms,eps)
        out2 = _read_RotPermList_1D(gpath+"get_rotation_perms_lists_RPlist.out."+str(case))
        for out in out1:
            out.perm.arrow_perm = None
        self._compare_outputs(out1,out2)

    def test_52(self):
        from phenum.grouptheory import _get_rotation_perms_lists
        case = 52
        A = list(map(list,zip(*_read_float_2D(gpath+"get_rotation_perms_lists_A.in."+str(case)))))
        HNF = _read_int_3D(gpath+"get_rotation_perms_lists_HNF.in."+str(case))
        L = _read_int_3D(gpath+"get_rotation_perms_lists_L.in."+str(case))
        SNF = _read_int_3D(gpath+"get_rotation_perms_lists_SNF.in."+str(case))
        Op = _read_fixOp_1D(gpath+"get_rotation_perms_lists_Op.in."+str(case))
        dperms = _read_RotPermList(gpath+"get_rotation_perms_lists_dperms.in."+str(case))
        RPlist = _read_RotPermList_1D(gpath+"get_rotation_perms_lists_RPlist.in."+str(case))
        eps = _read_float(gpath+"get_rotation_perms_lists_eps.in."+str(case))
        out1 = _get_rotation_perms_lists(A,HNF,L,SNF,Op,RPlist,dperms,eps)
        out2 = _read_RotPermList_1D(gpath+"get_rotation_perms_lists_RPlist.out."+str(case))
        for out in out1:
            out.perm.arrow_perm = None
        self._compare_outputs(out1,out2)

    def test_53(self):
        from phenum.grouptheory import _get_rotation_perms_lists
        case = 53
        A = list(map(list,zip(*_read_float_2D(gpath+"get_rotation_perms_lists_A.in."+str(case)))))
        HNF = _read_int_3D(gpath+"get_rotation_perms_lists_HNF.in."+str(case))
        L = _read_int_3D(gpath+"get_rotation_perms_lists_L.in."+str(case))
        SNF = _read_int_3D(gpath+"get_rotation_perms_lists_SNF.in."+str(case))
        Op = _read_fixOp_1D(gpath+"get_rotation_perms_lists_Op.in."+str(case))
        dperms = _read_RotPermList(gpath+"get_rotation_perms_lists_dperms.in."+str(case))
        RPlist = _read_RotPermList_1D(gpath+"get_rotation_perms_lists_RPlist.in."+str(case))
        eps = _read_float(gpath+"get_rotation_perms_lists_eps.in."+str(case))
        out1 = _get_rotation_perms_lists(A,HNF,L,SNF,Op,RPlist,dperms,eps)
        out2 = _read_RotPermList_1D(gpath+"get_rotation_perms_lists_RPlist.out."+str(case))
        for out in out1:
            out.perm.arrow_perm = None
        self._compare_outputs(out1,out2)

    def test_54(self):
        from phenum.grouptheory import _get_rotation_perms_lists
        case = 54
        A = list(map(list,zip(*_read_float_2D(gpath+"get_rotation_perms_lists_A.in."+str(case)))))
        HNF = _read_int_3D(gpath+"get_rotation_perms_lists_HNF.in."+str(case))
        L = _read_int_3D(gpath+"get_rotation_perms_lists_L.in."+str(case))
        SNF = _read_int_3D(gpath+"get_rotation_perms_lists_SNF.in."+str(case))
        Op = _read_fixOp_1D(gpath+"get_rotation_perms_lists_Op.in."+str(case))
        dperms = _read_RotPermList(gpath+"get_rotation_perms_lists_dperms.in."+str(case))
        RPlist = _read_RotPermList_1D(gpath+"get_rotation_perms_lists_RPlist.in."+str(case))
        eps = _read_float(gpath+"get_rotation_perms_lists_eps.in."+str(case))
        out1 = _get_rotation_perms_lists(A,HNF,L,SNF,Op,RPlist,dperms,eps)
        out2 = _read_RotPermList_1D(gpath+"get_rotation_perms_lists_RPlist.out."+str(case))
        for out in out1:
            out.perm.arrow_perm = None
        self._compare_outputs(out1,out2)

    def test_55(self):
        from phenum.grouptheory import _get_rotation_perms_lists
        case = 55
        A = list(map(list,zip(*_read_float_2D(gpath+"get_rotation_perms_lists_A.in."+str(case)))))
        HNF = _read_int_3D(gpath+"get_rotation_perms_lists_HNF.in."+str(case))
        L = _read_int_3D(gpath+"get_rotation_perms_lists_L.in."+str(case))
        SNF = _read_int_3D(gpath+"get_rotation_perms_lists_SNF.in."+str(case))
        Op = _read_fixOp_1D(gpath+"get_rotation_perms_lists_Op.in."+str(case))
        dperms = _read_RotPermList(gpath+"get_rotation_perms_lists_dperms.in."+str(case))
        RPlist = _read_RotPermList_1D(gpath+"get_rotation_perms_lists_RPlist.in."+str(case))
        eps = _read_float(gpath+"get_rotation_perms_lists_eps.in."+str(case))
        out1 = _get_rotation_perms_lists(A,HNF,L,SNF,Op,RPlist,dperms,eps)
        out2 = _read_RotPermList_1D(gpath+"get_rotation_perms_lists_RPlist.out."+str(case))
        for out in out1:
            out.perm.arrow_perm = None
        self._compare_outputs(out1,out2)

    def test_56(self):
        from phenum.grouptheory import _get_rotation_perms_lists
        case = 56
        A = list(map(list,zip(*_read_float_2D(gpath+"get_rotation_perms_lists_A.in."+str(case)))))
        HNF = _read_int_3D(gpath+"get_rotation_perms_lists_HNF.in."+str(case))
        L = _read_int_3D(gpath+"get_rotation_perms_lists_L.in."+str(case))
        SNF = _read_int_3D(gpath+"get_rotation_perms_lists_SNF.in."+str(case))
        Op = _read_fixOp_1D(gpath+"get_rotation_perms_lists_Op.in."+str(case))
        dperms = _read_RotPermList(gpath+"get_rotation_perms_lists_dperms.in."+str(case))
        RPlist = _read_RotPermList_1D(gpath+"get_rotation_perms_lists_RPlist.in."+str(case))
        eps = _read_float(gpath+"get_rotation_perms_lists_eps.in."+str(case))
        out1 = _get_rotation_perms_lists(A,HNF,L,SNF,Op,RPlist,dperms,eps)
        out2 = _read_RotPermList_1D(gpath+"get_rotation_perms_lists_RPlist.out."+str(case))
        for out in out1:
            out.perm.arrow_perm = None
        self._compare_outputs(out1,out2)

    def test_57(self):
        from phenum.grouptheory import _get_rotation_perms_lists
        case = 57
        A = list(map(list,zip(*_read_float_2D(gpath+"get_rotation_perms_lists_A.in."+str(case)))))
        HNF = _read_int_3D(gpath+"get_rotation_perms_lists_HNF.in."+str(case))
        L = _read_int_3D(gpath+"get_rotation_perms_lists_L.in."+str(case))
        SNF = _read_int_3D(gpath+"get_rotation_perms_lists_SNF.in."+str(case))
        Op = _read_fixOp_1D(gpath+"get_rotation_perms_lists_Op.in."+str(case))
        dperms = _read_RotPermList(gpath+"get_rotation_perms_lists_dperms.in."+str(case))
        RPlist = _read_RotPermList_1D(gpath+"get_rotation_perms_lists_RPlist.in."+str(case))
        eps = _read_float(gpath+"get_rotation_perms_lists_eps.in."+str(case))
        out1 = _get_rotation_perms_lists(A,HNF,L,SNF,Op,RPlist,dperms,eps)
        out2 = _read_RotPermList_1D(gpath+"get_rotation_perms_lists_RPlist.out."+str(case))
        for out in out1:
            out.perm.arrow_perm = None
        self._compare_outputs(out1,out2)

    def test_58(self):
        from phenum.grouptheory import _get_rotation_perms_lists
        case = 58
        A = list(map(list,zip(*_read_float_2D(gpath+"get_rotation_perms_lists_A.in."+str(case)))))
        HNF = _read_int_3D(gpath+"get_rotation_perms_lists_HNF.in."+str(case))
        L = _read_int_3D(gpath+"get_rotation_perms_lists_L.in."+str(case))
        SNF = _read_int_3D(gpath+"get_rotation_perms_lists_SNF.in."+str(case))
        Op = _read_fixOp_1D(gpath+"get_rotation_perms_lists_Op.in."+str(case))
        dperms = _read_RotPermList(gpath+"get_rotation_perms_lists_dperms.in."+str(case))
        RPlist = _read_RotPermList_1D(gpath+"get_rotation_perms_lists_RPlist.in."+str(case))
        eps = _read_float(gpath+"get_rotation_perms_lists_eps.in."+str(case))
        out1 = _get_rotation_perms_lists(A,HNF,L,SNF,Op,RPlist,dperms,eps)
        out2 = _read_RotPermList_1D(gpath+"get_rotation_perms_lists_RPlist.out."+str(case))
        for out in out1:
            out.perm.arrow_perm = None
        self._compare_outputs(out1,out2)

class TestRM3DOperations(ut.TestCase):
    """Tests of the _rm_3D_operations subroutine."""

    def test_1(self):
        from phenum.grouptheory import _rm_3D_operations
        with pytest.raises(ValueError):
            _rm_3D_operations([[1,1,0],[1,1,1],[0,1,1]],[0],[0],1E-7)
            
    
class TestGetSymGroup(ut.TestCase):
    """ Tests of the get_sym_group subroutine."""
        
    def test_1(self):
        from phenum.grouptheory import get_sym_group
        par_lat = [[0.5, 0.5, 0.0], [0.5, 0.0, 0.5], [0.0, 0.5, 0.5]]
        bas_vacs = [[0.0, 0.0, 0.0]]
        HNF = [[1, 0, 0], [0, 1, 0], [2, 3, 6]]
        LatDim = 3
        out = _read_output("arrow_group.out.1")
        symm = get_sym_group(par_lat,bas_vacs,HNF,LatDim)
        agroup = []
        for i in range(len(symm.perm.site_perm)):
            agroup.append([symm.perm.site_perm[i],[int(j) for j in symm.perm.arrow_perm[i]]])
        for perm in agroup:
            self.assertIn(perm,out)

    def test_2(self):
        from phenum.grouptheory import get_sym_group
        par_lat = [[0.5, 0.5, 0.0], [0.5, 0.0, 0.5], [0.0, 0.5, 0.5]]
        bas_vacs = [[0.0, 0.0, 0.0]]
        HNF = [[1, 0, 0], [0, 1, 0], [0, 5, 6]]
        LatDim = 3
        out = _read_output("arrow_group.out.2")
        symm = get_sym_group(par_lat,bas_vacs,HNF,LatDim)
        agroup = []
        for i in range(len(symm.perm.site_perm)):
            agroup.append([symm.perm.site_perm[i],[int(j) for j in symm.perm.arrow_perm[i]]])
        for perm in agroup:
            self.assertIn(perm,out)

    def test_3(self):
        from phenum.grouptheory import get_sym_group
        par_lat = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
        bas_vacs = [[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]]
        HNF = [[1, 0, 0], [0, 1, 0], [0, 1, 2]]
        LatDim = 3
        out = _read_output("arrow_group.out.3")
        symm = get_sym_group(par_lat,bas_vacs,HNF,LatDim)
        agroup = []
        for i in range(len(symm.perm.site_perm)):
            agroup.append([symm.perm.site_perm[i],[int(j) for j in symm.perm.arrow_perm[i]]])
        for perm in agroup:
            self.assertIn(perm,out)

    def test_4(self):
        from phenum.grouptheory import get_sym_group
        par_lat = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
        bas_vacs = [[0.0, 0.0, 0.0]]
        HNF = [[1, 0, 0], [0, 1, 0], [0, 0, 7]]
        LatDim = 3
        out = _read_output("arrow_group.out.4")
        symm = get_sym_group(par_lat,bas_vacs,HNF,LatDim)
        agroup = []
        for i in range(len(symm.perm.site_perm)):
            agroup.append([symm.perm.site_perm[i],[int(j) for j in symm.perm.arrow_perm[i]]])
        for perm in agroup:
            self.assertIn(perm,out)

    def test_5(self):
        from phenum.grouptheory import get_sym_group
        par_lat = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
        bas_vacs = [[0.0, 0.0, 0.0]]
        HNF = [[1, 0, 0], [1, 2, 0], [1, 0, 2]]
        LatDim = 3
        out = _read_output("arrow_group.out.5")
        symm = get_sym_group(par_lat,bas_vacs,HNF,LatDim)
        agroup = []
        for i in range(len(symm.perm.site_perm)):
            agroup.append([symm.perm.site_perm[i],[int(j) for j in symm.perm.arrow_perm[i]]])
        for perm in agroup:
            self.assertIn(perm,out)

    def test_6(self):
        from phenum.grouptheory import get_sym_group
        par_lat = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
        bas_vacs = [[0.0, 0.0, 0.0]]
        HNF = [[1, 0, 0], [0, 2, 0], [0, 0, 2]]
        LatDim = 3
        out = _read_output("arrow_group.out.6")
        symm = get_sym_group(par_lat,bas_vacs,HNF,LatDim)
        agroup = []
        for i in range(len(symm.perm.site_perm)):
            agroup.append([symm.perm.site_perm[i],[int(j) for j in symm.perm.arrow_perm[i]]])
        for perm in agroup:
            self.assertIn(perm,out)

    def test_7(self):
        from phenum.grouptheory import get_sym_group
        par_lat = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
        bas_vacs = [[0.0, 0.0, 0.0], [0.5, 0.5, 0.5], [0.25, 0.25, 0.75]]
        HNF = [[1, 0, 0], [0, 1, 0], [0, 1, 2]]
        LatDim = 3
        out = _read_output("arrow_group.out.7")
        symm = get_sym_group(par_lat,bas_vacs,HNF,LatDim)
        agroup = []
        for i in range(len(symm.perm.site_perm)):
            agroup.append([symm.perm.site_perm[i],[int(j) for j in symm.perm.arrow_perm[i]]])
        for perm in agroup:
            self.assertIn(perm,out)

    def test_8(self):
        from phenum.grouptheory import get_sym_group
        par_lat = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
        bas_vacs = [[0.0, 0.0, 0.0], [0.5, 0.5, 0.5], [0.25, 0.25, 0.75]]
        HNF = [[1, 0, 0], [0, 1, 0], [0, 0, 3]]
        LatDim = 3
        out = _read_output("arrow_group.out.8")
        symm = get_sym_group(par_lat,bas_vacs,HNF,LatDim)
        agroup = []
        for i in range(len(symm.perm.site_perm)):
            agroup.append([symm.perm.site_perm[i],[int(j) for j in symm.perm.arrow_perm[i]]])
        for perm in agroup:
            self.assertIn(perm,out)

    def test_9(self):
        from phenum.grouptheory import get_sym_group
        par_lat = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
        bas_vacs = [[0.0, 0.0, 0.0], [0.5, 0.5, 0.5], [0.25, 0.25, 0.75]]
        HNF = [[1, 0, 0], [0, 1, 0], [0, 2, 3]]
        LatDim = 3
        out = _read_output("arrow_group.out.9")
        symm = get_sym_group(par_lat,bas_vacs,HNF,LatDim)
        agroup = []
        for i in range(len(symm.perm.site_perm)):
            agroup.append([symm.perm.site_perm[i],[int(j) for j in symm.perm.arrow_perm[i]]])
        for perm in agroup:
            self.assertIn(perm,out)

    def test_10(self):
        from phenum.grouptheory import get_sym_group
        par_lat = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
        bas_vacs = [[0.0, 0.0, 0.0]]
        HNF = [[2, 0, 0], [0, 2, 0], [0, 0, 2]]
        LatDim = 3
        out = _read_output("arrow_group.out.10")
        symm = get_sym_group(par_lat,bas_vacs,HNF,LatDim)
        agroup = []
        for i in range(len(symm.perm.site_perm)):
            agroup.append([symm.perm.site_perm[i],[int(j) for j in symm.perm.arrow_perm[i]]])
        for perm in agroup:
            self.assertIn(perm,out)

    def test_11(self):
        from phenum.grouptheory import get_sym_group
        par_lat = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
        bas_vacs = [[2.0, 2.0, 2.0]]
        HNF = [[2, 0, 0], [0, 2, 0], [0, 0, 2]]
        LatDim = 3
        out = _read_output("arrow_group.out.10")
        symm = get_sym_group(par_lat,bas_vacs,HNF,LatDim)
        agroup = []
        for i in range(len(symm.perm.site_perm)):
            agroup.append([symm.perm.site_perm[i],[int(j) for j in symm.perm.arrow_perm[i]]])
        for perm in agroup:
            self.assertIn(perm,out)
