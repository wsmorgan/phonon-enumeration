"""Methods for testing the subroutines in the symmetry module."""
import unittest as ut
import numpy as np

gpath =  "tests/symmetry/"

def _read_pg(fname):
    pg = []
    array = []
    with open(fname,"r") as f:
        for line in f:
            temp = line.strip().split()
            if len(temp)==3:
                array.append([float(i) for i in temp])
            else:
                pg.append(array)
                array = []
    return pg                

def _read_float_3D(fname):
    array = []
    parray = []
    lc = 0
    dc = 0
    with open(fname,"r") as f1:
        for line in f1:
            lc +=1
            if lc == 2:
                d1 = int(line.strip().split()[1])
                d2 = int(line.strip().split()[2])
                d3 = int(line.strip().split()[3])
            elif lc > 3:
                if "#" not in line:
                    dc +=1
                    parray.append([float(i) for i in line.strip().split()])

                    if dc == 3:
                        array.append(list(map(list,zip(*parray))))
                        parray = []
                        dc = 0
    array2 = []
    for i in range(d3):
        array2.append(list(map(list,zip(*[array[0][i],array[1][i],array[2][i]]))))
    return array2

def _read_float_2D(fname):
    array = []
    with open(fname,"r") as f1:
        for line in f1:
            if "#" not in line:
                array.append([float(i) for i in line.strip().split()])
    return array

def _read_float_1D(fname):
    array = []
    with open(fname,"r") as f1:
        for line in f1:
            if "#" not in line:
                array = [float(i) for i in line.strip().split()]
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

def _read_output(test):
    values = []
    with open("tests/symmetry/"+test) as f:
        for line in f:
            values.append(eval(line))
    return values

def _read_spaceGroup(case):
    sg_ops = _read_float_3D(gpath+"get_spaceGroup_sg_op.out."+str(case))
    sg_fracts = list(map(list,zip(*_read_float_2D(gpath+"get_spaceGroup_sg_fract.out."+str(case)))))
        
    return [sg_ops,sg_fracts]
    
class TestGetConcsForSize(ut.TestCase):
    """Tests of the get_concs_for_size subroutine."""

    def test_1(self):
        from phenum.symmetry import get_concs_for_size
        size =  4
        nspecies =  2
        nB =  3
        res_concs =  True
        concs =  [[1, 1, 4], [1, 1, 4]]
        out =  []
        self.assertEqual(get_concs_for_size(size,nspecies,res_concs,nB,concs),out)

    def test_2(self):
        from phenum.symmetry import get_concs_for_size        
        size =  15
        nspecies =  3
        nB =  2
        res_concs =  True
        concs =  [[2, 8, 15], [5, 5, 15], [4, 8, 15]]
        out =  [[4, 10, 16], [5, 10, 15], [6, 10, 14], [7, 10, 13], [8, 10, 12], [9, 10, 11], [10, 10, 10], [11, 10, 9], [12, 10, 8]]
        self.assertEqual(get_concs_for_size(size,nspecies,res_concs,nB,concs),out)

    def test_3(self):
        from phenum.symmetry import get_concs_for_size
        size =  18
        nspecies =  4
        nB =  2
        res_concs =  True
        concs =  [[5, 5, 18], [6, 13, 18], [4, 16, 18], [5, 8, 18]]
        out =  []
        self.assertEqual(get_concs_for_size(size,nspecies,res_concs,nB,concs),out)

    def test_4(self):
        from phenum.symmetry import get_concs_for_size
        size =  16
        nspecies =  2
        nB =  3
        res_concs =  False
        concs =  []
        out =  [[0, 48], [1, 47], [2, 46], [3, 45], [4, 44], [5, 43], [6, 42], [7, 41], [8, 40], [9, 39], [10, 38], [11, 37], [12, 36], [13, 35], [14, 34], [15, 33], [16, 32], [17, 31], [18, 30], [19, 29], [20, 28], [21, 27], [22, 26], [23, 25], [24, 24], [25, 23], [26, 22], [27, 21], [28, 20], [29, 19], [30, 18], [31, 17], [32, 16], [33, 15], [34, 14], [35, 13], [36, 12], [37, 11], [38, 10], [39, 9], [40, 8], [41, 7], [42, 6], [43, 5], [44, 4], [45, 3], [46, 2], [47, 1], [48, 0]]
        self.assertEqual(get_concs_for_size(size,nspecies,res_concs,nB,concs),out)

    def test_5(self):
        from phenum.symmetry import get_concs_for_size
        size =  13
        nspecies =  3
        nB =  2
        res_concs =  True
        concs =  [[3, 7, 13], [5, 12, 13], [1, 3, 13]]
        out =  [[6, 14, 6], [6, 15, 5], [6, 16, 4], [6, 17, 3], [6, 18, 2], [7, 13, 6], [7, 14, 5], [7, 15, 4], [7, 16, 3], [7, 17, 2], [8, 12, 6], [8, 13, 5], [8, 14, 4], [8, 15, 3], [8, 16, 2], [9, 11, 6], [9, 12, 5], [9, 13, 4], [9, 14, 3], [9, 15, 2], [10, 10, 6], [10, 11, 5], [10, 12, 4], [10, 13, 3], [10, 14, 2], [11, 10, 5], [11, 11, 4], [11, 12, 3], [11, 13, 2], [12, 10, 4], [12, 11, 3], [12, 12, 2], [13, 10, 3], [13, 11, 2], [14, 10, 2]]
        self.assertEqual(get_concs_for_size(size,nspecies,res_concs,nB,concs),out)

    def test_6(self):
        from phenum.symmetry import get_concs_for_size
        size =  13
        nspecies =  4
        nB =  3
        res_concs =  True
        concs =  [[2, 5, 13], [5, 12, 13], [3, 4, 13], [4, 6, 13]]
        out =  []
        self.assertEqual(get_concs_for_size(size,nspecies,res_concs,nB,concs),out)

    def test_7(self):
        from phenum.symmetry import get_concs_for_size
        size =  16
        nspecies =  3
        nB =  2
        res_concs =  False
        concs =  []
        out =  [[0, 0, 32], [0, 1, 31], [0, 2, 30], [0, 3, 29], [0, 4, 28], [0, 5, 27], [0, 6, 26], [0, 7, 25], [0, 8, 24], [0, 9, 23], [0, 10, 22], [0, 11, 21], [0, 12, 20], [0, 13, 19], [0, 14, 18], [0, 15, 17], [0, 16, 16], [0, 17, 15], [0, 18, 14], [0, 19, 13], [0, 20, 12], [0, 21, 11], [0, 22, 10], [0, 23, 9], [0, 24, 8], [0, 25, 7], [0, 26, 6], [0, 27, 5], [0, 28, 4], [0, 29, 3], [0, 30, 2], [0, 31, 1], [0, 32, 0], [1, 0, 31], [1, 1, 30], [1, 2, 29], [1, 3, 28], [1, 4, 27], [1, 5, 26], [1, 6, 25], [1, 7, 24], [1, 8, 23], [1, 9, 22], [1, 10, 21], [1, 11, 20], [1, 12, 19], [1, 13, 18], [1, 14, 17], [1, 15, 16], [1, 16, 15], [1, 17, 14], [1, 18, 13], [1, 19, 12], [1, 20, 11], [1, 21, 10], [1, 22, 9], [1, 23, 8], [1, 24, 7], [1, 25, 6], [1, 26, 5], [1, 27, 4], [1, 28, 3], [1, 29, 2], [1, 30, 1], [1, 31, 0], [2, 0, 30], [2, 1, 29], [2, 2, 28], [2, 3, 27], [2, 4, 26], [2, 5, 25], [2, 6, 24], [2, 7, 23], [2, 8, 22], [2, 9, 21], [2, 10, 20], [2, 11, 19], [2, 12, 18], [2, 13, 17], [2, 14, 16], [2, 15, 15], [2, 16, 14], [2, 17, 13], [2, 18, 12], [2, 19, 11], [2, 20, 10], [2, 21, 9], [2, 22, 8], [2, 23, 7], [2, 24, 6], [2, 25, 5], [2, 26, 4], [2, 27, 3], [2, 28, 2], [2, 29, 1], [2, 30, 0], [3, 0, 29], [3, 1, 28], [3, 2, 27], [3, 3, 26], [3, 4, 25], [3, 5, 24], [3, 6, 23], [3, 7, 22], [3, 8, 21], [3, 9, 20], [3, 10, 19], [3, 11, 18], [3, 12, 17], [3, 13, 16], [3, 14, 15], [3, 15, 14], [3, 16, 13], [3, 17, 12], [3, 18, 11], [3, 19, 10], [3, 20, 9], [3, 21, 8], [3, 22, 7], [3, 23, 6], [3, 24, 5], [3, 25, 4], [3, 26, 3], [3, 27, 2], [3, 28, 1], [3, 29, 0], [4, 0, 28], [4, 1, 27], [4, 2, 26], [4, 3, 25], [4, 4, 24], [4, 5, 23], [4, 6, 22], [4, 7, 21], [4, 8, 20], [4, 9, 19], [4, 10, 18], [4, 11, 17], [4, 12, 16], [4, 13, 15], [4, 14, 14], [4, 15, 13], [4, 16, 12], [4, 17, 11], [4, 18, 10], [4, 19, 9], [4, 20, 8], [4, 21, 7], [4, 22, 6], [4, 23, 5], [4, 24, 4], [4, 25, 3], [4, 26, 2], [4, 27, 1], [4, 28, 0], [5, 0, 27], [5, 1, 26], [5, 2, 25], [5, 3, 24], [5, 4, 23], [5, 5, 22], [5, 6, 21], [5, 7, 20], [5, 8, 19], [5, 9, 18], [5, 10, 17], [5, 11, 16], [5, 12, 15], [5, 13, 14], [5, 14, 13], [5, 15, 12], [5, 16, 11], [5, 17, 10], [5, 18, 9], [5, 19, 8], [5, 20, 7], [5, 21, 6], [5, 22, 5], [5, 23, 4], [5, 24, 3], [5, 25, 2], [5, 26, 1], [5, 27, 0], [6, 0, 26], [6, 1, 25], [6, 2, 24], [6, 3, 23], [6, 4, 22], [6, 5, 21], [6, 6, 20], [6, 7, 19], [6, 8, 18], [6, 9, 17], [6, 10, 16], [6, 11, 15], [6, 12, 14], [6, 13, 13], [6, 14, 12], [6, 15, 11], [6, 16, 10], [6, 17, 9], [6, 18, 8], [6, 19, 7], [6, 20, 6], [6, 21, 5], [6, 22, 4], [6, 23, 3], [6, 24, 2], [6, 25, 1], [6, 26, 0], [7, 0, 25], [7, 1, 24], [7, 2, 23], [7, 3, 22], [7, 4, 21], [7, 5, 20], [7, 6, 19], [7, 7, 18], [7, 8, 17], [7, 9, 16], [7, 10, 15], [7, 11, 14], [7, 12, 13], [7, 13, 12], [7, 14, 11], [7, 15, 10], [7, 16, 9], [7, 17, 8], [7, 18, 7], [7, 19, 6], [7, 20, 5], [7, 21, 4], [7, 22, 3], [7, 23, 2], [7, 24, 1], [7, 25, 0], [8, 0, 24], [8, 1, 23], [8, 2, 22], [8, 3, 21], [8, 4, 20], [8, 5, 19], [8, 6, 18], [8, 7, 17], [8, 8, 16], [8, 9, 15], [8, 10, 14], [8, 11, 13], [8, 12, 12], [8, 13, 11], [8, 14, 10], [8, 15, 9], [8, 16, 8], [8, 17, 7], [8, 18, 6], [8, 19, 5], [8, 20, 4], [8, 21, 3], [8, 22, 2], [8, 23, 1], [8, 24, 0], [9, 0, 23], [9, 1, 22], [9, 2, 21], [9, 3, 20], [9, 4, 19], [9, 5, 18], [9, 6, 17], [9, 7, 16], [9, 8, 15], [9, 9, 14], [9, 10, 13], [9, 11, 12], [9, 12, 11], [9, 13, 10], [9, 14, 9], [9, 15, 8], [9, 16, 7], [9, 17, 6], [9, 18, 5], [9, 19, 4], [9, 20, 3], [9, 21, 2], [9, 22, 1], [9, 23, 0], [10, 0, 22], [10, 1, 21], [10, 2, 20], [10, 3, 19], [10, 4, 18], [10, 5, 17], [10, 6, 16], [10, 7, 15], [10, 8, 14], [10, 9, 13], [10, 10, 12], [10, 11, 11], [10, 12, 10], [10, 13, 9], [10, 14, 8], [10, 15, 7], [10, 16, 6], [10, 17, 5], [10, 18, 4], [10, 19, 3], [10, 20, 2], [10, 21, 1], [10, 22, 0], [11, 0, 21], [11, 1, 20], [11, 2, 19], [11, 3, 18], [11, 4, 17], [11, 5, 16], [11, 6, 15], [11, 7, 14], [11, 8, 13], [11, 9, 12], [11, 10, 11], [11, 11, 10], [11, 12, 9], [11, 13, 8], [11, 14, 7], [11, 15, 6], [11, 16, 5], [11, 17, 4], [11, 18, 3], [11, 19, 2], [11, 20, 1], [11, 21, 0], [12, 0, 20], [12, 1, 19], [12, 2, 18], [12, 3, 17], [12, 4, 16], [12, 5, 15], [12, 6, 14], [12, 7, 13], [12, 8, 12], [12, 9, 11], [12, 10, 10], [12, 11, 9], [12, 12, 8], [12, 13, 7], [12, 14, 6], [12, 15, 5], [12, 16, 4], [12, 17, 3], [12, 18, 2], [12, 19, 1], [12, 20, 0], [13, 0, 19], [13, 1, 18], [13, 2, 17], [13, 3, 16], [13, 4, 15], [13, 5, 14], [13, 6, 13], [13, 7, 12], [13, 8, 11], [13, 9, 10], [13, 10, 9], [13, 11, 8], [13, 12, 7], [13, 13, 6], [13, 14, 5], [13, 15, 4], [13, 16, 3], [13, 17, 2], [13, 18, 1], [13, 19, 0], [14, 0, 18], [14, 1, 17], [14, 2, 16], [14, 3, 15], [14, 4, 14], [14, 5, 13], [14, 6, 12], [14, 7, 11], [14, 8, 10], [14, 9, 9], [14, 10, 8], [14, 11, 7], [14, 12, 6], [14, 13, 5], [14, 14, 4], [14, 15, 3], [14, 16, 2], [14, 17, 1], [14, 18, 0], [15, 0, 17], [15, 1, 16], [15, 2, 15], [15, 3, 14], [15, 4, 13], [15, 5, 12], [15, 6, 11], [15, 7, 10], [15, 8, 9], [15, 9, 8], [15, 10, 7], [15, 11, 6], [15, 12, 5], [15, 13, 4], [15, 14, 3], [15, 15, 2], [15, 16, 1], [15, 17, 0], [16, 0, 16], [16, 1, 15], [16, 2, 14], [16, 3, 13], [16, 4, 12], [16, 5, 11], [16, 6, 10], [16, 7, 9], [16, 8, 8], [16, 9, 7], [16, 10, 6], [16, 11, 5], [16, 12, 4], [16, 13, 3], [16, 14, 2], [16, 15, 1], [16, 16, 0], [17, 0, 15], [17, 1, 14], [17, 2, 13], [17, 3, 12], [17, 4, 11], [17, 5, 10], [17, 6, 9], [17, 7, 8], [17, 8, 7], [17, 9, 6], [17, 10, 5], [17, 11, 4], [17, 12, 3], [17, 13, 2], [17, 14, 1], [17, 15, 0], [18, 0, 14], [18, 1, 13], [18, 2, 12], [18, 3, 11], [18, 4, 10], [18, 5, 9], [18, 6, 8], [18, 7, 7], [18, 8, 6], [18, 9, 5], [18, 10, 4], [18, 11, 3], [18, 12, 2], [18, 13, 1], [18, 14, 0], [19, 0, 13], [19, 1, 12], [19, 2, 11], [19, 3, 10], [19, 4, 9], [19, 5, 8], [19, 6, 7], [19, 7, 6], [19, 8, 5], [19, 9, 4], [19, 10, 3], [19, 11, 2], [19, 12, 1], [19, 13, 0], [20, 0, 12], [20, 1, 11], [20, 2, 10], [20, 3, 9], [20, 4, 8], [20, 5, 7], [20, 6, 6], [20, 7, 5], [20, 8, 4], [20, 9, 3], [20, 10, 2], [20, 11, 1], [20, 12, 0], [21, 0, 11], [21, 1, 10], [21, 2, 9], [21, 3, 8], [21, 4, 7], [21, 5, 6], [21, 6, 5], [21, 7, 4], [21, 8, 3], [21, 9, 2], [21, 10, 1], [21, 11, 0], [22, 0, 10], [22, 1, 9], [22, 2, 8], [22, 3, 7], [22, 4, 6], [22, 5, 5], [22, 6, 4], [22, 7, 3], [22, 8, 2], [22, 9, 1], [22, 10, 0], [23, 0, 9], [23, 1, 8], [23, 2, 7], [23, 3, 6], [23, 4, 5], [23, 5, 4], [23, 6, 3], [23, 7, 2], [23, 8, 1], [23, 9, 0], [24, 0, 8], [24, 1, 7], [24, 2, 6], [24, 3, 5], [24, 4, 4], [24, 5, 3], [24, 6, 2], [24, 7, 1], [24, 8, 0], [25, 0, 7], [25, 1, 6], [25, 2, 5], [25, 3, 4], [25, 4, 3], [25, 5, 2], [25, 6, 1], [25, 7, 0], [26, 0, 6], [26, 1, 5], [26, 2, 4], [26, 3, 3], [26, 4, 2], [26, 5, 1], [26, 6, 0], [27, 0, 5], [27, 1, 4], [27, 2, 3], [27, 3, 2], [27, 4, 1], [27, 5, 0], [28, 0, 4], [28, 1, 3], [28, 2, 2], [28, 3, 1], [28, 4, 0], [29, 0, 3], [29, 1, 2], [29, 2, 1], [29, 3, 0], [30, 0, 2], [30, 1, 1], [30, 2, 0], [31, 0, 1], [31, 1, 0], [32, 0, 0]]
        self.assertEqual(get_concs_for_size(size,nspecies,res_concs,nB,concs),out)

    def test_8(self):
        from phenum.symmetry import get_concs_for_size
        size =  18
        nspecies =  3
        nB =  4
        res_concs =  True
        concs =  [[3, 7, 18], [1, 8, 18], [6, 7, 18]]
        out =  [[12, 32, 28], [13, 31, 28], [13, 32, 27], [14, 30, 28], [14, 31, 27], [14, 32, 26], [15, 29, 28], [15, 30, 27], [15, 31, 26], [15, 32, 25], [16, 28, 28], [16, 29, 27], [16, 30, 26], [16, 31, 25], [16, 32, 24], [17, 27, 28], [17, 28, 27], [17, 29, 26], [17, 30, 25], [17, 31, 24], [18, 26, 28], [18, 27, 27], [18, 28, 26], [18, 29, 25], [18, 30, 24], [19, 25, 28], [19, 26, 27], [19, 27, 26], [19, 28, 25], [19, 29, 24], [20, 24, 28], [20, 25, 27], [20, 26, 26], [20, 27, 25], [20, 28, 24], [21, 23, 28], [21, 24, 27], [21, 25, 26], [21, 26, 25], [21, 27, 24], [22, 22, 28], [22, 23, 27], [22, 24, 26], [22, 25, 25], [22, 26, 24], [23, 21, 28], [23, 22, 27], [23, 23, 26], [23, 24, 25], [23, 25, 24], [24, 20, 28], [24, 21, 27], [24, 22, 26], [24, 23, 25], [24, 24, 24], [25, 19, 28], [25, 20, 27], [25, 21, 26], [25, 22, 25], [25, 23, 24], [26, 18, 28], [26, 19, 27], [26, 20, 26], [26, 21, 25], [26, 22, 24], [27, 17, 28], [27, 18, 27], [27, 19, 26], [27, 20, 25], [27, 21, 24], [28, 16, 28], [28, 17, 27], [28, 18, 26], [28, 19, 25], [28, 20, 24]]
        self.assertEqual(get_concs_for_size(size,nspecies,res_concs,nB,concs),out)

    def test_9(self):
        from phenum.symmetry import get_concs_for_size
        size =  15
        nspecies =  2
        nB =  3
        res_concs =  False
        concs =  []
        out =  [[0, 45], [1, 44], [2, 43], [3, 42], [4, 41], [5, 40], [6, 39], [7, 38], [8, 37], [9, 36], [10, 35], [11, 34], [12, 33], [13, 32], [14, 31], [15, 30], [16, 29], [17, 28], [18, 27], [19, 26], [20, 25], [21, 24], [22, 23], [23, 22], [24, 21], [25, 20], [26, 19], [27, 18], [28, 17], [29, 16], [30, 15], [31, 14], [32, 13], [33, 12], [34, 11], [35, 10], [36, 9], [37, 8], [38, 7], [39, 6], [40, 5], [41, 4], [42, 3], [43, 2], [44, 1], [45, 0]]
        self.assertEqual(get_concs_for_size(size,nspecies,res_concs,nB,concs),out)

    def test_10(self):
        from phenum.symmetry import get_concs_for_size
        size =  14
        nspecies =  2
        nB =  4
        res_concs =  True
        concs =  [[3, 3, 14], [2, 13, 14]]
        out =  [[12, 44]]
        self.assertEqual(get_concs_for_size(size,nspecies,res_concs,nB,concs),out)

class TestGetSpaceGroup(ut.TestCase):
    """Tests of the get_spaceGroup subroutine."""

    def _compare_space_group(self,out1,out2):
        ops1 = out1[0]
        ops2 = out2[0]
        fract1 = out1[1]
        fract2 = out2[1]

        if len(ops1) == len(ops2):
            for i in range(len(ops1)):
                for j in range(3):
                    for k in range(3):
                        self.assertAlmostEqual(ops1[i][j][k],ops2[i][j][k],places=12)
        else:
            self.assertEqual(len(ops1),len(ops2))

        if len(fract1) == len(fract2):
            for i in range(len(ops1)):
                for j in range(3):
                    self.assertAlmostEqual(fract1[i][j],fract2[i][j],places=12)
        else:
            self.assertEqual(len(fract1),len(fract2))

    def test_1(self):
        from phenum.symmetry import get_spaceGroup
        par_lat =  [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
        atomType =  [1]
        bas_vecs =  [[0.0, 0.0, 0.0]]
        eps =  1e-10
        lattcoords =  False
        out =  ([[[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, -1.0]], [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], [[-1.0, 0.0, 0.0], [0.0, 0.0, -1.0], [0.0, -1.0, 0.0]], [[-1.0, 0.0, 0.0], [0.0, 0.0, -1.0], [0.0, 1.0, 0.0]], [[-1.0, 0.0, 0.0], [0.0, 0.0, 1.0], [0.0, -1.0, 0.0]], [[-1.0, 0.0, 0.0], [0.0, 0.0, 1.0], [0.0, 1.0, 0.0]], [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, -1.0]], [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], [[0.0, -1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, -1.0]], [[0.0, -1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], [[0.0, -1.0, 0.0], [0.0, 0.0, -1.0], [-1.0, 0.0, 0.0]], [[0.0, -1.0, 0.0], [0.0, 0.0, -1.0], [1.0, 0.0, 0.0]], [[0.0, -1.0, 0.0], [0.0, 0.0, 1.0], [-1.0, 0.0, 0.0]], [[0.0, -1.0, 0.0], [0.0, 0.0, 1.0], [1.0, 0.0, 0.0]], [[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, -1.0]], [[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], [[0.0, 0.0, -1.0], [-1.0, 0.0, 0.0], [0.0, -1.0, 0.0]], [[0.0, 0.0, -1.0], [-1.0, 0.0, 0.0], [0.0, 1.0, 0.0]], [[0.0, 0.0, -1.0], [0.0, -1.0, 0.0], [-1.0, 0.0, 0.0]], [[0.0, 0.0, -1.0], [0.0, -1.0, 0.0], [1.0, 0.0, 0.0]], [[0.0, 0.0, -1.0], [0.0, 1.0, 0.0], [-1.0, 0.0, 0.0]], [[0.0, 0.0, -1.0], [0.0, 1.0, 0.0], [1.0, 0.0, 0.0]], [[0.0, 0.0, -1.0], [1.0, 0.0, 0.0], [0.0, -1.0, 0.0]], [[0.0, 0.0, -1.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]], [[0.0, 0.0, 1.0], [-1.0, 0.0, 0.0], [0.0, -1.0, 0.0]], [[0.0, 0.0, 1.0], [-1.0, 0.0, 0.0], [0.0, 1.0, 0.0]], [[0.0, 0.0, 1.0], [0.0, -1.0, 0.0], [-1.0, 0.0, 0.0]], [[0.0, 0.0, 1.0], [0.0, -1.0, 0.0], [1.0, 0.0, 0.0]], [[0.0, 0.0, 1.0], [0.0, 1.0, 0.0], [-1.0, 0.0, 0.0]], [[0.0, 0.0, 1.0], [0.0, 1.0, 0.0], [1.0, 0.0, 0.0]], [[0.0, 0.0, 1.0], [1.0, 0.0, 0.0], [0.0, -1.0, 0.0]], [[0.0, 0.0, 1.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]], [[0.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, -1.0]], [[0.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], [[0.0, 1.0, 0.0], [0.0, 0.0, -1.0], [-1.0, 0.0, 0.0]], [[0.0, 1.0, 0.0], [0.0, 0.0, -1.0], [1.0, 0.0, 0.0]], [[0.0, 1.0, 0.0], [0.0, 0.0, 1.0], [-1.0, 0.0, 0.0]], [[0.0, 1.0, 0.0], [0.0, 0.0, 1.0], [1.0, 0.0, 0.0]], [[0.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, -1.0]], [[0.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], [[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, -1.0]], [[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], [[1.0, 0.0, 0.0], [0.0, 0.0, -1.0], [0.0, -1.0, 0.0]], [[1.0, 0.0, 0.0], [0.0, 0.0, -1.0], [0.0, 1.0, 0.0]], [[1.0, 0.0, 0.0], [0.0, 0.0, 1.0], [0.0, -1.0, 0.0]], [[1.0, 0.0, 0.0], [0.0, 0.0, 1.0], [0.0, 1.0, 0.0]], [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, -1.0]], [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]], [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]])
        self.assertEqual(get_spaceGroup(par_lat,atomType,bas_vecs,eps,lattcoords),out)

    def test_5(self):
        from phenum.symmetry import get_spaceGroup
        par_lat =  [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
        atomType =  [1, 1]
        bas_vecs =  [[0.0, 0.0, 0.0], [0.25, 0.25, 0.75]]
        eps =  1e-10
        lattcoords =  False
        out =  ([[[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, -1.0]], [[-1.0, 0.0, 0.0], [0.0, 0.0, 1.0], [0.0, 1.0, 0.0]], [[0.0, -1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, -1.0]], [[0.0, -1.0, 0.0], [0.0, 0.0, 1.0], [1.0, 0.0, 0.0]], [[0.0, 0.0, -1.0], [0.0, 1.0, 0.0], [-1.0, 0.0, 0.0]], [[0.0, 0.0, -1.0], [1.0, 0.0, 0.0], [0.0, -1.0, 0.0]], [[0.0, 0.0, 1.0], [-1.0, 0.0, 0.0], [0.0, 1.0, 0.0]], [[0.0, 0.0, 1.0], [0.0, -1.0, 0.0], [1.0, 0.0, 0.0]], [[0.0, 1.0, 0.0], [0.0, 0.0, -1.0], [-1.0, 0.0, 0.0]], [[0.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], [[1.0, 0.0, 0.0], [0.0, 0.0, -1.0], [0.0, -1.0, 0.0]], [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]], [[0.25, 0.25, 0.75], [0.25, 0.25, 0.75], [0.25, 0.25, 0.75], [0.25, 0.25, 0.75], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.25, 0.25, 0.75], [0.25, 0.25, 0.75], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]])
        self.assertEqual(get_spaceGroup(par_lat,atomType,bas_vecs,eps,lattcoords),out)

    def test_7(self):
        from phenum.symmetry import get_spaceGroup
        par_lat =  [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
        atomType =  [1, 1, 1]
        bas_vecs =  [[0.0, 0.0, 0.0], [0.25, 0.25, 0.75], [0.5, 0.5, 0.25]]
        eps =  1e-10
        lattcoords =  False
        out =  ([[[0.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]], [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]])
        self.assertEqual(get_spaceGroup(par_lat,atomType,bas_vecs,eps,lattcoords),out)

    def test_9(self):
        from phenum.symmetry import get_spaceGroup
        par_lat =  [[0.0, 0.5, 0.5], [0.5, 0.5, 0.0], [0.5, 0.0, 0.5]]
        atomType =  [1]
        bas_vecs =  [[0.0, 0.0, 0.0]]
        eps =  1e-10
        lattcoords =  False
        out =  ([[[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, -1.0]], [[1.0, 0.0, 0.0], [0.0, 0.0, -1.0], [0.0, -1.0, 0.0]], [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, -1.0]], [[-1.0, 0.0, 0.0], [0.0, 0.0, -1.0], [0.0, -1.0, 0.0]], [[0.0, 0.0, -1.0], [0.0, -1.0, 0.0], [1.0, 0.0, 0.0]], [[0.0, 0.0, -1.0], [1.0, 0.0, 0.0], [0.0, -1.0, 0.0]], [[0.0, 0.0, 1.0], [0.0, -1.0, 0.0], [1.0, 0.0, 0.0]], [[0.0, 0.0, 1.0], [1.0, 0.0, 0.0], [0.0, -1.0, 0.0]], [[0.0, -1.0, 0.0], [0.0, 0.0, -1.0], [1.0, 0.0, 0.0]], [[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, -1.0]], [[0.0, 1.0, 0.0], [0.0, 0.0, -1.0], [1.0, 0.0, 0.0]], [[0.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, -1.0]], [[0.0, 0.0, -1.0], [0.0, -1.0, 0.0], [-1.0, 0.0, 0.0]], [[0.0, 0.0, 1.0], [0.0, -1.0, 0.0], [-1.0, 0.0, 0.0]], [[0.0, 0.0, -1.0], [-1.0, 0.0, 0.0], [0.0, -1.0, 0.0]], [[0.0, 0.0, 1.0], [-1.0, 0.0, 0.0], [0.0, -1.0, 0.0]], [[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], [[1.0, 0.0, 0.0], [0.0, 0.0, 1.0], [0.0, -1.0, 0.0]], [[-1.0, 0.0, 0.0], [0.0, 0.0, 1.0], [0.0, -1.0, 0.0]], [[0.0, -1.0, 0.0], [0.0, 0.0, -1.0], [-1.0, 0.0, 0.0]], [[0.0, -1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, -1.0]], [[0.0, 1.0, 0.0], [0.0, 0.0, -1.0], [-1.0, 0.0, 0.0]], [[0.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, -1.0]], [[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], [[0.0, -1.0, 0.0], [0.0, 0.0, 1.0], [1.0, 0.0, 0.0]], [[0.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], [[0.0, 1.0, 0.0], [0.0, 0.0, 1.0], [1.0, 0.0, 0.0]], [[1.0, 0.0, 0.0], [0.0, 0.0, -1.0], [0.0, 1.0, 0.0]], [[-1.0, 0.0, 0.0], [0.0, 0.0, -1.0], [0.0, 1.0, 0.0]], [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, -1.0]], [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, -1.0]], [[0.0, 0.0, -1.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]], [[0.0, 0.0, 1.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]], [[0.0, 0.0, -1.0], [0.0, 1.0, 0.0], [1.0, 0.0, 0.0]], [[0.0, 0.0, 1.0], [0.0, 1.0, 0.0], [1.0, 0.0, 0.0]], [[0.0, -1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], [[0.0, -1.0, 0.0], [0.0, 0.0, 1.0], [-1.0, 0.0, 0.0]], [[0.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], [[0.0, 1.0, 0.0], [0.0, 0.0, 1.0], [-1.0, 0.0, 0.0]], [[0.0, 0.0, -1.0], [-1.0, 0.0, 0.0], [0.0, 1.0, 0.0]], [[0.0, 0.0, -1.0], [0.0, 1.0, 0.0], [-1.0, 0.0, 0.0]], [[0.0, 0.0, 1.0], [-1.0, 0.0, 0.0], [0.0, 1.0, 0.0]], [[0.0, 0.0, 1.0], [0.0, 1.0, 0.0], [-1.0, 0.0, 0.0]], [[1.0, 0.0, 0.0], [0.0, 0.0, 1.0], [0.0, 1.0, 0.0]], [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], [[-1.0, 0.0, 0.0], [0.0, 0.0, 1.0], [0.0, 1.0, 0.0]], [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]], [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]])
        self.assertEqual(get_spaceGroup(par_lat,atomType,bas_vecs,eps,lattcoords),out)

    def test_getsg11(self):
        from phenum.symmetry import get_spaceGroup
        case = 1
        par_lat =  list(map(list,zip(*_read_float_2D(gpath+"get_spaceGroup_aVecs.in."+str(case)))))
        atomType =  _read_int_1D(gpath+"get_spaceGroup_atomType.in."+str(case))
        bas_vecs =  list(map(list,zip(*_read_float_2D(gpath+"get_spaceGroup_input_pos.in."+str(case)))))
        eps =  _read_float(gpath+"get_spaceGroup_eps.in."+str(case))
        lattcoords =  _read_logical(gpath+"get_spaceGroup_lattcoords.in."+str(case))

        out = _read_spaceGroup(case)
        ops, fract = get_spaceGroup(par_lat,atomType,bas_vecs,eps,lattcoords)
        self._compare_space_group(out,[ops,fract])

    def test_getsg12(self):
        from phenum.symmetry import get_spaceGroup
        case = 2
        par_lat =  list(map(list,zip(*_read_float_2D(gpath+"get_spaceGroup_aVecs.in."+str(case)))))
        atomType =  _read_int_1D(gpath+"get_spaceGroup_atomType.in."+str(case))
        bas_vecs =  list(map(list,zip(*_read_float_2D(gpath+"get_spaceGroup_input_pos.in."+str(case)))))
        eps =  _read_float(gpath+"get_spaceGroup_eps.in."+str(case))
        lattcoords =  _read_logical(gpath+"get_spaceGroup_lattcoords.in."+str(case))

        out = _read_spaceGroup(case)
        ops, fract = get_spaceGroup(par_lat,atomType,bas_vecs,eps,lattcoords)
        self._compare_space_group(out,[ops,fract])

    def test_getsg13(self):
        from phenum.symmetry import get_spaceGroup
        case = 3
        par_lat =  list(map(list,zip(*_read_float_2D(gpath+"get_spaceGroup_aVecs.in."+str(case)))))
        atomType =  _read_int_1D(gpath+"get_spaceGroup_atomType.in."+str(case))
        bas_vecs =  list(map(list,zip(*_read_float_2D(gpath+"get_spaceGroup_input_pos.in."+str(case)))))
        eps =  _read_float(gpath+"get_spaceGroup_eps.in."+str(case))
        lattcoords =  _read_logical(gpath+"get_spaceGroup_lattcoords.in."+str(case))

        out = _read_spaceGroup(case)
        ops, fract = get_spaceGroup(par_lat,atomType,bas_vecs,eps,lattcoords)
        self._compare_space_group(out,[ops,fract])

    def test_getsg14(self):
        from phenum.symmetry import get_spaceGroup
        case = 4
        par_lat =  list(map(list,zip(*_read_float_2D(gpath+"get_spaceGroup_aVecs.in."+str(case)))))
        atomType =  _read_int_1D(gpath+"get_spaceGroup_atomType.in."+str(case))
        bas_vecs =  list(map(list,zip(*_read_float_2D(gpath+"get_spaceGroup_input_pos.in."+str(case)))))
        eps =  _read_float(gpath+"get_spaceGroup_eps.in."+str(case))
        lattcoords =  _read_logical(gpath+"get_spaceGroup_lattcoords.in."+str(case))

        out = _read_spaceGroup(case)
        ops, fract = get_spaceGroup(par_lat,atomType,bas_vecs,eps,lattcoords)
        self._compare_space_group(out,[ops,fract])

    def test_getsg15(self):
        from phenum.symmetry import get_spaceGroup
        case = 5
        par_lat =  list(map(list,zip(*_read_float_2D(gpath+"get_spaceGroup_aVecs.in."+str(case)))))
        atomType =  _read_int_1D(gpath+"get_spaceGroup_atomType.in."+str(case))
        bas_vecs =  list(map(list,zip(*_read_float_2D(gpath+"get_spaceGroup_input_pos.in."+str(case)))))
        eps =  _read_float(gpath+"get_spaceGroup_eps.in."+str(case))
        lattcoords =  _read_logical(gpath+"get_spaceGroup_lattcoords.in."+str(case))

        out = _read_spaceGroup(case)
        ops, fract = get_spaceGroup(par_lat,atomType,bas_vecs,eps,lattcoords)
        self._compare_space_group(out,[ops,fract])

    def test_getsg16(self):
        from phenum.symmetry import get_spaceGroup
        case = 6
        par_lat =  list(map(list,zip(*_read_float_2D(gpath+"get_spaceGroup_aVecs.in."+str(case)))))
        atomType =  _read_int_1D(gpath+"get_spaceGroup_atomType.in."+str(case))
        bas_vecs =  list(map(list,zip(*_read_float_2D(gpath+"get_spaceGroup_input_pos.in."+str(case)))))
        eps =  _read_float(gpath+"get_spaceGroup_eps.in."+str(case))
        lattcoords =  _read_logical(gpath+"get_spaceGroup_lattcoords.in."+str(case))

        out = _read_spaceGroup(case)
        ops, fract = get_spaceGroup(par_lat,atomType,bas_vecs,eps,lattcoords)
        self._compare_space_group(out,[ops,fract])

    def test_getsg17(self):
        from phenum.symmetry import get_spaceGroup
        case = 7
        par_lat =  list(map(list,zip(*_read_float_2D(gpath+"get_spaceGroup_aVecs.in."+str(case)))))
        atomType =  _read_int_1D(gpath+"get_spaceGroup_atomType.in."+str(case))
        bas_vecs =  list(map(list,zip(*_read_float_2D(gpath+"get_spaceGroup_input_pos.in."+str(case)))))
        eps =  _read_float(gpath+"get_spaceGroup_eps.in."+str(case))
        lattcoords =  _read_logical(gpath+"get_spaceGroup_lattcoords.in."+str(case))

        out = _read_spaceGroup(case)

        ops, fract = get_spaceGroup(par_lat,atomType,bas_vecs,eps=eps,lattcoords=lattcoords)

        self._compare_space_group(out,[ops,fract])

    def test_getsg18(self):
        from phenum.symmetry import get_spaceGroup
        case = 8
        par_lat =  list(map(list,zip(*_read_float_2D(gpath+"get_spaceGroup_aVecs.in."+str(case)))))
        atomType =  _read_int_1D(gpath+"get_spaceGroup_atomType.in."+str(case))
        bas_vecs =  list(map(list,zip(*_read_float_2D(gpath+"get_spaceGroup_input_pos.in."+str(case)))))
        eps =  _read_float(gpath+"get_spaceGroup_eps.in."+str(case))
        lattcoords =  _read_logical(gpath+"get_spaceGroup_lattcoords.in."+str(case))

        out = _read_spaceGroup(case)
        ops, fract = get_spaceGroup(par_lat,atomType,bas_vecs,eps,lattcoords)
        self._compare_space_group(out,[ops,fract])

    def test_getsg19(self):
        from phenum.symmetry import get_spaceGroup
        case = 9
        par_lat =  list(map(list,zip(*_read_float_2D(gpath+"get_spaceGroup_aVecs.in."+str(case)))))
        atomType =  _read_int_1D(gpath+"get_spaceGroup_atomType.in."+str(case))
        bas_vecs =  list(map(list,zip(*_read_float_2D(gpath+"get_spaceGroup_input_pos.in."+str(case)))))
        eps =  _read_float(gpath+"get_spaceGroup_eps.in."+str(case))
        lattcoords =  _read_logical(gpath+"get_spaceGroup_lattcoords.in."+str(case))

        out = _read_spaceGroup(case)
        ops, fract = get_spaceGroup(par_lat,atomType,bas_vecs,eps,lattcoords)
        self._compare_space_group(out,[ops,fract])
            
    def test_getsg20(self):
        from phenum.symmetry import get_spaceGroup
        from numpy.testing import assert_allclose
        from numpy import array 
        case = 10
        par_lat =  list(map(list,zip(*_read_float_2D(gpath+"get_spaceGroup_aVecs.in."+str(case)))))
        atomType =  _read_int_1D(gpath+"get_spaceGroup_atomType.in."+str(case))
        bas_vecs =  list(map(list,zip(*_read_float_2D(gpath+"get_spaceGroup_input_pos.in."+str(case)))))
        eps =  _read_float(gpath+"get_spaceGroup_eps.in."+str(case))
        lattcoords =  _read_logical(gpath+"get_spaceGroup_lattcoords.in."+str(case))

        out = _read_spaceGroup(case)
        ops, fract = get_spaceGroup(par_lat,atomType,bas_vecs,eps,lattcoords)
        self._compare_space_group(out,[ops,fract])
        
    def test_getsg21(self):
        from phenum.symmetry import get_spaceGroup, _get_transformations
        par_lat =  [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
        (prim_to_cart, cart_to_prim) = _get_transformations(par_lat)
        atomType =  [1, 1, 1]
        bas_vecs =  [[0.0, 0.0, 0.0], [0.25, 0.25, 0.75], [0.5, 0.5, 0.25]]
        bas_vecs = [np.matmul(cart_to_prim, i).tolist() for i in bas_vecs]
        eps =  1e-10
        lattcoords =  True
        out =  ([[[0.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]], [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]])
        self.assertEqual(get_spaceGroup(par_lat,atomType,bas_vecs,eps,lattcoords),out)

class TestBringIntoCell(ut.TestCase):
    """Tests of the bring_into_cell subroutine."""

    def test_1(self):
        from phenum.symmetry import bring_into_cell
        bas_vecs =  [0.0, 0.0, 0.0]
        cart_to_latt =  [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
        latt_to_cart =  [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
        eps =  1e-10
        out =  [0.0, 0.0, 0.0]
        self.assertEqual(bring_into_cell(bas_vecs,cart_to_latt,latt_to_cart,eps),out)

    def test_2(self):
        from phenum.symmetry import bring_into_cell
        bas_vecs =  [0.0, 0.0, 0.0]
        cart_to_latt =  [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
        latt_to_cart =  [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
        eps =  1e-10
        out =  [0.0, 0.0, 0.0]
        self.assertEqual(bring_into_cell(bas_vecs,cart_to_latt,latt_to_cart,eps),out)

    def test_3(self):
        from phenum.symmetry import bring_into_cell
        bas_vecs =  [0.0, 0.0, 0.0]
        cart_to_latt =  [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
        latt_to_cart =  [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
        eps =  1e-10
        out =  [0.0, 0.0, 0.0]
        self.assertEqual(bring_into_cell(bas_vecs,cart_to_latt,latt_to_cart,eps),out)        

    def test_4(self):
        from phenum.symmetry import bring_into_cell
        bas_vecs =  [0.0, 0.0, 0.0]
        cart_to_latt =  [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
        latt_to_cart =  [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
        eps =  1e-10
        out =  [0.0, 0.0, 0.0]
        self.assertEqual(bring_into_cell(bas_vecs,cart_to_latt,latt_to_cart,eps),out)

    def test_5(self):
        from phenum.symmetry import bring_into_cell
        bas_vecs =  [0.0, 0.0, 0.0]
        cart_to_latt =  [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
        latt_to_cart =  [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
        eps =  1e-10
        out =  [0.0, 0.0, 0.0]
        self.assertEqual(bring_into_cell(bas_vecs,cart_to_latt,latt_to_cart,eps),out)

    def test_6(self):
        from phenum.symmetry import bring_into_cell
        bas_vecs =  [0.25, 0.25, 0.75]
        cart_to_latt =  [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
        latt_to_cart =  [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
        eps =  1e-10
        out =  [0.25, 0.25, 0.75]
        self.assertEqual(bring_into_cell(bas_vecs,cart_to_latt,latt_to_cart,eps),out)

    def test_7(self):
        from phenum.symmetry import bring_into_cell
        bas_vecs =  [0.0, 0.0, 0.0]
        cart_to_latt =  [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
        latt_to_cart =  [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
        eps =  1e-10
        out =  [0.0, 0.0, 0.0]
        self.assertEqual(bring_into_cell(bas_vecs,cart_to_latt,latt_to_cart,eps),out)

    def test_8(self):
        from phenum.symmetry import bring_into_cell
        bas_vecs =  [0.25, 0.25, 0.75]
        cart_to_latt =  [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
        latt_to_cart =  [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
        eps =  1e-10
        out =  [0.25, 0.25, 0.75]
        self.assertEqual(bring_into_cell(bas_vecs,cart_to_latt,latt_to_cart,eps),out)

    def test_9(self):
        from phenum.symmetry import bring_into_cell
        bas_vecs =  [0.0, 0.0, 0.0]
        cart_to_latt =  [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
        latt_to_cart =  [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
        eps =  1e-10
        out =  [0.0, 0.0, 0.0]
        self.assertEqual(bring_into_cell(bas_vecs,cart_to_latt,latt_to_cart,eps),out)

    def test_10(self):
        from phenum.symmetry import bring_into_cell
        bas_vecs =  [0.25, 0.25, 0.75]
        cart_to_latt =  [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
        latt_to_cart =  [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
        eps =  1e-10
        out =  [0.25, 0.25, 0.75]
        self.assertEqual(bring_into_cell(bas_vecs,cart_to_latt,latt_to_cart,eps),out)

    def test_11(self):
        from phenum.symmetry import bring_into_cell
        bas_vecs =  [0.5, 0.5, 0.25]
        cart_to_latt =  [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
        latt_to_cart =  [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
        eps =  1e-10
        out =  [0.5, 0.5, 0.25]
        self.assertEqual(bring_into_cell(bas_vecs,cart_to_latt,latt_to_cart,eps),out)

    def test_11(self):
        from phenum.symmetry import bring_into_cell
        bas_vecs =  [0.0, 0.0, 0.0]
        cart_to_latt =  [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
        latt_to_cart =  [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
        eps =  1e-10
        out =  [0.0, 0.0, 0.0]
        self.assertEqual(bring_into_cell(bas_vecs,cart_to_latt,latt_to_cart,eps),out)

    def test_12(self):
        from phenum.symmetry import bring_into_cell
        bas_vecs =  [0.25, 0.25, 0.75]
        cart_to_latt =  [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
        latt_to_cart =  [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
        eps =  1e-10
        out =  [0.25, 0.25, 0.75]
        self.assertEqual(bring_into_cell(bas_vecs,cart_to_latt,latt_to_cart,eps),out)

    def test_13(self):
        from phenum.symmetry import bring_into_cell
        bas_vecs =  [0.5, 0.5, 0.25]
        cart_to_latt =  [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
        latt_to_cart =  [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
        eps =  1e-10
        out =  [0.5, 0.5, 0.25]
        self.assertEqual(bring_into_cell(bas_vecs,cart_to_latt,latt_to_cart,eps),out)

    def test_14(self):
        from phenum.symmetry import bring_into_cell
        bas_vecs =  [0.0, 0.0, 0.0]
        cart_to_latt =  [[-1.0, 1.0, 1.0], [1.0, 1.0, -1.0], [1.0, -1.0, 1.0]]
        latt_to_cart =  [[0.0, 0.5, 0.5], [0.5, 0.5, 0.0], [0.5, 0.0, 0.5]]
        eps =  1e-10
        out =  [0.0, 0.0, 0.0]
        self.assertEqual(bring_into_cell(bas_vecs,cart_to_latt,latt_to_cart,eps),out)

    def test_15(self):
        from phenum.symmetry import bring_into_cell
        bas_vecs =  [0.0, 0.0, 0.0]
        cart_to_latt =  [[-1.0, 1.0, 1.0], [1.0, 1.0, -1.0], [1.0, -1.0, 1.0]]
        latt_to_cart =  [[0.0, 0.5, 0.5], [0.5, 0.5, 0.0], [0.5, 0.0, 0.5]]
        eps =  1e-10
        out =  [0.0, 0.0, 0.0]
        self.assertEqual(bring_into_cell(bas_vecs,cart_to_latt,latt_to_cart,eps),out)

    def test_16(self):
        from phenum.symmetry import bring_into_cell
        bas_vecs = [1.0,0.0,1.0]
        cart_to_latt = np.transpose([[1.0, -1.0, 0.1], [1.0, 1.0, -0.1], [-1.0, 1.0, 0.1]])
        latt_to_cart = np.transpose([[0.5, 0.5, 0.0], [0.0, 0.5, 0.5], [5.0, 0.0, 5.0]])
        eps = 1e-3
        out = [1.0,0.0,1.0]
        self.assertEqual(bring_into_cell(bas_vecs,cart_to_latt,latt_to_cart,eps),out)

    def test_17(self):
        from phenum.symmetry import bring_into_cell
        bas_vecs = [2.0,0.0,2.0]
        cart_to_latt = np.transpose([[1.0, -1.0, 0.1], [1.0, 1.0, -0.1], [-1.0, 1.0, 0.1]])
        latt_to_cart = np.transpose([[0.5, 0.5, 0.0], [0.0, 0.5, 0.5], [5.0, 0.0, 5.0]])
        eps = 1e-3
        out = [2.0,0.0,2.0]
        self.assertEqual(bring_into_cell(bas_vecs,cart_to_latt,latt_to_cart,eps),out)

    def test_18(self):
        from phenum.symmetry import bring_into_cell
        bas_vecs = [3.0,0.0,3.0]
        cart_to_latt = np.transpose([[1.0, -1.0, 0.1], [1.0, 1.0, -0.1], [-1.0, 1.0, 0.1]])
        latt_to_cart = np.transpose([[0.5, 0.5, 0.0], [0.0, 0.5, 0.5], [5.0, 0.0, 5.0]])
        eps = 1e-3
        out = [3.0000000000000004, 0.0, 3.0000000000000004]
        self.assertEqual(bring_into_cell(bas_vecs,cart_to_latt,latt_to_cart,eps),out)

    def test_19(self):
        from phenum.symmetry import bring_into_cell
        bas_vecs = [0.0, 0.5779502399999998, 1.6329931599999998]
        cart_to_latt = np.transpose([[1.0, 0.0, 0.0], [-0.5773502717125849, 1.1547005434251698, 0.0], [0.0, 0.0, 0.6123724213915893]])
        latt_to_cart = np.transpose([[1.0, 0.0, 0.0], [0.5, 0.8660254, 0.0], [0.0, 0.0, 1.6329932]])
        eps = 1e-3
        out = [1.0000000000000000, 0.57795023999999980, -4.0000000416390380E-008]
        self.assertEqual(bring_into_cell(bas_vecs,cart_to_latt,latt_to_cart,eps),out)
       
class TestGetLatticePointGroup(ut.TestCase):
    def test_getpg1(self):
        from phenum.symmetry import get_lattice_pointGroup
        avecs = [[1,1,0],[1,0,1],[0,1,1]]
        eps = 1E-6
        out = _read_pg(gpath+"fcc_pg.out")

        pg = get_lattice_pointGroup(avecs,eps)
        present = []
        for p in pg:
            for g in out:
                if np.allclose(p,g):
                    present.append(True)
                    
        self.assertEqual(len(out),len(present))

    def test_getpg2(self):
        from phenum.symmetry import get_lattice_pointGroup
        avecs = [[1,0,0],[0.5,-0.86602540378444,0],[0,0,2]]
        eps = 1E-6
        out = _read_pg(gpath+"hex_pg.out")

        pg = get_lattice_pointGroup(avecs,eps)
        present = []
        for p in pg:
            for g in out:
                if np.allclose(p,g):
                    present.append(True)
                    
        self.assertEqual(len(out),len(present))

    def test_getpg3(self):
        from phenum.symmetry import get_lattice_pointGroup
        avecs = [[1,0,0],[0,1,0],[0,0,1]]
        eps = 1E-6
        out = _read_pg(gpath+"sc_pg.out")

        pg = get_lattice_pointGroup(avecs,eps)
        present = []
        for p in pg:
            for g in out:
                if np.allclose(p,g):
                    present.append(True)
                    
        self.assertEqual(len(out),len(present))

    def test_getpg4(self):
        from phenum.symmetry import get_lattice_pointGroup
        avecs = [[1,-1,1],[-1,1,1],[1,1,-1]]
        eps = 1E-6
        out = _read_pg(gpath+"bcc_pg.out")

        pg = get_lattice_pointGroup(avecs,eps)
        present = []
        for p in pg:
            for g in out:
                if np.allclose(p,g):
                    present.append(True)
                    
        self.assertEqual(len(out),len(present))

    def test_getpg5(self):
        from phenum.symmetry import get_lattice_pointGroup
        avecs = [[1,2,2],[2,1,2],[2,2,1]]
        eps = 1E-6
        out = _read_pg(gpath+"trig_pg.out")

        pg = get_lattice_pointGroup(avecs,eps)
        present = []
        for p in pg:
            for g in out:
                if np.allclose(p,g):
                    present.append(True)
                    
        self.assertEqual(len(out),len(present))

    def test_getpg6(self):
        from phenum.symmetry import get_lattice_pointGroup
        avecs = [[1,0,0],[0.15,1,0],[0.25,0,1]]
        eps = 1E-6
        out = _read_pg(gpath+"tric_pg.out")

        pg = get_lattice_pointGroup(avecs,eps)
        present = []
        for p in pg:
            for g in out:
                if np.allclose(p,g):
                    present.append(True)
                    
        self.assertEqual(len(out),len(present))

    def test_getpg7(self):
        from phenum.symmetry import get_lattice_pointGroup
        avecs = [[1,0,0],[0,1,0],[0,0,2]]
        eps = 1E-6
        out = _read_pg(gpath+"st_pg.out")

        pg = get_lattice_pointGroup(avecs,eps)
        present = []
        for p in pg:
            for g in out:
                if np.allclose(p,g):
                    present.append(True)
                    
        self.assertEqual(len(out),len(present))

    def test_getpg8(self):
        from phenum.symmetry import get_lattice_pointGroup
        avecs = [[-0.5,0.5,1],[0.5,-0.5,1],[0.5,0.5,-1]]
        eps = 1E-6
        out = _read_pg(gpath+"bct_pg.out")

        pg = get_lattice_pointGroup(avecs,eps)
        present = []
        for p in pg:
            for g in out:
                if np.allclose(p,g):
                    present.append(True)
                    
        self.assertEqual(len(out),len(present))

    def test_getpg9(self):
        from phenum.symmetry import get_lattice_pointGroup
        avecs = [[1,0,0],[0,2,0],[0,0,3]]
        eps = 1E-6
        out = _read_pg(gpath+"so_pg.out")

        pg = get_lattice_pointGroup(avecs,eps)
        present = []
        for p in pg:
            for g in out:
                if np.allclose(p,g):
                    present.append(True)
                    
        self.assertEqual(len(out),len(present))

    def test_getpg10(self):
        from phenum.symmetry import get_lattice_pointGroup
        avecs = [[0.5,1,0],[0.5,-1,0],[0,0,3]]
        eps = 1E-6
        out = _read_pg(gpath+"cco_pg.out")

        pg = get_lattice_pointGroup(avecs,eps)
        present = []
        for p in pg:
            for g in out:
                if np.allclose(p,g):
                    present.append(True)
                    
        self.assertEqual(len(out),len(present))

    def test_getpg11(self):
        from phenum.symmetry import get_lattice_pointGroup
        avecs = [[-0.5,1,1.5],[0.5,-1,1.5],[0.5,1,-1.5]]
        eps = 1E-6
        out = _read_pg(gpath+"bco_pg.out")

        pg = get_lattice_pointGroup(avecs,eps)
        present = []
        for p in pg:
            for g in out:
                if np.allclose(p,g):
                    present.append(True)
                    
        self.assertEqual(len(out),len(present))

    def test_getpg13(self):
        from phenum.symmetry import get_lattice_pointGroup
        avecs = [[0.5,1,0],[0.5,0,1.5],[0,1,1.5]]
        eps = 1E-6
        out = _read_pg(gpath+"fco_pg.out")

        pg = get_lattice_pointGroup(avecs,eps)
        present = []
        for p in pg:
            for g in out:
                if np.allclose(p,g):
                    present.append(True)
                    
        self.assertEqual(len(out),len(present))

    def test_getpg14(self):
        from phenum.symmetry import get_lattice_pointGroup
        avecs = [[1,0,0],[0,1,0],[0.25,0,1]]
        eps = 1E-6
        out = _read_pg(gpath+"sm_pg.out")

        pg = get_lattice_pointGroup(avecs,eps)
        present = []
        for p in pg:
            for g in out:
                if np.allclose(p,g):
                    present.append(True)
                    
        self.assertEqual(len(out),len(present))

    def test_getpg12(self):
        from phenum.symmetry import get_lattice_pointGroup
        avecs = [[0.5,0.5,0],[0.5,-0.5,0],[0.25,0,1]]
        eps = 1E-6
        out = _read_pg(gpath+"ccm_pg.out")

        pg = get_lattice_pointGroup(avecs,eps)
        present = []
        for p in pg:
            for g in out:
                if np.allclose(p,g):
                    present.append(True)
                    
        self.assertEqual(len(out),len(present))

    def test_getpg15(self):
        from phenum.symmetry import get_lattice_pointGroup
        avecs = [[1,1,0],[1,0,1],[0,0,2]]
        eps = 1E-6
        out = _read_pg(gpath+"fcc2_pg.out")

        pg = get_lattice_pointGroup(avecs,eps)
        present = []
        for p in pg:
            for g in out:
                if np.allclose(p,g):
                    present.append(True)
                    
        self.assertEqual(len(out),len(present))

class TestGetTransformations(ut.TestCase):
    """Tests of the _get_transformations subroutine."""

    def _compare_outputs(self,out1,out2):

        ptc1 = out1[0]
        ctp1 = out1[1]
        ptc2 = out2[0]
        ctp2 = out2[1]
        
        for i in range(3):
            for j in range(3):
                self.assertAlmostEqual(ptc1[i][j],ptc2[i][j])
                self.assertAlmostEqual(ctp1[i][j],ctp2[i][j])
    
    def _trans_out(self,case):
        ctp = np.transpose(_read_float_2D(gpath+"get_transformations_cart_to_prim.out."+str(case)))
        ptc = np.transpose(_read_float_2D(gpath+"get_transformations_prim_to_cart.out."+str(case)))
        return (ptc,ctp)
    
    def test_1(self):
        from phenum.symmetry import _get_transformations
        case = 1
        aVecs = _read_float_2D(gpath+"get_transformations_aVecs.in."+str(case))

        out = self._trans_out(case)
        self._compare_outputs(_get_transformations(aVecs),out)
    
    def test_2(self):
        from phenum.symmetry import _get_transformations
        case = 2
        aVecs = _read_float_2D(gpath+"get_transformations_aVecs.in."+str(case))

        out = self._trans_out(case)
        self._compare_outputs(_get_transformations(aVecs),out)
    
    def test_3(self):
        from phenum.symmetry import _get_transformations
        case = 3
        aVecs = _read_float_2D(gpath+"get_transformations_aVecs.in."+str(case))

        out = self._trans_out(case)
        self._compare_outputs(_get_transformations(aVecs),out)
    
    def test_4(self):
        from phenum.symmetry import _get_transformations
        case = 4
        aVecs = _read_float_2D(gpath+"get_transformations_aVecs.in."+str(case))

        out = self._trans_out(case)
        self._compare_outputs(_get_transformations(aVecs),out)
    
    def test_5(self):
        from phenum.symmetry import _get_transformations
        case = 5
        aVecs = _read_float_2D(gpath+"get_transformations_aVecs.in."+str(case))

        out = self._trans_out(case)
        self._compare_outputs(_get_transformations(aVecs),out)
    
    def test_6(self):
        from phenum.symmetry import _get_transformations
        case = 6
        aVecs = _read_float_2D(gpath+"get_transformations_aVecs.in."+str(case))

        out = self._trans_out(case)
        self._compare_outputs(_get_transformations(aVecs),out)
    
    def test_7(self):
        from phenum.symmetry import _get_transformations
        case = 7
        aVecs = _read_float_2D(gpath+"get_transformations_aVecs.in."+str(case))

        out = self._trans_out(case)
        self._compare_outputs(_get_transformations(aVecs),out)
    
    def test_8(self):
        from phenum.symmetry import _get_transformations
        case = 8
        aVecs = _read_float_2D(gpath+"get_transformations_aVecs.in."+str(case))

        out = self._trans_out(case)
        self._compare_outputs(_get_transformations(aVecs),out)
    
    def test_9(self):
        from phenum.symmetry import _get_transformations
        case = 9
        aVecs = _read_float_2D(gpath+"get_transformations_aVecs.in."+str(case))

        out = self._trans_out(case)
        self._compare_outputs(_get_transformations(aVecs),out)
    
    def test_10(self):
        from phenum.symmetry import _get_transformations
        case = 10
        aVecs = _read_float_2D(gpath+"get_transformations_aVecs.in."+str(case))

        out = self._trans_out(case)
        self._compare_outputs(_get_transformations(aVecs),out)
    
    def test_11(self):
        from phenum.symmetry import _get_transformations
        case = 11
        aVecs = _read_float_2D(gpath+"get_transformations_aVecs.in."+str(case))

        out = self._trans_out(case)
        self._compare_outputs(_get_transformations(aVecs),out)
    
    def test_12(self):
        from phenum.symmetry import _get_transformations
        case = 12
        aVecs = _read_float_2D(gpath+"get_transformations_aVecs.in."+str(case))

        out = self._trans_out(case)
        self._compare_outputs(_get_transformations(aVecs),out)
    
    def test_13(self):
        from phenum.symmetry import _get_transformations
        case = 13
        aVecs = _read_float_2D(gpath+"get_transformations_aVecs.in."+str(case))

        out = self._trans_out(case)
        self._compare_outputs(_get_transformations(aVecs),out)
    
    def test_14(self):
        from phenum.symmetry import _get_transformations
        case = 14
        aVecs = _read_float_2D(gpath+"get_transformations_aVecs.in."+str(case))

        out = self._trans_out(case)
        self._compare_outputs(_get_transformations(aVecs),out)
    
    def test_15(self):
        from phenum.symmetry import _get_transformations
        case = 15
        aVecs = _read_float_2D(gpath+"get_transformations_aVecs.in."+str(case))

        out = self._trans_out(case)
        self._compare_outputs(_get_transformations(aVecs),out)
    
    def test_16(self):
        from phenum.symmetry import _get_transformations
        case = 16
        aVecs = _read_float_2D(gpath+"get_transformations_aVecs.in."+str(case))

        out = self._trans_out(case)
        self._compare_outputs(_get_transformations(aVecs),out)
    
    def test_17(self):
        from phenum.symmetry import _get_transformations
        case = 17
        aVecs = _read_float_2D(gpath+"get_transformations_aVecs.in."+str(case))

        out = self._trans_out(case)
        self._compare_outputs(_get_transformations(aVecs),out)
    
    def test_18(self):
        from phenum.symmetry import _get_transformations
        case = 18
        aVecs = _read_float_2D(gpath+"get_transformations_aVecs.in."+str(case))

        out = self._trans_out(case)
        self._compare_outputs(_get_transformations(aVecs),out)
    
    def test_19(self):
        from phenum.symmetry import _get_transformations
        case = 19
        aVecs = _read_float_2D(gpath+"get_transformations_aVecs.in."+str(case))

        out = self._trans_out(case)
        self._compare_outputs(_get_transformations(aVecs),out)
    
    def test_20(self):
        from phenum.symmetry import _get_transformations
        case = 20
        aVecs = _read_float_2D(gpath+"get_transformations_aVecs.in."+str(case))

        out = self._trans_out(case)
        self._compare_outputs(_get_transformations(aVecs),out)
    
class TestDoesMappingExist(ut.TestCase):
    """Tests of the _does_mapping_exist subroutine."""

    def test_1(self):
        from phenum.symmetry import _does_mapping_exist
        case = 1
        v = _read_float_1D(gpath+"does_mapping_exist_v.in."+str(case))
        this_type = _read_int(gpath+"does_mapping_exist_this_type.in."+str(case))
        atom_pos = list(map(list,zip(*_read_float_2D(gpath+"does_mapping_exist_atom_pos.in."+str(case)))))
        atomType = _read_int_1D(gpath+"does_mapping_exist_atomType.in."+str(case))
        eps = _read_float(gpath+"does_mapping_exist_eps.in."+str(case))
        out = _read_logical(gpath+"does_mapping_exist_mapped.out."+str(case))
        self.assertEqual(_does_mapping_exist(v,this_type,atom_pos,atomType,eps),out)
        
    def test_2(self):
        from phenum.symmetry import _does_mapping_exist
        case = 2
        v = _read_float_1D(gpath+"does_mapping_exist_v.in."+str(case))
        this_type = _read_int(gpath+"does_mapping_exist_this_type.in."+str(case))
        atom_pos = list(map(list,zip(*_read_float_2D(gpath+"does_mapping_exist_atom_pos.in."+str(case)))))
        atomType = _read_int_1D(gpath+"does_mapping_exist_atomType.in."+str(case))
        eps = _read_float(gpath+"does_mapping_exist_eps.in."+str(case))
        out = _read_logical(gpath+"does_mapping_exist_mapped.out."+str(case))
        self.assertEqual(_does_mapping_exist(v,this_type,atom_pos,atomType,eps),out)
        
    def test_3(self):
        from phenum.symmetry import _does_mapping_exist
        case = 3
        v = _read_float_1D(gpath+"does_mapping_exist_v.in."+str(case))
        this_type = _read_int(gpath+"does_mapping_exist_this_type.in."+str(case))
        atom_pos = list(map(list,zip(*_read_float_2D(gpath+"does_mapping_exist_atom_pos.in."+str(case)))))
        atomType = _read_int_1D(gpath+"does_mapping_exist_atomType.in."+str(case))
        eps = _read_float(gpath+"does_mapping_exist_eps.in."+str(case))
        out = _read_logical(gpath+"does_mapping_exist_mapped.out."+str(case))
        self.assertEqual(_does_mapping_exist(v,this_type,atom_pos,atomType,eps),out)
        
    def test_4(self):
        from phenum.symmetry import _does_mapping_exist
        case = 4
        v = _read_float_1D(gpath+"does_mapping_exist_v.in."+str(case))
        this_type = _read_int(gpath+"does_mapping_exist_this_type.in."+str(case))
        atom_pos = list(map(list,zip(*_read_float_2D(gpath+"does_mapping_exist_atom_pos.in."+str(case)))))
        atomType = _read_int_1D(gpath+"does_mapping_exist_atomType.in."+str(case))
        eps = _read_float(gpath+"does_mapping_exist_eps.in."+str(case))
        out = _read_logical(gpath+"does_mapping_exist_mapped.out."+str(case))
        self.assertEqual(_does_mapping_exist(v,this_type,atom_pos,atomType,eps),out)
        
    def test_5(self):
        from phenum.symmetry import _does_mapping_exist
        case = 5
        v = _read_float_1D(gpath+"does_mapping_exist_v.in."+str(case))
        this_type = _read_int(gpath+"does_mapping_exist_this_type.in."+str(case))
        atom_pos = list(map(list,zip(*_read_float_2D(gpath+"does_mapping_exist_atom_pos.in."+str(case)))))
        atomType = _read_int_1D(gpath+"does_mapping_exist_atomType.in."+str(case))
        eps = _read_float(gpath+"does_mapping_exist_eps.in."+str(case))
        out = _read_logical(gpath+"does_mapping_exist_mapped.out."+str(case))
        self.assertEqual(_does_mapping_exist(v,this_type,atom_pos,atomType,eps),out)
        
    def test_6(self):
        from phenum.symmetry import _does_mapping_exist
        case = 6
        v = _read_float_1D(gpath+"does_mapping_exist_v.in."+str(case))
        this_type = _read_int(gpath+"does_mapping_exist_this_type.in."+str(case))
        atom_pos = list(map(list,zip(*_read_float_2D(gpath+"does_mapping_exist_atom_pos.in."+str(case)))))
        atomType = _read_int_1D(gpath+"does_mapping_exist_atomType.in."+str(case))
        eps = _read_float(gpath+"does_mapping_exist_eps.in."+str(case))
        out = _read_logical(gpath+"does_mapping_exist_mapped.out."+str(case))
        self.assertEqual(_does_mapping_exist(v,this_type,atom_pos,atomType,eps),out)
        
    def test_7(self):
        from phenum.symmetry import _does_mapping_exist
        case = 7
        v = _read_float_1D(gpath+"does_mapping_exist_v.in."+str(case))
        this_type = _read_int(gpath+"does_mapping_exist_this_type.in."+str(case))
        atom_pos = list(map(list,zip(*_read_float_2D(gpath+"does_mapping_exist_atom_pos.in."+str(case)))))
        atomType = _read_int_1D(gpath+"does_mapping_exist_atomType.in."+str(case))
        eps = _read_float(gpath+"does_mapping_exist_eps.in."+str(case))
        out = _read_logical(gpath+"does_mapping_exist_mapped.out."+str(case))
        self.assertEqual(_does_mapping_exist(v,this_type,atom_pos,atomType,eps),out)
        
    def test_8(self):
        from phenum.symmetry import _does_mapping_exist
        case = 8
        v = _read_float_1D(gpath+"does_mapping_exist_v.in."+str(case))
        this_type = _read_int(gpath+"does_mapping_exist_this_type.in."+str(case))
        atom_pos = list(map(list,zip(*_read_float_2D(gpath+"does_mapping_exist_atom_pos.in."+str(case)))))
        atomType = _read_int_1D(gpath+"does_mapping_exist_atomType.in."+str(case))
        eps = _read_float(gpath+"does_mapping_exist_eps.in."+str(case))
        out = _read_logical(gpath+"does_mapping_exist_mapped.out."+str(case))
        self.assertEqual(_does_mapping_exist(v,this_type,atom_pos,atomType,eps),out)
        
    def test_9(self):
        from phenum.symmetry import _does_mapping_exist
        case = 9
        v = _read_float_1D(gpath+"does_mapping_exist_v.in."+str(case))
        this_type = _read_int(gpath+"does_mapping_exist_this_type.in."+str(case))
        atom_pos = list(map(list,zip(*_read_float_2D(gpath+"does_mapping_exist_atom_pos.in."+str(case)))))
        atomType = _read_int_1D(gpath+"does_mapping_exist_atomType.in."+str(case))
        eps = _read_float(gpath+"does_mapping_exist_eps.in."+str(case))
        out = _read_logical(gpath+"does_mapping_exist_mapped.out."+str(case))
        self.assertEqual(_does_mapping_exist(v,this_type,atom_pos,atomType,eps),out)
        
    def test_10(self):
        from phenum.symmetry import _does_mapping_exist
        case = 10
        v = _read_float_1D(gpath+"does_mapping_exist_v.in."+str(case))
        this_type = _read_int(gpath+"does_mapping_exist_this_type.in."+str(case))
        atom_pos = list(map(list,zip(*_read_float_2D(gpath+"does_mapping_exist_atom_pos.in."+str(case)))))
        atomType = _read_int_1D(gpath+"does_mapping_exist_atomType.in."+str(case))
        eps = _read_float(gpath+"does_mapping_exist_eps.in."+str(case))
        out = _read_logical(gpath+"does_mapping_exist_mapped.out."+str(case))
        self.assertEqual(_does_mapping_exist(v,this_type,atom_pos,atomType,eps),out)
        
    def test_11(self):
        from phenum.symmetry import _does_mapping_exist
        case = 11
        v = _read_float_1D(gpath+"does_mapping_exist_v.in."+str(case))
        this_type = _read_int(gpath+"does_mapping_exist_this_type.in."+str(case))
        atom_pos = list(map(list,zip(*_read_float_2D(gpath+"does_mapping_exist_atom_pos.in."+str(case)))))
        atomType = _read_int_1D(gpath+"does_mapping_exist_atomType.in."+str(case))
        eps = _read_float(gpath+"does_mapping_exist_eps.in."+str(case))
        out = _read_logical(gpath+"does_mapping_exist_mapped.out."+str(case))
        self.assertEqual(_does_mapping_exist(v,this_type,atom_pos,atomType,eps),out)
        
    def test_12(self):
        from phenum.symmetry import _does_mapping_exist
        case = 12
        v = _read_float_1D(gpath+"does_mapping_exist_v.in."+str(case))
        this_type = _read_int(gpath+"does_mapping_exist_this_type.in."+str(case))
        atom_pos = list(map(list,zip(*_read_float_2D(gpath+"does_mapping_exist_atom_pos.in."+str(case)))))
        atomType = _read_int_1D(gpath+"does_mapping_exist_atomType.in."+str(case))
        eps = _read_float(gpath+"does_mapping_exist_eps.in."+str(case))
        out = _read_logical(gpath+"does_mapping_exist_mapped.out."+str(case))
        self.assertEqual(_does_mapping_exist(v,this_type,atom_pos,atomType,eps),out)
        
    def test_13(self):
        from phenum.symmetry import _does_mapping_exist
        case = 13
        v = _read_float_1D(gpath+"does_mapping_exist_v.in."+str(case))
        this_type = _read_int(gpath+"does_mapping_exist_this_type.in."+str(case))
        atom_pos = list(map(list,zip(*_read_float_2D(gpath+"does_mapping_exist_atom_pos.in."+str(case)))))
        atomType = _read_int_1D(gpath+"does_mapping_exist_atomType.in."+str(case))
        eps = _read_float(gpath+"does_mapping_exist_eps.in."+str(case))
        out = _read_logical(gpath+"does_mapping_exist_mapped.out."+str(case))
        self.assertEqual(_does_mapping_exist(v,this_type,atom_pos,atomType,eps),out)
        
    def test_14(self):
        from phenum.symmetry import _does_mapping_exist
        case = 14
        v = _read_float_1D(gpath+"does_mapping_exist_v.in."+str(case))
        this_type = _read_int(gpath+"does_mapping_exist_this_type.in."+str(case))
        atom_pos = list(map(list,zip(*_read_float_2D(gpath+"does_mapping_exist_atom_pos.in."+str(case)))))
        atomType = _read_int_1D(gpath+"does_mapping_exist_atomType.in."+str(case))
        eps = _read_float(gpath+"does_mapping_exist_eps.in."+str(case))
        out = _read_logical(gpath+"does_mapping_exist_mapped.out."+str(case))
        self.assertEqual(_does_mapping_exist(v,this_type,atom_pos,atomType,eps),out)
        
    def test_15(self):
        from phenum.symmetry import _does_mapping_exist
        case = 15
        v = _read_float_1D(gpath+"does_mapping_exist_v.in."+str(case))
        this_type = _read_int(gpath+"does_mapping_exist_this_type.in."+str(case))
        atom_pos = list(map(list,zip(*_read_float_2D(gpath+"does_mapping_exist_atom_pos.in."+str(case)))))
        atomType = _read_int_1D(gpath+"does_mapping_exist_atomType.in."+str(case))
        eps = _read_float(gpath+"does_mapping_exist_eps.in."+str(case))
        out = _read_logical(gpath+"does_mapping_exist_mapped.out."+str(case))
        self.assertEqual(_does_mapping_exist(v,this_type,atom_pos,atomType,eps),out)
        
    def test_16(self):
        from phenum.symmetry import _does_mapping_exist
        case = 16
        v = _read_float_1D(gpath+"does_mapping_exist_v.in."+str(case))
        this_type = _read_int(gpath+"does_mapping_exist_this_type.in."+str(case))
        atom_pos = list(map(list,zip(*_read_float_2D(gpath+"does_mapping_exist_atom_pos.in."+str(case)))))
        atomType = _read_int_1D(gpath+"does_mapping_exist_atomType.in."+str(case))
        eps = _read_float(gpath+"does_mapping_exist_eps.in."+str(case))
        out = _read_logical(gpath+"does_mapping_exist_mapped.out."+str(case))
        self.assertEqual(_does_mapping_exist(v,this_type,atom_pos,atomType,eps),out)
        
    def test_17(self):
        from phenum.symmetry import _does_mapping_exist
        case = 17
        v = _read_float_1D(gpath+"does_mapping_exist_v.in."+str(case))
        this_type = _read_int(gpath+"does_mapping_exist_this_type.in."+str(case))
        atom_pos = list(map(list,zip(*_read_float_2D(gpath+"does_mapping_exist_atom_pos.in."+str(case)))))
        atomType = _read_int_1D(gpath+"does_mapping_exist_atomType.in."+str(case))
        eps = _read_float(gpath+"does_mapping_exist_eps.in."+str(case))
        out = _read_logical(gpath+"does_mapping_exist_mapped.out."+str(case))
        self.assertEqual(_does_mapping_exist(v,this_type,atom_pos,atomType,eps),out)
        
    def test_18(self):
        from phenum.symmetry import _does_mapping_exist
        case = 18
        v = _read_float_1D(gpath+"does_mapping_exist_v.in."+str(case))
        this_type = _read_int(gpath+"does_mapping_exist_this_type.in."+str(case))
        atom_pos = list(map(list,zip(*_read_float_2D(gpath+"does_mapping_exist_atom_pos.in."+str(case)))))
        atomType = _read_int_1D(gpath+"does_mapping_exist_atomType.in."+str(case))
        eps = _read_float(gpath+"does_mapping_exist_eps.in."+str(case))
        out = _read_logical(gpath+"does_mapping_exist_mapped.out."+str(case))
        self.assertEqual(_does_mapping_exist(v,this_type,atom_pos,atomType,eps),out)
        
    def test_19(self):
        from phenum.symmetry import _does_mapping_exist
        case = 19
        v = _read_float_1D(gpath+"does_mapping_exist_v.in."+str(case))
        this_type = _read_int(gpath+"does_mapping_exist_this_type.in."+str(case))
        atom_pos = list(map(list,zip(*_read_float_2D(gpath+"does_mapping_exist_atom_pos.in."+str(case)))))
        atomType = _read_int_1D(gpath+"does_mapping_exist_atomType.in."+str(case))
        eps = _read_float(gpath+"does_mapping_exist_eps.in."+str(case))
        out = _read_logical(gpath+"does_mapping_exist_mapped.out."+str(case))
        self.assertEqual(_does_mapping_exist(v,this_type,atom_pos,atomType,eps),out)
        
    def test_20(self):
        from phenum.symmetry import _does_mapping_exist
        case = 20
        v = _read_float_1D(gpath+"does_mapping_exist_v.in."+str(case))
        this_type = _read_int(gpath+"does_mapping_exist_this_type.in."+str(case))
        atom_pos = list(map(list,zip(*_read_float_2D(gpath+"does_mapping_exist_atom_pos.in."+str(case)))))
        atomType = _read_int_1D(gpath+"does_mapping_exist_atomType.in."+str(case))
        eps = _read_float(gpath+"does_mapping_exist_eps.in."+str(case))
        out = _read_logical(gpath+"does_mapping_exist_mapped.out."+str(case))
        self.assertEqual(_does_mapping_exist(v,this_type,atom_pos,atomType,eps),out)
        
