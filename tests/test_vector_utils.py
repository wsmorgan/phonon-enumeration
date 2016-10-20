"""Methods for testing the subroutines in the vector_utils module."""
import unittest as ut
import pytest

class TestMapEnumStrToRealSpace(ut.TestCase):
    """Tests of the map_enumStr_to_real_space subroutine."""

    def _read_array_int_2d(self,fname):
        array = []
        with open(fname,"r") as f1:
            for line in f1:
                if "#" not in line:
                    array.append([int(i) for i in line.strip().split()])
        return array

    def _read_array_int_1d(self,fname):
        with open(fname,"r") as f1:
            line = f1.readline()
            if "#" in line:
                line = f1.readline()
            
        array = [int(i) for i in line.strip().split()]
        return array

    def _read_int(self,fname):
        with open(fname,"r") as f1:
            line = f1.readline()
            if "#" in line:
                line = f1.readline()

        return int(line.strip())

    def _read_array_float(self,fname):
        with open(fname,"r") as f1:
            line = f1.readline()
            if "#" in line:
                line = f1.readline()
            
        array = [float(i) for i in line.strip().split()]
        return array

    def _read_float(self,fname):
        with open(fname,"r") as f1:
            line = f1.readline()
            if "#" in line:
                line = f1.readline()

        return float(line.strip())
    
    def _read_array_float_2d(self,fname):
        array = []
        with open(fname,"r") as f1:
            for line in f1:
                if "#" not in line:
                    array.append([float(i) for i in line.strip().split()])
        return array
    
    def _read_array_float_trans(self,fname):
        array = []
        with open(fname,"r") as f1:
            for line in f1:
                if "#" not in line:
                    array.append([float(i) for i in line.strip().split()])
        return list(map(list,zip(*array)))

    def _read_labeling(self,fname):
        with open(fname,"r") as f1:
            line = f1.readline()
            if "#" in line:
                line = f1.readline()

        return line.strip()
                
    def _read_system_data(self,case):
        path = "tests/vector_utils/map_enumStr_to_real_space_"
        system_data={'bulksurf': 'bulk',
                     'k': self._read_int(path+"k.in."+str(case)),
                     'nD': 0,
                     'eps': self._read_float(path+"eps.in."+str(case)),
                     'title': 'Random structure enumeration: 2016-07-19 16:09:58.264761',
                     'plattice': self._read_array_float_2d(path+"pLV.in."+str(case)),
                     'dvecs': self._read_array_float_trans(path+"pBas.in."+str(case))}
        system_data["nD"] = len(system_data["dvecs"])
        return system_data

    def _read_struct_data(self,case):
        path = "tests/vector_utils/map_enumStr_to_real_space_"
        structure_data={'pgOps': 1,
                        'hnf_degen': 1,
                        'hnfN': 1,
                        'HNF': self._read_array_int_2d(path+"HNF.in."+str(case)),
                        'L': self._read_array_int_2d(path+"L.in."+str(case)),
                        'labeling': self._read_labeling(path+"labeling.in."+str(case)),
                        'n': self._read_int(path+"n.in."+str(case)),
                        'strN': 10,
                        'diag': self._read_array_int_1d(path+"S.in."+str(case)),
                        'directions': '000',
                        'sizeN': 5,
                        'tot_degen': 1}
        return structure_data
    
    def _read_mink(self,fname):
        with open(fname,"r") as f1:
            line = f1.readline()
            if "#" in line:
                line = f1.readline()

        if "t" in line.lower():
            out = True
        else:
            out = False
        return out

    def _read_space_data(self,case):
        path = "tests/vector_utils/map_enumStr_to_real_space_"
        space_data ={'aBas': self._read_array_float_trans(path+"aBas.out."+str(case)),
                     'sLV': self._read_array_float_trans(path+"sLV.out."+str(case)),
                     'gIndx': [i-1 for i in self._read_array_int_1d(path+"gIndx.out."+str(case))],
                     'spin':  self._read_array_int_1d(path+"spin.out."+str(case)),
                     'x': self._read_array_float(path+"x.out."+str(case))}
        return space_data

    def test_map1(self):
        from phenum.vector_utils import map_enumStr_to_real_space
        case = 1
        path = "tests/vector_utils/"
        minkowskiReduce = self._read_mink(path + "map_enumStr_to_real_space_minkowskiReduce.in."+str(case))
        system_data = self._read_system_data(case)
        structure_data = self._read_struct_data(case)
        out_test = map_enumStr_to_real_space(system_data,structure_data,minkowskiReduce)
        out = self._read_space_data(case)
        self.assertEqual(out,out_test)

    def test_map2(self):
        from phenum.vector_utils import map_enumStr_to_real_space
        case = 2
        path = "tests/vector_utils/"
        minkowskiReduce = self._read_mink(path + "map_enumStr_to_real_space_minkowskiReduce.in."+str(case))
        system_data = self._read_system_data(case)
        structure_data = self._read_struct_data(case)
        out_test = map_enumStr_to_real_space(system_data,structure_data,minkowskiReduce)
        out = self._read_space_data(case)
        self.assertEqual(out,out_test)

    def test_map3(self):
        from phenum.vector_utils import map_enumStr_to_real_space
        case = 3
        path = "tests/vector_utils/"
        minkowskiReduce = self._read_mink(path + "map_enumStr_to_real_space_minkowskiReduce.in."+str(case))
        system_data = self._read_system_data(case)
        structure_data = self._read_struct_data(case)
        out_test = map_enumStr_to_real_space(system_data,structure_data,minkowskiReduce)
        out = self._read_space_data(case)
        self.assertEqual(out,out_test)

    def test_map4(self):
        from phenum.vector_utils import map_enumStr_to_real_space
        case = 4
        path = "tests/vector_utils/"
        minkowskiReduce = self._read_mink(path + "map_enumStr_to_real_space_minkowskiReduce.in."+str(case))
        system_data = self._read_system_data(case)
        structure_data = self._read_struct_data(case)
        out_test = map_enumStr_to_real_space(system_data,structure_data,minkowskiReduce)
        out = self._read_space_data(case)
        self.assertEqual(out,out_test)

    def test_map5(self):
        from phenum.vector_utils import map_enumStr_to_real_space
        case = 5
        path = "tests/vector_utils/"
        minkowskiReduce = self._read_mink(path + "map_enumStr_to_real_space_minkowskiReduce.in."+str(case))
        system_data = self._read_system_data(case)
        structure_data = self._read_struct_data(case)
        out_test = map_enumStr_to_real_space(system_data,structure_data,minkowskiReduce)
        out = self._read_space_data(case)
        self.assertEqual(out,out_test)

    def test_map6(self):
        from phenum.vector_utils import map_enumStr_to_real_space
        case = 6
        path = "tests/vector_utils/"
        minkowskiReduce = self._read_mink(path + "map_enumStr_to_real_space_minkowskiReduce.in."+str(case))
        system_data = self._read_system_data(case)
        structure_data = self._read_struct_data(case)
        out_test = map_enumStr_to_real_space(system_data,structure_data,minkowskiReduce)
        out = self._read_space_data(case)
        self.assertEqual(out,out_test)

    def test_map7(self):
        from phenum.vector_utils import map_enumStr_to_real_space
        case = 7
        path = "tests/vector_utils/"
        minkowskiReduce = self._read_mink(path + "map_enumStr_to_real_space_minkowskiReduce.in."+str(case))
        system_data = self._read_system_data(case)
        structure_data = self._read_struct_data(case)
        out_test = map_enumStr_to_real_space(system_data,structure_data,minkowskiReduce)
        out = self._read_space_data(case)
        self.assertEqual(out,out_test)

    def test_map8(self):
        from phenum.vector_utils import map_enumStr_to_real_space
        case = 8
        path = "tests/vector_utils/"
        minkowskiReduce = self._read_mink(path + "map_enumStr_to_real_space_minkowskiReduce.in."+str(case))
        system_data = self._read_system_data(case)
        structure_data = self._read_struct_data(case)
        out_test = map_enumStr_to_real_space(system_data,structure_data,minkowskiReduce)
        out = self._read_space_data(case)
        self.assertEqual(out,out_test)

    def test_map9(self):
        from phenum.vector_utils import map_enumStr_to_real_space
        case = 9
        path = "tests/vector_utils/"
        minkowskiReduce = self._read_mink(path + "map_enumStr_to_real_space_minkowskiReduce.in."+str(case))
        system_data = self._read_system_data(case)
        structure_data = self._read_struct_data(case)
        out_test = map_enumStr_to_real_space(system_data,structure_data,minkowskiReduce)
        out = self._read_space_data(case)
        self.assertEqual(out,out_test)

    def test_map10(self):
        from phenum.vector_utils import map_enumStr_to_real_space
        case = 10
        path = "tests/vector_utils/"
        minkowskiReduce = self._read_mink(path + "map_enumStr_to_real_space_minkowskiReduce.in."+str(case))
        system_data = self._read_system_data(case)
        structure_data = self._read_struct_data(case)
        out_test = map_enumStr_to_real_space(system_data,structure_data,minkowskiReduce)
        out = self._read_space_data(case)
        self.assertEqual(out,out_test)

    def test_map11(self):
        from phenum.vector_utils import map_enumStr_to_real_space
        case = 11
        path = "tests/vector_utils/"
        minkowskiReduce = self._read_mink(path + "map_enumStr_to_real_space_minkowskiReduce.in."+str(case))
        system_data = self._read_system_data(case)
        structure_data = self._read_struct_data(case)
        out_test = map_enumStr_to_real_space(system_data,structure_data,minkowskiReduce)
        out = self._read_space_data(case)
        self.assertEqual(out,out_test)

    def test_map12(self):
        from phenum.vector_utils import map_enumStr_to_real_space
        case = 12
        path = "tests/vector_utils/"
        minkowskiReduce = self._read_mink(path + "map_enumStr_to_real_space_minkowskiReduce.in."+str(case))
        system_data = self._read_system_data(case)
        structure_data = self._read_struct_data(case)
        out_test = map_enumStr_to_real_space(system_data,structure_data,minkowskiReduce)
        out = self._read_space_data(case)
        self.assertEqual(out,out_test)

    def test_map13(self):
        from phenum.vector_utils import map_enumStr_to_real_space
        case = 13
        path = "tests/vector_utils/"
        minkowskiReduce = self._read_mink(path + "map_enumStr_to_real_space_minkowskiReduce.in."+str(case))
        system_data = self._read_system_data(case)
        structure_data = self._read_struct_data(case)
        out_test = map_enumStr_to_real_space(system_data,structure_data,minkowskiReduce)
        out = self._read_space_data(case)
        self.assertEqual(out,out_test)

    def test_map14(self):
        from phenum.vector_utils import map_enumStr_to_real_space
        case = 14
        path = "tests/vector_utils/"
        minkowskiReduce = self._read_mink(path + "map_enumStr_to_real_space_minkowskiReduce.in."+str(case))
        system_data = self._read_system_data(case)
        structure_data = self._read_struct_data(case)
        out_test = map_enumStr_to_real_space(system_data,structure_data,minkowskiReduce)
        out = self._read_space_data(case)
        self.assertEqual(out,out_test)

    def test_map15(self):
        from phenum.vector_utils import map_enumStr_to_real_space
        case = 15
        path = "tests/vector_utils/"
        minkowskiReduce = self._read_mink(path + "map_enumStr_to_real_space_minkowskiReduce.in."+str(case))
        system_data = self._read_system_data(case)
        structure_data = self._read_struct_data(case)
        out_test = map_enumStr_to_real_space(system_data,structure_data,minkowskiReduce)
        out = self._read_space_data(case)
        self.assertEqual(out,out_test)

    def test_map16(self):
        from phenum.vector_utils import map_enumStr_to_real_space
        case = 16
        path = "tests/vector_utils/"
        minkowskiReduce = self._read_mink(path + "map_enumStr_to_real_space_minkowskiReduce.in."+str(case))
        system_data = self._read_system_data(case)
        structure_data = self._read_struct_data(case)
        out_test = map_enumStr_to_real_space(system_data,structure_data,minkowskiReduce)
        out = self._read_space_data(case)
        self.assertEqual(out,out_test)

    def test_map17(self):
        from phenum.vector_utils import map_enumStr_to_real_space
        case = 17
        path = "tests/vector_utils/"
        minkowskiReduce = self._read_mink(path + "map_enumStr_to_real_space_minkowskiReduce.in."+str(case))
        system_data = self._read_system_data(case)
        structure_data = self._read_struct_data(case)
        out_test = map_enumStr_to_real_space(system_data,structure_data,minkowskiReduce)
        out = self._read_space_data(case)
        self.assertEqual(out,out_test)

    def test_map18(self):
        from phenum.vector_utils import map_enumStr_to_real_space
        case = 18
        path = "tests/vector_utils/"
        minkowskiReduce = self._read_mink(path + "map_enumStr_to_real_space_minkowskiReduce.in."+str(case))
        system_data = self._read_system_data(case)
        structure_data = self._read_struct_data(case)
        out_test = map_enumStr_to_real_space(system_data,structure_data,minkowskiReduce)
        out = self._read_space_data(case)
        self.assertEqual(out,out_test)

    def test_map19(self):
        from phenum.vector_utils import map_enumStr_to_real_space
        case = 19
        path = "tests/vector_utils/"
        minkowskiReduce = self._read_mink(path + "map_enumStr_to_real_space_minkowskiReduce.in."+str(case))
        system_data = self._read_system_data(case)
        structure_data = self._read_struct_data(case)
        out_test = map_enumStr_to_real_space(system_data,structure_data,minkowskiReduce)
        out = self._read_space_data(case)
        self.assertEqual(out,out_test)

    def test_map20(self):
        from phenum.vector_utils import map_enumStr_to_real_space
        case = 20
        path = "tests/vector_utils/"
        minkowskiReduce = self._read_mink(path + "map_enumStr_to_real_space_minkowskiReduce.in."+str(case))
        system_data = self._read_system_data(case)
        structure_data = self._read_struct_data(case)
        out_test = map_enumStr_to_real_space(system_data,structure_data,minkowskiReduce)
        out = self._read_space_data(case)
        self.assertEqual(out,out_test)

    def test_map20(self):
        from phenum.vector_utils import map_enumStr_to_real_space
        case = 30
        path = "tests/vector_utils/"
        minkowskiReduce = self._read_mink(path + "map_enumStr_to_real_space_minkowskiReduce.in."+str(case))
        system_data = self._read_system_data(case)
        structure_data = self._read_struct_data(case)
        out_test = map_enumStr_to_real_space(system_data,structure_data,minkowskiReduce)
        out = self._read_space_data(case)
        self.assertEqual(out,out_test)
        
class TestReduceCInABC(ut.TestCase):
    """Tests of the reduce_C_in_ABC subroutine."""

    def _read_array_float(self,file):
        with open(file,"r") as f1:
            line = f1.readline()
            if "#" in line:
                line = f1.readline()
            
        array = [float(i) for i in line.strip().split()]
        return array

    def _read_float(self,file):
        with open(file,"r") as f1:
            line = f1.readline()
            if "#" in line:
                line = f1.readline()
        val = float(line.strip())
        return val

    def _read_output(self,case):
        path = "tests/vector_utils/"
        A = self._read_array_float(path + "reduce_C_in_ABC_A.out."+str(case))
        B = self._read_array_float(path +"reduce_C_in_ABC_B.out."+str(case))
        C = self._read_array_float(path+"reduce_C_in_ABC_C.out."+str(case))
        
        return [A, B, C]

    def _compare(self,out1,out2):
        for i in range(3):
            for j in range(3):
                self.assertAlmostEqual(out1[i][j],out2[i][j],places=12)
                
    def test_reduce1(self):
        from phenum.vector_utils import _reduce_C_in_ABC
        case = 1
        path = "tests/vector_utils/"
        A = self._read_array_float(path+"reduce_C_in_ABC_A.in."+str(case))
        B = self._read_array_float(path+"reduce_C_in_ABC_B.in."+str(case))
        C = self._read_array_float(path+"reduce_C_in_ABC_C.in."+str(case))
        eps = self._read_float(path+"reduce_C_in_ABC_eps.in."+str(case))

        AF, BF, CF = _reduce_C_in_ABC(A,B,C,eps)

        self._compare([AF,BF,CF],self._read_output(case))

    def test_reduce2(self):
        from phenum.vector_utils import _reduce_C_in_ABC
        case = 2
        path = "tests/vector_utils/"
        A = self._read_array_float(path+"reduce_C_in_ABC_A.in."+str(case))
        B = self._read_array_float(path+"reduce_C_in_ABC_B.in."+str(case))
        C = self._read_array_float(path+"reduce_C_in_ABC_C.in."+str(case))
        eps = self._read_float(path+"reduce_C_in_ABC_eps.in."+str(case))

        AF, BF, CF = _reduce_C_in_ABC(A,B,C,eps)

        self._compare([AF,BF,CF],self._read_output(case))

    def test_reduce3(self):
        from phenum.vector_utils import _reduce_C_in_ABC
        case = 3
        path = "tests/vector_utils/"
        A = self._read_array_float(path+"reduce_C_in_ABC_A.in."+str(case))
        B = self._read_array_float(path+"reduce_C_in_ABC_B.in."+str(case))
        C = self._read_array_float(path+"reduce_C_in_ABC_C.in."+str(case))
        eps = self._read_float(path+"reduce_C_in_ABC_eps.in."+str(case))

        AF, BF, CF = _reduce_C_in_ABC(A,B,C,eps)

        self._compare([AF,BF,CF],self._read_output(case))

    def test_reduce4(self):
        from phenum.vector_utils import _reduce_C_in_ABC
        case = 4
        path = "tests/vector_utils/"
        A = self._read_array_float(path+"reduce_C_in_ABC_A.in."+str(case))
        B = self._read_array_float(path+"reduce_C_in_ABC_B.in."+str(case))
        C = self._read_array_float(path+"reduce_C_in_ABC_C.in."+str(case))
        eps = self._read_float(path+"reduce_C_in_ABC_eps.in."+str(case))

        AF, BF, CF = _reduce_C_in_ABC(A,B,C,eps)

        self._compare([AF,BF,CF],self._read_output(case))

    def test_reduce5(self):
        from phenum.vector_utils import _reduce_C_in_ABC
        case = 5
        path = "tests/vector_utils/"
        A = self._read_array_float(path+"reduce_C_in_ABC_A.in."+str(case))
        B = self._read_array_float(path+"reduce_C_in_ABC_B.in."+str(case))
        C = self._read_array_float(path+"reduce_C_in_ABC_C.in."+str(case))
        eps = self._read_float(path+"reduce_C_in_ABC_eps.in."+str(case))

        AF, BF, CF = _reduce_C_in_ABC(A,B,C,eps)

        self._compare([AF,BF,CF],self._read_output(case))

    def test_reduce6(self):
        from phenum.vector_utils import _reduce_C_in_ABC
        case = 6
        path = "tests/vector_utils/"
        A = self._read_array_float(path+"reduce_C_in_ABC_A.in."+str(case))
        B = self._read_array_float(path+"reduce_C_in_ABC_B.in."+str(case))
        C = self._read_array_float(path+"reduce_C_in_ABC_C.in."+str(case))
        eps = self._read_float(path+"reduce_C_in_ABC_eps.in."+str(case))

        AF, BF, CF = _reduce_C_in_ABC(A,B,C,eps)

        self._compare([AF,BF,CF],self._read_output(case))

    def test_reduce7(self):
        from phenum.vector_utils import _reduce_C_in_ABC
        case = 7
        path = "tests/vector_utils/"
        A = self._read_array_float(path+"reduce_C_in_ABC_A.in."+str(case))
        B = self._read_array_float(path+"reduce_C_in_ABC_B.in."+str(case))
        C = self._read_array_float(path+"reduce_C_in_ABC_C.in."+str(case))
        eps = self._read_float(path+"reduce_C_in_ABC_eps.in."+str(case))

        AF, BF, CF = _reduce_C_in_ABC(A,B,C,eps)

        self._compare([AF,BF,CF],self._read_output(case))

    def test_reduce8(self):
        from phenum.vector_utils import _reduce_C_in_ABC
        case = 8
        path = "tests/vector_utils/"
        A = self._read_array_float(path+"reduce_C_in_ABC_A.in."+str(case))
        B = self._read_array_float(path+"reduce_C_in_ABC_B.in."+str(case))
        C = self._read_array_float(path+"reduce_C_in_ABC_C.in."+str(case))
        eps = self._read_float(path+"reduce_C_in_ABC_eps.in."+str(case))

        AF, BF, CF = _reduce_C_in_ABC(A,B,C,eps)

        self._compare([AF,BF,CF],self._read_output(case))

    def test_reduce9(self):
        from phenum.vector_utils import _reduce_C_in_ABC
        case = 9
        path = "tests/vector_utils/"
        A = self._read_array_float(path+"reduce_C_in_ABC_A.in."+str(case))
        B = self._read_array_float(path+"reduce_C_in_ABC_B.in."+str(case))
        C = self._read_array_float(path+"reduce_C_in_ABC_C.in."+str(case))
        eps = self._read_float(path+"reduce_C_in_ABC_eps.in."+str(case))

        AF, BF, CF = _reduce_C_in_ABC(A,B,C,eps)

        self._compare([AF,BF,CF],self._read_output(case))

    def test_reduce10(self):
        from phenum.vector_utils import _reduce_C_in_ABC
        case = 10
        path = "tests/vector_utils/"
        A = self._read_array_float(path+"reduce_C_in_ABC_A.in."+str(case))
        B = self._read_array_float(path+"reduce_C_in_ABC_B.in."+str(case))
        C = self._read_array_float(path+"reduce_C_in_ABC_C.in."+str(case))
        eps = self._read_float(path+"reduce_C_in_ABC_eps.in."+str(case))

        AF, BF, CF = _reduce_C_in_ABC(A,B,C,eps)

        self._compare([AF,BF,CF],self._read_output(case))

    def test_reduce11(self):
        from phenum.vector_utils import _reduce_C_in_ABC
        case = 11
        path = "tests/vector_utils/"
        A = self._read_array_float(path+"reduce_C_in_ABC_A.in."+str(case))
        B = self._read_array_float(path+"reduce_C_in_ABC_B.in."+str(case))
        C = self._read_array_float(path+"reduce_C_in_ABC_C.in."+str(case))
        eps = self._read_float(path+"reduce_C_in_ABC_eps.in."+str(case))

        AF, BF, CF = _reduce_C_in_ABC(A,B,C,eps)

        self._compare([AF,BF,CF],self._read_output(case))

    def test_reduce12(self):
        from phenum.vector_utils import _reduce_C_in_ABC
        case = 12
        path = "tests/vector_utils/"
        A = self._read_array_float(path+"reduce_C_in_ABC_A.in."+str(case))
        B = self._read_array_float(path+"reduce_C_in_ABC_B.in."+str(case))
        C = self._read_array_float(path+"reduce_C_in_ABC_C.in."+str(case))
        eps = self._read_float(path+"reduce_C_in_ABC_eps.in."+str(case))

        AF, BF, CF = _reduce_C_in_ABC(A,B,C,eps)

        self._compare([AF,BF,CF],self._read_output(case))

    def test_reduce13(self):
        from phenum.vector_utils import _reduce_C_in_ABC
        case = 13
        path = "tests/vector_utils/"
        A = self._read_array_float(path+"reduce_C_in_ABC_A.in."+str(case))
        B = self._read_array_float(path+"reduce_C_in_ABC_B.in."+str(case))
        C = self._read_array_float(path+"reduce_C_in_ABC_C.in."+str(case))
        eps = self._read_float(path+"reduce_C_in_ABC_eps.in."+str(case))

        AF, BF, CF = _reduce_C_in_ABC(A,B,C,eps)

        self._compare([AF,BF,CF],self._read_output(case))

    def test_reduce14(self):
        from phenum.vector_utils import _reduce_C_in_ABC
        case = 14
        path = "tests/vector_utils/"
        A = self._read_array_float(path+"reduce_C_in_ABC_A.in."+str(case))
        B = self._read_array_float(path+"reduce_C_in_ABC_B.in."+str(case))
        C = self._read_array_float(path+"reduce_C_in_ABC_C.in."+str(case))
        eps = self._read_float(path+"reduce_C_in_ABC_eps.in."+str(case))

        AF, BF, CF = _reduce_C_in_ABC(A,B,C,eps)

        self._compare([AF,BF,CF],self._read_output(case))

    def test_reduce15(self):
        from phenum.vector_utils import _reduce_C_in_ABC
        case = 15
        path = "tests/vector_utils/"
        A = self._read_array_float(path+"reduce_C_in_ABC_A.in."+str(case))
        B = self._read_array_float(path+"reduce_C_in_ABC_B.in."+str(case))
        C = self._read_array_float(path+"reduce_C_in_ABC_C.in."+str(case))
        eps = self._read_float(path+"reduce_C_in_ABC_eps.in."+str(case))

        AF, BF, CF = _reduce_C_in_ABC(A,B,C,eps)

        self._compare([AF,BF,CF],self._read_output(case))

    def test_reduce16(self):
        from phenum.vector_utils import _reduce_C_in_ABC
        case = 16
        path = "tests/vector_utils/"
        A = self._read_array_float(path+"reduce_C_in_ABC_A.in."+str(case))
        B = self._read_array_float(path+"reduce_C_in_ABC_B.in."+str(case))
        C = self._read_array_float(path+"reduce_C_in_ABC_C.in."+str(case))
        eps = self._read_float(path+"reduce_C_in_ABC_eps.in."+str(case))

        AF, BF, CF = _reduce_C_in_ABC(A,B,C,eps)

        self._compare([AF,BF,CF],self._read_output(case))

    def test_reduce17(self):
        from phenum.vector_utils import _reduce_C_in_ABC
        case = 17
        path = "tests/vector_utils/"
        A = self._read_array_float(path+"reduce_C_in_ABC_A.in."+str(case))
        B = self._read_array_float(path+"reduce_C_in_ABC_B.in."+str(case))
        C = self._read_array_float(path+"reduce_C_in_ABC_C.in."+str(case))
        eps = self._read_float(path+"reduce_C_in_ABC_eps.in."+str(case))

        AF, BF, CF = _reduce_C_in_ABC(A,B,C,eps)

        self._compare([AF,BF,CF],self._read_output(case))

    def test_reduce18(self):
        from phenum.vector_utils import _reduce_C_in_ABC
        case = 18
        path = "tests/vector_utils/"
        A = self._read_array_float(path+"reduce_C_in_ABC_A.in."+str(case))
        B = self._read_array_float(path+"reduce_C_in_ABC_B.in."+str(case))
        C = self._read_array_float(path+"reduce_C_in_ABC_C.in."+str(case))
        eps = self._read_float(path+"reduce_C_in_ABC_eps.in."+str(case))

        AF, BF, CF = _reduce_C_in_ABC(A,B,C,eps)

        self._compare([AF,BF,CF],self._read_output(case))

    def test_reduce19(self):
        from phenum.vector_utils import _reduce_C_in_ABC
        case = 19
        path = "tests/vector_utils/"
        A = self._read_array_float(path+"reduce_C_in_ABC_A.in."+str(case))
        B = self._read_array_float(path+"reduce_C_in_ABC_B.in."+str(case))
        C = self._read_array_float(path+"reduce_C_in_ABC_C.in."+str(case))
        eps = self._read_float(path+"reduce_C_in_ABC_eps.in."+str(case))

        AF, BF, CF = _reduce_C_in_ABC(A,B,C,eps)

        self._compare([AF,BF,CF],self._read_output(case))

    def test_reduce20(self):
        from phenum.vector_utils import _reduce_C_in_ABC
        case = 20
        path = "tests/vector_utils/"
        A = self._read_array_float(path+"reduce_C_in_ABC_A.in."+str(case))
        B = self._read_array_float(path+"reduce_C_in_ABC_B.in."+str(case))
        C = self._read_array_float(path+"reduce_C_in_ABC_C.in."+str(case))
        eps = self._read_float(path+"reduce_C_in_ABC_eps.in."+str(case))

        AF, BF, CF = _reduce_C_in_ABC(A,B,C,eps)

        self._compare([AF,BF,CF],self._read_output(case))

class TestGaussionReduceTwoVectors(ut.TestCase):
    """Tests of the _gaussion_reduce_two_vectors subroutine."""
    

    def _read_array_float(self,file):
        with open(file,"r") as f1:
            line = f1.readline()
            if "#" in line:
                line = f1.readline()
            
        array = [float(i) for i in line.strip().split()]
        return array

    def _read_float(self,file):
        with open(file,"r") as f1:
            line = f1.readline()
            if "#" in line:
                line = f1.readline()
        val = float(line.strip())
        return val

    def _read_output(self,case):
        path = "tests/vector_utils/"
        V = self._read_array_float(path + "gaussian_reduce_two_vectors_V.out."+str(case))
        U = self._read_array_float(path +"gaussian_reduce_two_vectors_U.out."+str(case))
        
        return [U, V]

    def _compare(self,out1,out2):
        for i in range(2):
            for j in range(3):
                self.assertAlmostEqual(out1[i][j],out2[i][j],places=12)

    def test_gaus1(self):
        from phenum.vector_utils import _gaussian_reduce_two_vectors
        case = 1
        path = "tests/vector_utils/"

        V = self._read_array_float(path + "gaussian_reduce_two_vectors_V.in."+str(case))
        U = self._read_array_float(path + "gaussian_reduce_two_vectors_U.in."+str(case))
        eps = self._read_float(path + "gaussian_reduce_two_vectors_eps.in."+str(case))
        
        UF, VF = _gaussian_reduce_two_vectors(U,V,eps)

        self._compare([UF,VF],self._read_output(1))

    def test_gaus2(self):
        from phenum.vector_utils import _gaussian_reduce_two_vectors
        case = 2
        path = "tests/vector_utils/"

        V = self._read_array_float(path + "gaussian_reduce_two_vectors_V.in."+str(case))
        U = self._read_array_float(path + "gaussian_reduce_two_vectors_U.in."+str(case))
        eps = self._read_float(path + "gaussian_reduce_two_vectors_eps.in."+str(case))
        
        UF, VF = _gaussian_reduce_two_vectors(U,V,eps)

        self._compare([UF,VF],self._read_output(case))

    def test_gaus2(self):
        from phenum.vector_utils import _gaussian_reduce_two_vectors
        case = 2
        path = "tests/vector_utils/"

        V = self._read_array_float(path + "gaussian_reduce_two_vectors_V.in."+str(case))
        U = self._read_array_float(path + "gaussian_reduce_two_vectors_U.in."+str(case))
        eps = self._read_float(path + "gaussian_reduce_two_vectors_eps.in."+str(case))
        
        UF, VF = _gaussian_reduce_two_vectors(U,V,eps)

        self._compare([UF,VF],self._read_output(case))

    def test_gaus3(self):
        from phenum.vector_utils import _gaussian_reduce_two_vectors
        case = 3
        path = "tests/vector_utils/"

        V = self._read_array_float(path + "gaussian_reduce_two_vectors_V.in."+str(case))
        U = self._read_array_float(path + "gaussian_reduce_two_vectors_U.in."+str(case))
        eps = self._read_float(path + "gaussian_reduce_two_vectors_eps.in."+str(case))
        
        UF, VF = _gaussian_reduce_two_vectors(U,V,eps)

        self._compare([UF,VF],self._read_output(case))

    def test_gaus4(self):
        from phenum.vector_utils import _gaussian_reduce_two_vectors
        case = 4
        path = "tests/vector_utils/"

        V = self._read_array_float(path + "gaussian_reduce_two_vectors_V.in."+str(case))
        U = self._read_array_float(path + "gaussian_reduce_two_vectors_U.in."+str(case))
        eps = self._read_float(path + "gaussian_reduce_two_vectors_eps.in."+str(case))
        
        UF, VF = _gaussian_reduce_two_vectors(U,V,eps)

        self._compare([UF,VF],self._read_output(case))

    def test_gaus5(self):
        from phenum.vector_utils import _gaussian_reduce_two_vectors
        case = 5
        path = "tests/vector_utils/"

        V = self._read_array_float(path + "gaussian_reduce_two_vectors_V.in."+str(case))
        U = self._read_array_float(path + "gaussian_reduce_two_vectors_U.in."+str(case))
        eps = self._read_float(path + "gaussian_reduce_two_vectors_eps.in."+str(case))
        
        UF, VF = _gaussian_reduce_two_vectors(U,V,eps)

        self._compare([UF,VF],self._read_output(case))

    def test_gaus6(self):
        from phenum.vector_utils import _gaussian_reduce_two_vectors
        case = 6
        path = "tests/vector_utils/"

        V = self._read_array_float(path + "gaussian_reduce_two_vectors_V.in."+str(case))
        U = self._read_array_float(path + "gaussian_reduce_two_vectors_U.in."+str(case))
        eps = self._read_float(path + "gaussian_reduce_two_vectors_eps.in."+str(case))
        
        UF, VF = _gaussian_reduce_two_vectors(U,V,eps)

        self._compare([UF,VF],self._read_output(case))

    def test_gaus7(self):
        from phenum.vector_utils import _gaussian_reduce_two_vectors
        case = 7
        path = "tests/vector_utils/"

        V = self._read_array_float(path + "gaussian_reduce_two_vectors_V.in."+str(case))
        U = self._read_array_float(path + "gaussian_reduce_two_vectors_U.in."+str(case))
        eps = self._read_float(path + "gaussian_reduce_two_vectors_eps.in."+str(case))
        
        UF, VF = _gaussian_reduce_two_vectors(U,V,eps)

        self._compare([UF,VF],self._read_output(case))

    def test_gaus8(self):
        from phenum.vector_utils import _gaussian_reduce_two_vectors
        case = 8
        path = "tests/vector_utils/"

        V = self._read_array_float(path + "gaussian_reduce_two_vectors_V.in."+str(case))
        U = self._read_array_float(path + "gaussian_reduce_two_vectors_U.in."+str(case))
        eps = self._read_float(path + "gaussian_reduce_two_vectors_eps.in."+str(case))
        
        UF, VF = _gaussian_reduce_two_vectors(U,V,eps)

        self._compare([UF,VF],self._read_output(case))

    def test_gaus9(self):
        from phenum.vector_utils import _gaussian_reduce_two_vectors
        case = 9
        path = "tests/vector_utils/"

        V = self._read_array_float(path + "gaussian_reduce_two_vectors_V.in."+str(case))
        U = self._read_array_float(path + "gaussian_reduce_two_vectors_U.in."+str(case))
        eps = self._read_float(path + "gaussian_reduce_two_vectors_eps.in."+str(case))
        
        UF, VF = _gaussian_reduce_two_vectors(U,V,eps)

        self._compare([UF,VF],self._read_output(case))

    def test_gaus10(self):
        from phenum.vector_utils import _gaussian_reduce_two_vectors
        case = 10
        path = "tests/vector_utils/"

        V = self._read_array_float(path + "gaussian_reduce_two_vectors_V.in."+str(case))
        U = self._read_array_float(path + "gaussian_reduce_two_vectors_U.in."+str(case))
        eps = self._read_float(path + "gaussian_reduce_two_vectors_eps.in."+str(case))
        
        UF, VF = _gaussian_reduce_two_vectors(U,V,eps)

        self._compare([UF,VF],self._read_output(case))

    def test_gaus11(self):
        from phenum.vector_utils import _gaussian_reduce_two_vectors
        case = 11
        path = "tests/vector_utils/"

        V = self._read_array_float(path + "gaussian_reduce_two_vectors_V.in."+str(case))
        U = self._read_array_float(path + "gaussian_reduce_two_vectors_U.in."+str(case))
        eps = self._read_float(path + "gaussian_reduce_two_vectors_eps.in."+str(case))
        
        UF, VF = _gaussian_reduce_two_vectors(U,V,eps)

        self._compare([UF,VF],self._read_output(case))

    def test_gaus12(self):
        from phenum.vector_utils import _gaussian_reduce_two_vectors
        case = 12
        path = "tests/vector_utils/"

        V = self._read_array_float(path + "gaussian_reduce_two_vectors_V.in."+str(case))
        U = self._read_array_float(path + "gaussian_reduce_two_vectors_U.in."+str(case))
        eps = self._read_float(path + "gaussian_reduce_two_vectors_eps.in."+str(case))
        
        UF, VF = _gaussian_reduce_two_vectors(U,V,eps)

        self._compare([UF,VF],self._read_output(case))

    def test_gaus13(self):
        from phenum.vector_utils import _gaussian_reduce_two_vectors
        case = 13
        path = "tests/vector_utils/"

        V = self._read_array_float(path + "gaussian_reduce_two_vectors_V.in."+str(case))
        U = self._read_array_float(path + "gaussian_reduce_two_vectors_U.in."+str(case))
        eps = self._read_float(path + "gaussian_reduce_two_vectors_eps.in."+str(case))
        
        UF, VF = _gaussian_reduce_two_vectors(U,V,eps)

        self._compare([UF,VF],self._read_output(case))

    def test_gaus14(self):
        from phenum.vector_utils import _gaussian_reduce_two_vectors
        case = 14
        path = "tests/vector_utils/"

        V = self._read_array_float(path + "gaussian_reduce_two_vectors_V.in."+str(case))
        U = self._read_array_float(path + "gaussian_reduce_two_vectors_U.in."+str(case))
        eps = self._read_float(path + "gaussian_reduce_two_vectors_eps.in."+str(case))
        
        UF, VF = _gaussian_reduce_two_vectors(U,V,eps)

        self._compare([UF,VF],self._read_output(case))

    def test_gaus15(self):
        from phenum.vector_utils import _gaussian_reduce_two_vectors
        case = 15
        path = "tests/vector_utils/"

        V = self._read_array_float(path + "gaussian_reduce_two_vectors_V.in."+str(case))
        U = self._read_array_float(path + "gaussian_reduce_two_vectors_U.in."+str(case))
        eps = self._read_float(path + "gaussian_reduce_two_vectors_eps.in."+str(case))
        
        UF, VF = _gaussian_reduce_two_vectors(U,V,eps)

        self._compare([UF,VF],self._read_output(case))

    def test_gaus16(self):
        from phenum.vector_utils import _gaussian_reduce_two_vectors
        case = 16
        path = "tests/vector_utils/"

        V = self._read_array_float(path + "gaussian_reduce_two_vectors_V.in."+str(case))
        U = self._read_array_float(path + "gaussian_reduce_two_vectors_U.in."+str(case))
        eps = self._read_float(path + "gaussian_reduce_two_vectors_eps.in."+str(case))
        
        UF, VF = _gaussian_reduce_two_vectors(U,V,eps)

        self._compare([UF,VF],self._read_output(case))

    def test_gaus17(self):
        from phenum.vector_utils import _gaussian_reduce_two_vectors
        case = 17
        path = "tests/vector_utils/"

        V = self._read_array_float(path + "gaussian_reduce_two_vectors_V.in."+str(case))
        U = self._read_array_float(path + "gaussian_reduce_two_vectors_U.in."+str(case))
        eps = self._read_float(path + "gaussian_reduce_two_vectors_eps.in."+str(case))
        
        UF, VF = _gaussian_reduce_two_vectors(U,V,eps)

        self._compare([UF,VF],self._read_output(case))

    def test_gaus18(self):
        from phenum.vector_utils import _gaussian_reduce_two_vectors
        case = 18
        path = "tests/vector_utils/"

        V = self._read_array_float(path + "gaussian_reduce_two_vectors_V.in."+str(case))
        U = self._read_array_float(path + "gaussian_reduce_two_vectors_U.in."+str(case))
        eps = self._read_float(path + "gaussian_reduce_two_vectors_eps.in."+str(case))
        
        UF, VF = _gaussian_reduce_two_vectors(U,V,eps)

        self._compare([UF,VF],self._read_output(case))

    def test_gaus19(self):
        from phenum.vector_utils import _gaussian_reduce_two_vectors
        case = 19
        path = "tests/vector_utils/"

        V = self._read_array_float(path + "gaussian_reduce_two_vectors_V.in."+str(case))
        U = self._read_array_float(path + "gaussian_reduce_two_vectors_U.in."+str(case))
        eps = self._read_float(path + "gaussian_reduce_two_vectors_eps.in."+str(case))
        
        UF, VF = _gaussian_reduce_two_vectors(U,V,eps)

        self._compare([UF,VF],self._read_output(case))

    def test_gaus20(self):
        from phenum.vector_utils import _gaussian_reduce_two_vectors
        case = 20
        path = "tests/vector_utils/"

        V = self._read_array_float(path + "gaussian_reduce_two_vectors_V.in."+str(case))
        U = self._read_array_float(path + "gaussian_reduce_two_vectors_U.in."+str(case))
        eps = self._read_float(path + "gaussian_reduce_two_vectors_eps.in."+str(case))
        
        UF, VF = _gaussian_reduce_two_vectors(U,V,eps)

        self._compare([UF,VF],self._read_output(case))

    def test_gaus21(self):
        from phenum.vector_utils import _gaussian_reduce_two_vectors
        case = 21
        path = "tests/vector_utils/"

        V = self._read_array_float(path + "gaussian_reduce_two_vectors_V.in."+str(case))
        U = self._read_array_float(path + "gaussian_reduce_two_vectors_U.in."+str(case))
        eps = self._read_float(path + "gaussian_reduce_two_vectors_eps.in."+str(case))
        
        UF, VF = _gaussian_reduce_two_vectors(U,V,eps)

        self._compare([UF,VF],self._read_output(case))

class TestMinkowskiReduceBasis(ut.TestCase):
    """Tests of the _minkowski_reduce_basis subroutine."""

    def _read_array_float(self,file):
        array = []        
        with open(file,"r") as f1:
            for line in f1:
                if "#" not in line:
                    array.append([float(i) for i in line.strip().split()])
        return list(map(list,zip(*array)))

    def _read_float(self,file):
        with open(file,"r") as f1:
            line = f1.readline()
            if "#" in line:
                line = f1.readline()
        val = float(line.strip())
        return val

    def _read_output(self,case):
        path = "tests/vector_utils/"
        out = self._read_array_float(path + "minkowski_reduce_basis_OUT.out."+str(case))
        
        return out

    def _compare(self,out1,out2):
        for i in range(3):
            for j in range(3):
                self.assertAlmostEqual(out1[i][j],out2[i][j],places=12)

    def test_minkr1(self):
        from phenum.vector_utils import _minkowski_reduce_basis
        case = 1
        path = "tests/vector_utils/"

        inp = self._read_array_float(path+"minkowski_reduce_basis_IN.in."+str(case))
        eps = self._read_float(path+"minkowski_reduce_basis_eps.in."+str(case))

        out = self._read_output(case)

        self._compare(_minkowski_reduce_basis(inp,eps),out)

    def test_minkr2(self):
        from phenum.vector_utils import _minkowski_reduce_basis
        case = 2
        path = "tests/vector_utils/"

        inp = self._read_array_float(path+"minkowski_reduce_basis_IN.in."+str(case))
        eps = self._read_float(path+"minkowski_reduce_basis_eps.in."+str(case))

        out = self._read_output(case)

        self._compare(_minkowski_reduce_basis(inp,eps),out)

    def test_minkr3(self):
        from phenum.vector_utils import _minkowski_reduce_basis
        case = 3
        path = "tests/vector_utils/"

        inp = self._read_array_float(path+"minkowski_reduce_basis_IN.in."+str(case))
        eps = self._read_float(path+"minkowski_reduce_basis_eps.in."+str(case))

        out = self._read_output(case)

        self._compare(_minkowski_reduce_basis(inp,eps),out)

    def test_minkr4(self):
        from phenum.vector_utils import _minkowski_reduce_basis
        case = 4
        path = "tests/vector_utils/"

        inp = self._read_array_float(path+"minkowski_reduce_basis_IN.in."+str(case))
        eps = self._read_float(path+"minkowski_reduce_basis_eps.in."+str(case))

        out = self._read_output(case)

        self._compare(_minkowski_reduce_basis(inp,eps),out)

    def test_minkr5(self):
        from phenum.vector_utils import _minkowski_reduce_basis
        case = 5
        path = "tests/vector_utils/"

        inp = self._read_array_float(path+"minkowski_reduce_basis_IN.in."+str(case))
        eps = self._read_float(path+"minkowski_reduce_basis_eps.in."+str(case))

        out = self._read_output(case)

        self._compare(_minkowski_reduce_basis(inp,eps),out)

    def test_minkr6(self):
        from phenum.vector_utils import _minkowski_reduce_basis
        case = 6
        path = "tests/vector_utils/"

        inp = self._read_array_float(path+"minkowski_reduce_basis_IN.in."+str(case))
        eps = self._read_float(path+"minkowski_reduce_basis_eps.in."+str(case))

        out = self._read_output(case)

        self._compare(_minkowski_reduce_basis(inp,eps),out)

    def test_minkr7(self):
        from phenum.vector_utils import _minkowski_reduce_basis
        case = 7
        path = "tests/vector_utils/"

        inp = self._read_array_float(path+"minkowski_reduce_basis_IN.in."+str(case))
        eps = self._read_float(path+"minkowski_reduce_basis_eps.in."+str(case))

        out = self._read_output(case)

        self._compare(_minkowski_reduce_basis(inp,eps),out)

    def test_minkr8(self):
        from phenum.vector_utils import _minkowski_reduce_basis
        case = 8
        path = "tests/vector_utils/"

        inp = self._read_array_float(path+"minkowski_reduce_basis_IN.in."+str(case))
        eps = self._read_float(path+"minkowski_reduce_basis_eps.in."+str(case))

        out = self._read_output(case)

        self._compare(_minkowski_reduce_basis(inp,eps),out)

    def test_minkr9(self):
        from phenum.vector_utils import _minkowski_reduce_basis
        case = 9
        path = "tests/vector_utils/"

        inp = self._read_array_float(path+"minkowski_reduce_basis_IN.in."+str(case))
        eps = self._read_float(path+"minkowski_reduce_basis_eps.in."+str(case))

        out = self._read_output(case)

        self._compare(_minkowski_reduce_basis(inp,eps),out)

    def test_minkr10(self):
        from phenum.vector_utils import _minkowski_reduce_basis
        case = 10
        path = "tests/vector_utils/"

        inp = self._read_array_float(path+"minkowski_reduce_basis_IN.in."+str(case))
        eps = self._read_float(path+"minkowski_reduce_basis_eps.in."+str(case))

        out = self._read_output(case)

        self._compare(_minkowski_reduce_basis(inp,eps),out)

    def test_minkr11(self):
        from phenum.vector_utils import _minkowski_reduce_basis
        case = 11
        path = "tests/vector_utils/"

        inp = self._read_array_float(path+"minkowski_reduce_basis_IN.in."+str(case))
        eps = self._read_float(path+"minkowski_reduce_basis_eps.in."+str(case))

        out = self._read_output(case)

        self._compare(_minkowski_reduce_basis(inp,eps),out)

    def test_minkr12(self):
        from phenum.vector_utils import _minkowski_reduce_basis
        case = 12
        path = "tests/vector_utils/"

        inp = self._read_array_float(path+"minkowski_reduce_basis_IN.in."+str(case))
        eps = self._read_float(path+"minkowski_reduce_basis_eps.in."+str(case))

        out = self._read_output(case)

        self._compare(_minkowski_reduce_basis(inp,eps),out)

    def test_minkr13(self):
        from phenum.vector_utils import _minkowski_reduce_basis
        case = 13
        path = "tests/vector_utils/"

        inp = self._read_array_float(path+"minkowski_reduce_basis_IN.in."+str(case))
        eps = self._read_float(path+"minkowski_reduce_basis_eps.in."+str(case))

        out = self._read_output(case)

        self._compare(_minkowski_reduce_basis(inp,eps),out)

    def test_minkr14(self):
        from phenum.vector_utils import _minkowski_reduce_basis
        case = 14
        path = "tests/vector_utils/"

        inp = self._read_array_float(path+"minkowski_reduce_basis_IN.in."+str(case))
        eps = self._read_float(path+"minkowski_reduce_basis_eps.in."+str(case))

        out = self._read_output(case)

        self._compare(_minkowski_reduce_basis(inp,eps),out)

    def test_minkr15(self):
        from phenum.vector_utils import _minkowski_reduce_basis
        case = 15
        path = "tests/vector_utils/"

        inp = self._read_array_float(path+"minkowski_reduce_basis_IN.in."+str(case))
        eps = self._read_float(path+"minkowski_reduce_basis_eps.in."+str(case))

        out = self._read_output(case)

        self._compare(_minkowski_reduce_basis(inp,eps),out)

    def test_minkr16(self):
        from phenum.vector_utils import _minkowski_reduce_basis
        case = 16
        path = "tests/vector_utils/"

        inp = self._read_array_float(path+"minkowski_reduce_basis_IN.in."+str(case))
        eps = self._read_float(path+"minkowski_reduce_basis_eps.in."+str(case))

        out = self._read_output(case)

        self._compare(_minkowski_reduce_basis(inp,eps),out)

    def test_minkr17(self):
        from phenum.vector_utils import _minkowski_reduce_basis
        case = 17
        path = "tests/vector_utils/"

        inp = self._read_array_float(path+"minkowski_reduce_basis_IN.in."+str(case))
        eps = self._read_float(path+"minkowski_reduce_basis_eps.in."+str(case))

        out = self._read_output(case)

        self._compare(_minkowski_reduce_basis(inp,eps),out)

    def test_minkr18(self):
        from phenum.vector_utils import _minkowski_reduce_basis
        case = 18
        path = "tests/vector_utils/"

        inp = self._read_array_float(path+"minkowski_reduce_basis_IN.in."+str(case))
        eps = self._read_float(path+"minkowski_reduce_basis_eps.in."+str(case))

        out = self._read_output(case)

        self._compare(_minkowski_reduce_basis(inp,eps),out)

    def test_minkr19(self):
        from phenum.vector_utils import _minkowski_reduce_basis
        case = 19
        path = "tests/vector_utils/"

        inp = self._read_array_float(path+"minkowski_reduce_basis_IN.in."+str(case))
        eps = self._read_float(path+"minkowski_reduce_basis_eps.in."+str(case))

        out = self._read_output(case)

        self._compare(_minkowski_reduce_basis(inp,eps),out)

    def test_minkr20(self):
        from phenum.vector_utils import _minkowski_reduce_basis
        case = 20
        path = "tests/vector_utils/"

        inp = self._read_array_float(path+"minkowski_reduce_basis_IN.in."+str(case))
        eps = self._read_float(path+"minkowski_reduce_basis_eps.in."+str(case))

        out = self._read_output(case)

        self._compare(_minkowski_reduce_basis(inp,eps),out)

    def test_minkr21(self):
        from phenum.vector_utils import _minkowski_reduce_basis
        case = 21
        path = "tests/vector_utils/"

        inp = self._read_array_float(path+"minkowski_reduce_basis_IN.in."+str(case))
        eps = self._read_float(path+"minkowski_reduce_basis_eps.in."+str(case))

        out = self._read_output(case)

        self._compare(_minkowski_reduce_basis(inp,eps),out)

    def test_minkr22(self):
        from phenum.vector_utils import _minkowski_reduce_basis
        case = 21
        path = "tests/vector_utils/"

        inp = [[1,0,0],[1,0,0],[1,0,0]]
        eps = self._read_float(path+"minkowski_reduce_basis_eps.in."+str(case))

        out = self._read_output(case)
        with pytest.raises(ValueError):
            self._compare(_minkowski_reduce_basis(inp,eps),out)
        
class TestMinkowskiConditionsCheck(ut.TestCase):
    """Tests of the _minkowski_conditions_check subroutine."""

    def _read_array_float(self,file):
        array = []
        
        with open(file,"r") as f1:
            for line in f1:
                if "#" not in line:
                    array.append([float(i) for i in line.strip().split()])
                
        return list(map(list,zip(*array)))

    def _read_float(self,file):
        with open(file,"r") as f1:
            for line in f1:
                if "#" in line:
                    pass
        
        val = float(line.strip())
        return val

    def _read_output(self,case):
        path = "tests/vector_utils/"
        with open(path + "minkowski_condition_check.out."+str(case),"r") as f1:
            line = f1.readline()
            if "#" in line:
                line = f1.readline()

        if line.strip() == "T":
            out = True
        else:
            out = False
            
        return out

    def test_minkc1(self):
        from phenum.vector_utils import _minkowski_conditions_check
        path = "tests/vector_utils/"
        case = 1

        basis = self._read_array_float(path+"minkowski_condition_check_basis.in."+str(case))
        eps = self._read_float(path+"minkowski_condition_check_eps.in."+str(case))

        out = self._read_output(case)

        self.assertEqual(out,_minkowski_conditions_check(basis,eps))

    def test_minkc2(self):
        from phenum.vector_utils import _minkowski_conditions_check
        path = "tests/vector_utils/"
        case = 2

        basis = self._read_array_float(path+"minkowski_condition_check_basis.in."+str(case))
        eps = self._read_float(path+"minkowski_condition_check_eps.in."+str(case))

        out = self._read_output(case)

        self.assertEqual(out,_minkowski_conditions_check(basis,eps))

    def test_minkc3(self):
        from phenum.vector_utils import _minkowski_conditions_check
        path = "tests/vector_utils/"
        case = 3

        basis = self._read_array_float(path+"minkowski_condition_check_basis.in."+str(case))
        eps = self._read_float(path+"minkowski_condition_check_eps.in."+str(case))

        out = self._read_output(case)

        self.assertEqual(out,_minkowski_conditions_check(basis,eps))

    def test_minkc4(self):
        from phenum.vector_utils import _minkowski_conditions_check
        path = "tests/vector_utils/"
        case = 4

        basis = self._read_array_float(path+"minkowski_condition_check_basis.in."+str(case))
        eps = self._read_float(path+"minkowski_condition_check_eps.in."+str(case))

        out = self._read_output(case)

        self.assertEqual(out,_minkowski_conditions_check(basis,eps))

    def test_minkc5(self):
        from phenum.vector_utils import _minkowski_conditions_check
        path = "tests/vector_utils/"
        case = 5

        basis = self._read_array_float(path+"minkowski_condition_check_basis.in."+str(case))
        eps = self._read_float(path+"minkowski_condition_check_eps.in."+str(case))

        out = self._read_output(case)

        self.assertEqual(out,_minkowski_conditions_check(basis,eps))

    def test_minkc6(self):
        from phenum.vector_utils import _minkowski_conditions_check
        path = "tests/vector_utils/"
        case = 6

        basis = self._read_array_float(path+"minkowski_condition_check_basis.in."+str(case))
        eps = self._read_float(path+"minkowski_condition_check_eps.in."+str(case))

        out = self._read_output(case)

        self.assertEqual(out,_minkowski_conditions_check(basis,eps))

    def test_minkc7(self):
        from phenum.vector_utils import _minkowski_conditions_check
        path = "tests/vector_utils/"
        case = 7

        basis = self._read_array_float(path+"minkowski_condition_check_basis.in."+str(case))
        eps = self._read_float(path+"minkowski_condition_check_eps.in."+str(case))

        out = self._read_output(case)

        self.assertEqual(out,_minkowski_conditions_check(basis,eps))

    def test_minkc8(self):
        from phenum.vector_utils import _minkowski_conditions_check
        path = "tests/vector_utils/"
        case = 8

        basis = self._read_array_float(path+"minkowski_condition_check_basis.in."+str(case))
        eps = self._read_float(path+"minkowski_condition_check_eps.in."+str(case))

        out = self._read_output(case)

        self.assertEqual(out,_minkowski_conditions_check(basis,eps))

    def test_minkc9(self):
        from phenum.vector_utils import _minkowski_conditions_check
        path = "tests/vector_utils/"
        case = 9

        basis = self._read_array_float(path+"minkowski_condition_check_basis.in."+str(case))
        eps = self._read_float(path+"minkowski_condition_check_eps.in."+str(case))

        out = self._read_output(case)

        self.assertEqual(out,_minkowski_conditions_check(basis,eps))

    def test_minkc10(self):
        from phenum.vector_utils import _minkowski_conditions_check
        path = "tests/vector_utils/"
        case = 10

        basis = self._read_array_float(path+"minkowski_condition_check_basis.in."+str(case))
        eps = self._read_float(path+"minkowski_condition_check_eps.in."+str(case))

        out = self._read_output(case)

        self.assertEqual(out,_minkowski_conditions_check(basis,eps))

    def test_minkc11(self):
        from phenum.vector_utils import _minkowski_conditions_check
        path = "tests/vector_utils/"
        case = 11

        basis = self._read_array_float(path+"minkowski_condition_check_basis.in."+str(case))
        eps = self._read_float(path+"minkowski_condition_check_eps.in."+str(case))

        out = self._read_output(case)

        self.assertEqual(out,_minkowski_conditions_check(basis,eps))

    def test_minkc12(self):
        from phenum.vector_utils import _minkowski_conditions_check
        path = "tests/vector_utils/"
        case = 12

        basis = self._read_array_float(path+"minkowski_condition_check_basis.in."+str(case))
        eps = self._read_float(path+"minkowski_condition_check_eps.in."+str(case))

        out = self._read_output(case)

        self.assertEqual(out,_minkowski_conditions_check(basis,eps))

    def test_minkc13(self):
        from phenum.vector_utils import _minkowski_conditions_check
        path = "tests/vector_utils/"
        case = 13

        basis = self._read_array_float(path+"minkowski_condition_check_basis.in."+str(case))
        eps = self._read_float(path+"minkowski_condition_check_eps.in."+str(case))

        out = self._read_output(case)

        self.assertEqual(out,_minkowski_conditions_check(basis,eps))

    def test_minkc14(self):
        from phenum.vector_utils import _minkowski_conditions_check
        path = "tests/vector_utils/"
        case = 14

        basis = self._read_array_float(path+"minkowski_condition_check_basis.in."+str(case))
        eps = self._read_float(path+"minkowski_condition_check_eps.in."+str(case))

        out = self._read_output(case)

        self.assertEqual(out,_minkowski_conditions_check(basis,eps))

    def test_minkc15(self):
        from phenum.vector_utils import _minkowski_conditions_check
        path = "tests/vector_utils/"
        case = 15

        basis = self._read_array_float(path+"minkowski_condition_check_basis.in."+str(case))
        eps = self._read_float(path+"minkowski_condition_check_eps.in."+str(case))

        out = self._read_output(case)

        self.assertEqual(out,_minkowski_conditions_check(basis,eps))

    def test_minkc16(self):
        from phenum.vector_utils import _minkowski_conditions_check
        path = "tests/vector_utils/"
        case = 16

        basis = self._read_array_float(path+"minkowski_condition_check_basis.in."+str(case))
        eps = self._read_float(path+"minkowski_condition_check_eps.in."+str(case))

        out = self._read_output(case)

        self.assertEqual(out,_minkowski_conditions_check(basis,eps))

    def test_minkc17(self):
        from phenum.vector_utils import _minkowski_conditions_check
        path = "tests/vector_utils/"
        case = 17

        basis = self._read_array_float(path+"minkowski_condition_check_basis.in."+str(case))
        eps = self._read_float(path+"minkowski_condition_check_eps.in."+str(case))

        out = self._read_output(case)

        self.assertEqual(out,_minkowski_conditions_check(basis,eps))

    def test_minkc18(self):
        from phenum.vector_utils import _minkowski_conditions_check
        path = "tests/vector_utils/"
        case = 18

        basis = self._read_array_float(path+"minkowski_condition_check_basis.in."+str(case))
        eps = self._read_float(path+"minkowski_condition_check_eps.in."+str(case))

        out = self._read_output(case)

        self.assertEqual(out,_minkowski_conditions_check(basis,eps))

    def test_minkc19(self):
        from phenum.vector_utils import _minkowski_conditions_check
        path = "tests/vector_utils/"
        case = 19

        basis = self._read_array_float(path+"minkowski_condition_check_basis.in."+str(case))
        eps = self._read_float(path+"minkowski_condition_check_eps.in."+str(case))

        out = self._read_output(case)

        self.assertEqual(out,_minkowski_conditions_check(basis,eps))

    def test_minkc20(self):
        from phenum.vector_utils import _minkowski_conditions_check
        path = "tests/vector_utils/"
        case = 20

        basis = self._read_array_float(path+"minkowski_condition_check_basis.in."+str(case))
        eps = self._read_float(path+"minkowski_condition_check_eps.in."+str(case))

        out = self._read_output(case)

        self.assertEqual(out,_minkowski_conditions_check(basis,eps))

    def test_minkc21(self):
        from phenum.vector_utils import _minkowski_conditions_check

        basis = [[3,3,3],[1,2,3],[0,0,1]]
        eps = 0.0
        out = False

        self.assertEqual(out,_minkowski_conditions_check(basis,eps))

    def test_minkc22(self):
        from phenum.vector_utils import _minkowski_conditions_check

        basis = [[0,0,-1],[3,3,3],[1,2,3]]
        eps = 0.0
        out = False

        self.assertEqual(out,_minkowski_conditions_check(basis,eps))

    def test_minkc23(self):
        from phenum.vector_utils import _minkowski_conditions_check

        basis = [[3,1,0],[3,2,0],[3,3,1]]
        eps = 0.0
        out = False

        self.assertEqual(out,_minkowski_conditions_check(basis,eps))

    def test_minkc24(self):
        from phenum.vector_utils import _minkowski_conditions_check

        basis = [[0,0,-1],[0,0,-1],[3,3,3]]
        eps = 0.0
        out = False

        self.assertEqual(out,_minkowski_conditions_check(basis,eps))

    def test_minkc25(self):
        from phenum.vector_utils import _minkowski_conditions_check

        basis = [[0,0,1],[0,0,-1],[3,3,3]]
        eps = 0.0
        out = False

        self.assertEqual(out,_minkowski_conditions_check(basis,eps))
        
class TestCartesian2Direct(ut.TestCase):
    """Test the cartesian2direct subroutine."""

    def _read_array_float(self,file):
        array = []
        with open(file,"r") as f1:
            for line in f1:
                if "#" not in line:
                    array.append([float(i) for i in line.strip().split()])
        return list(map(list,zip(*array)))

    def _read_float(self,file):
        with open(file,"r") as f1:
            line = f1.readline()
            if "#" in line:
                line = f1.readline()
        val = float(line.strip())
        return val

    def _compare(self,out1,out2):
        if len(out1) == len(out2):
            for i in range(len(out1)):
                for j in range(3):
                    self.assertAlmostEqual(out1[i][j],out2[i][j],places=12)
        else:
            self.assertEqual(len(out1),len(out2))

    def test_c2d1(self):
        from phenum.vector_utils import cartesian2direct
        case = 1
        path = "tests/vector_utils/"
        
        sLV = self._read_array_float(path+"cartesian2direct_sLV.in."+str(case))
        aBas = self._read_array_float(path+"cartesian2direct_aBas.in."+str(case))
        eps = self._read_float(path+"cartesian2direct_eps.in."+str(case))

        out = self._read_array_float(path+"cartesian2direct_aBas.out."+str(case))

        self._compare(out,cartesian2direct(sLV,aBas,eps))

    def test_c2d2(self):
        from phenum.vector_utils import cartesian2direct
        case = 2
        path = "tests/vector_utils/"
        
        sLV = self._read_array_float(path+"cartesian2direct_sLV.in."+str(case))
        aBas = self._read_array_float(path+"cartesian2direct_aBas.in."+str(case))
        eps = self._read_float(path+"cartesian2direct_eps.in."+str(case))

        out = self._read_array_float(path+"cartesian2direct_aBas.out."+str(case))

        self._compare(out,cartesian2direct(sLV,aBas,eps))

    def test_c2d3(self):
        from phenum.vector_utils import cartesian2direct
        case = 3
        path = "tests/vector_utils/"
        
        sLV = self._read_array_float(path+"cartesian2direct_sLV.in."+str(case))
        aBas = self._read_array_float(path+"cartesian2direct_aBas.in."+str(case))
        eps = self._read_float(path+"cartesian2direct_eps.in."+str(case))

        out = self._read_array_float(path+"cartesian2direct_aBas.out."+str(case))

        self._compare(out,cartesian2direct(sLV,aBas,eps))

    def test_c2d4(self):
        from phenum.vector_utils import cartesian2direct
        case = 4
        path = "tests/vector_utils/"
        
        sLV = self._read_array_float(path+"cartesian2direct_sLV.in."+str(case))
        aBas = self._read_array_float(path+"cartesian2direct_aBas.in."+str(case))
        eps = self._read_float(path+"cartesian2direct_eps.in."+str(case))

        out = self._read_array_float(path+"cartesian2direct_aBas.out."+str(case))

        self._compare(out,cartesian2direct(sLV,aBas,eps))

    def test_c2d5(self):
        from phenum.vector_utils import cartesian2direct
        case = 5
        path = "tests/vector_utils/"
        
        sLV = self._read_array_float(path+"cartesian2direct_sLV.in."+str(case))
        aBas = self._read_array_float(path+"cartesian2direct_aBas.in."+str(case))
        eps = self._read_float(path+"cartesian2direct_eps.in."+str(case))

        out = self._read_array_float(path+"cartesian2direct_aBas.out."+str(case))

        self._compare(out,cartesian2direct(sLV,aBas,eps))

    def test_c2d6(self):
        from phenum.vector_utils import cartesian2direct
        case = 6
        path = "tests/vector_utils/"
        
        sLV = self._read_array_float(path+"cartesian2direct_sLV.in."+str(case))
        aBas = self._read_array_float(path+"cartesian2direct_aBas.in."+str(case))
        eps = self._read_float(path+"cartesian2direct_eps.in."+str(case))

        out = self._read_array_float(path+"cartesian2direct_aBas.out."+str(case))

        self._compare(out,cartesian2direct(sLV,aBas,eps))

    def test_c2d7(self):
        from phenum.vector_utils import cartesian2direct
        case = 7
        path = "tests/vector_utils/"
        
        sLV = self._read_array_float(path+"cartesian2direct_sLV.in."+str(case))
        aBas = self._read_array_float(path+"cartesian2direct_aBas.in."+str(case))
        eps = self._read_float(path+"cartesian2direct_eps.in."+str(case))

        out = self._read_array_float(path+"cartesian2direct_aBas.out."+str(case))

        self._compare(out,cartesian2direct(sLV,aBas,eps))

    def test_c2d8(self):
        from phenum.vector_utils import cartesian2direct
        case = 8
        path = "tests/vector_utils/"
        
        sLV = self._read_array_float(path+"cartesian2direct_sLV.in."+str(case))
        aBas = self._read_array_float(path+"cartesian2direct_aBas.in."+str(case))
        eps = self._read_float(path+"cartesian2direct_eps.in."+str(case))

        out = self._read_array_float(path+"cartesian2direct_aBas.out."+str(case))

        self._compare(out,cartesian2direct(sLV,aBas,eps))

    def test_c2d9(self):
        from phenum.vector_utils import cartesian2direct
        case = 9
        path = "tests/vector_utils/"
        
        sLV = self._read_array_float(path+"cartesian2direct_sLV.in."+str(case))
        aBas = self._read_array_float(path+"cartesian2direct_aBas.in."+str(case))
        eps = self._read_float(path+"cartesian2direct_eps.in."+str(case))

        out = self._read_array_float(path+"cartesian2direct_aBas.out."+str(case))

        self._compare(out,cartesian2direct(sLV,aBas,eps))

    def test_c2d10(self):
        from phenum.vector_utils import cartesian2direct
        case = 10
        path = "tests/vector_utils/"
        
        sLV = self._read_array_float(path+"cartesian2direct_sLV.in."+str(case))
        aBas = self._read_array_float(path+"cartesian2direct_aBas.in."+str(case))
        eps = self._read_float(path+"cartesian2direct_eps.in."+str(case))

        out = self._read_array_float(path+"cartesian2direct_aBas.out."+str(case))

        self._compare(out,cartesian2direct(sLV,aBas,eps))

    def test_c2d11(self):
        from phenum.vector_utils import cartesian2direct
        case = 11
        path = "tests/vector_utils/"
        
        sLV = self._read_array_float(path+"cartesian2direct_sLV.in."+str(case))
        aBas = self._read_array_float(path+"cartesian2direct_aBas.in."+str(case))
        eps = self._read_float(path+"cartesian2direct_eps.in."+str(case))

        out = self._read_array_float(path+"cartesian2direct_aBas.out."+str(case))

        self._compare(out,cartesian2direct(sLV,aBas,eps))

    def test_c2d12(self):
        from phenum.vector_utils import cartesian2direct
        case = 12
        path = "tests/vector_utils/"
        
        sLV = self._read_array_float(path+"cartesian2direct_sLV.in."+str(case))
        aBas = self._read_array_float(path+"cartesian2direct_aBas.in."+str(case))
        eps = self._read_float(path+"cartesian2direct_eps.in."+str(case))

        out = self._read_array_float(path+"cartesian2direct_aBas.out."+str(case))

        self._compare(out,cartesian2direct(sLV,aBas,eps))

    def test_c2d13(self):
        from phenum.vector_utils import cartesian2direct
        case = 13
        path = "tests/vector_utils/"
        
        sLV = self._read_array_float(path+"cartesian2direct_sLV.in."+str(case))
        aBas = self._read_array_float(path+"cartesian2direct_aBas.in."+str(case))
        eps = self._read_float(path+"cartesian2direct_eps.in."+str(case))

        out = self._read_array_float(path+"cartesian2direct_aBas.out."+str(case))

        self._compare(out,cartesian2direct(sLV,aBas,eps))

    def test_c2d14(self):
        from phenum.vector_utils import cartesian2direct
        case = 14
        path = "tests/vector_utils/"
        
        sLV = self._read_array_float(path+"cartesian2direct_sLV.in."+str(case))
        aBas = self._read_array_float(path+"cartesian2direct_aBas.in."+str(case))
        eps = self._read_float(path+"cartesian2direct_eps.in."+str(case))

        out = self._read_array_float(path+"cartesian2direct_aBas.out."+str(case))

        self._compare(out,cartesian2direct(sLV,aBas,eps))

    def test_c2d15(self):
        from phenum.vector_utils import cartesian2direct
        case = 15
        path = "tests/vector_utils/"
        
        sLV = self._read_array_float(path+"cartesian2direct_sLV.in."+str(case))
        aBas = self._read_array_float(path+"cartesian2direct_aBas.in."+str(case))
        eps = self._read_float(path+"cartesian2direct_eps.in."+str(case))

        out = self._read_array_float(path+"cartesian2direct_aBas.out."+str(case))

        self._compare(out,cartesian2direct(sLV,aBas,eps))

    def test_c2d16(self):
        from phenum.vector_utils import cartesian2direct
        case = 16
        path = "tests/vector_utils/"
        
        sLV = self._read_array_float(path+"cartesian2direct_sLV.in."+str(case))
        aBas = self._read_array_float(path+"cartesian2direct_aBas.in."+str(case))
        eps = self._read_float(path+"cartesian2direct_eps.in."+str(case))

        out = self._read_array_float(path+"cartesian2direct_aBas.out."+str(case))

        self._compare(out,cartesian2direct(sLV,aBas,eps))

    def test_c2d17(self):
        from phenum.vector_utils import cartesian2direct
        case = 17
        path = "tests/vector_utils/"
        
        sLV = self._read_array_float(path+"cartesian2direct_sLV.in."+str(case))
        aBas = self._read_array_float(path+"cartesian2direct_aBas.in."+str(case))
        eps = self._read_float(path+"cartesian2direct_eps.in."+str(case))

        out = self._read_array_float(path+"cartesian2direct_aBas.out."+str(case))

        self._compare(out,cartesian2direct(sLV,aBas,eps))

    def test_c2d18(self):
        from phenum.vector_utils import cartesian2direct
        case = 18
        path = "tests/vector_utils/"
        
        sLV = self._read_array_float(path+"cartesian2direct_sLV.in."+str(case))
        aBas = self._read_array_float(path+"cartesian2direct_aBas.in."+str(case))
        eps = self._read_float(path+"cartesian2direct_eps.in."+str(case))

        out = self._read_array_float(path+"cartesian2direct_aBas.out."+str(case))

        self._compare(out,cartesian2direct(sLV,aBas,eps))

    def test_c2d19(self):
        from phenum.vector_utils import cartesian2direct
        case = 19
        path = "tests/vector_utils/"
        
        sLV = self._read_array_float(path+"cartesian2direct_sLV.in."+str(case))
        aBas = self._read_array_float(path+"cartesian2direct_aBas.in."+str(case))
        eps = self._read_float(path+"cartesian2direct_eps.in."+str(case))

        out = self._read_array_float(path+"cartesian2direct_aBas.out."+str(case))

        self._compare(out,cartesian2direct(sLV,aBas,eps))

    def test_c2d20(self):
        from phenum.vector_utils import cartesian2direct
        case = 20
        path = "tests/vector_utils/"
        
        sLV = self._read_array_float(path+"cartesian2direct_sLV.in."+str(case))
        aBas = self._read_array_float(path+"cartesian2direct_aBas.in."+str(case))
        eps = self._read_float(path+"cartesian2direct_eps.in."+str(case))

        out = self._read_array_float(path+"cartesian2direct_aBas.out."+str(case))

        self._compare(out,cartesian2direct(sLV,aBas,eps))
