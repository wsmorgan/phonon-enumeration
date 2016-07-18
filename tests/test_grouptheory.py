# """Methods for testing the subroutines in the grouptheory module."""
# import unittest as ut

# def _read_output(test):
#     values = []
#     with open("tests/grouptheory/"+test) as f:
#         for line in f:
#             values.append(eval(line))
#     return values

# class TestGetSymGroup(ut.TestCase):
#     """ Tests of the get_sym_group subroutine."""

#     def test1(self):
#         from phenum.grouptheory import get_sym_group
#         par_lat = [[0.5, 0.5, 0.0], [0.5, 0.0, 0.5], [0.0, 0.5, 0.5]]
#         bas_vacs = [[0.0, 0.0, 0.0]]
#         HNF = [[1, 0, 0], [0, 1, 0], [2, 3, 6]]
#         LatDim = 3
#         out = _read_output("arrow_group.out.1")
#         symm = get_sym_group(par_lat,bas_vacs,HNF,LatDim)
#         agroup = []
#         for i in range(len(symm.perm.site_perm)):
#             agroup.append([symm.perm.site_perm[i],symm.perm.arrow_perm[i]])            
#         self.assertEqual(agroup,out)

#     def test2(self):
#         from phenum.grouptheory import get_sym_group
#         par_lat = [[0.5, 0.5, 0.0], [0.5, 0.0, 0.5], [0.0, 0.5, 0.5]]
#         bas_vacs = [[0.0, 0.0, 0.0]]
#         HNF = [[1, 0, 0], [0, 1, 0], [0, 5, 6]]
#         LatDim = 3
#         out = _read_output("arrow_group.out.2")
#         symm = get_sym_group(par_lat,bas_vacs,HNF,LatDim)
#         agroup = []
#         for i in range(len(symm.perm.site_perm)):
#             agroup.append([symm.perm.site_perm[i],symm.perm.arrow_perm[i]])            
#         self.assertEqual(agroup,out)

#     def test3(self):
#         from phenum.grouptheory import get_sym_group
#         par_lat = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
#         bas_vacs = [[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]]
#         HNF = [[1, 0, 0], [0, 1, 0], [0, 1, 2]]
#         LatDim = 3
#         out = _read_output("arrow_group.out.3")
#         symm = get_sym_group(par_lat,bas_vacs,HNF,LatDim)
#         agroup = []
#         for i in range(len(symm.perm.site_perm)):
#             agroup.append([symm.perm.site_perm[i],symm.perm.arrow_perm[i]])            
#         self.assertEqual(agroup,out)

#     def test4(self):
#         from phenum.grouptheory import get_sym_group
#         par_lat = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
#         bas_vacs = [[0.0, 0.0, 0.0]]
#         HNF = [[1, 0, 0], [0, 1, 0], [0, 0, 7]]
#         LatDim = 3
#         out = _read_output("arrow_group.out.4")
#         symm = get_sym_group(par_lat,bas_vacs,HNF,LatDim)
#         agroup = []
#         for i in range(len(symm.perm.site_perm)):
#             agroup.append([symm.perm.site_perm[i],symm.perm.arrow_perm[i]])            
#         self.assertEqual(agroup,out)

#     def test5(self):
#         from phenum.grouptheory import get_sym_group
#         par_lat = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
#         bas_vacs = [[0.0, 0.0, 0.0]]
#         HNF = [[1, 0, 0], [1, 2, 0], [1, 0, 2]]
#         LatDim = 3
#         out = _read_output("arrow_group.out.5")
#         symm = get_sym_group(par_lat,bas_vacs,HNF,LatDim)
#         agroup = []
#         for i in range(len(symm.perm.site_perm)):
#             agroup.append([symm.perm.site_perm[i],symm.perm.arrow_perm[i]])            
#         self.assertEqual(agroup,out)

#     def test6(self):
#         from phenum.grouptheory import get_sym_group
#         par_lat = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
#         bas_vacs = [[0.0, 0.0, 0.0]]
#         HNF = [[1, 0, 0], [0, 2, 0], [0, 0, 2]]
#         LatDim = 3
#         out = _read_output("arrow_group.out.6")
#         symm = get_sym_group(par_lat,bas_vacs,HNF,LatDim)
#         agroup = []
#         for i in range(len(symm.perm.site_perm)):
#             agroup.append([symm.perm.site_perm[i],symm.perm.arrow_perm[i]])            
#         self.assertEqual(agroup,out)

#     def test7(self):
#         from phenum.grouptheory import get_sym_group
#         par_lat = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
#         bas_vacs = [[0.0, 0.0, 0.0], [0.5, 0.5, 0.5], [0.25, 0.25, 0.75]]
#         HNF = [[1, 0, 0], [0, 1, 0], [0, 1, 2]]
#         LatDim = 3
#         out = _read_output("arrow_group.out.7")
#         symm = get_sym_group(par_lat,bas_vacs,HNF,LatDim)
#         agroup = []
#         for i in range(len(symm.perm.site_perm)):
#             agroup.append([symm.perm.site_perm[i],symm.perm.arrow_perm[i]])            
#         self.assertEqual(agroup,out)

#     def test8(self):
#         from phenum.grouptheory import get_sym_group
#         par_lat = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
#         bas_vacs = [[0.0, 0.0, 0.0], [0.5, 0.5, 0.5], [0.25, 0.25, 0.75]]
#         HNF = [[1, 0, 0], [0, 1, 0], [0, 0, 3]]
#         LatDim = 3
#         out = _read_output("arrow_group.out.8")
#         symm = get_sym_group(par_lat,bas_vacs,HNF,LatDim)
#         agroup = []
#         for i in range(len(symm.perm.site_perm)):
#             agroup.append([symm.perm.site_perm[i],symm.perm.arrow_perm[i]])            
#         self.assertEqual(agroup,out)

#     def test9(self):
#         from phenum.grouptheory import get_sym_group
#         par_lat = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
#         bas_vacs = [[0.0, 0.0, 0.0], [0.5, 0.5, 0.5], [0.25, 0.25, 0.75]]
#         HNF = [[1, 0, 0], [0, 1, 0], [0, 2, 3]]
#         LatDim = 3
#         out = _read_output("arrow_group.out.9")
#         symm = get_sym_group(par_lat,bas_vacs,HNF,LatDim)
#         agroup = []
#         for i in range(len(symm.perm.site_perm)):
#             agroup.append([symm.perm.site_perm[i],symm.perm.arrow_perm[i]])            
#         self.assertEqual(agroup,out)

#     def test10(self):
#         from phenum.grouptheory import get_sym_group
#         par_lat = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
#         bas_vacs = [[0.0, 0.0, 0.0]]
#         HNF = [[2, 0, 0], [0, 2, 0], [0, 0, 2]]
#         LatDim = 3
#         out = _read_output("arrow_group.out.10")
#         symm = get_sym_group(par_lat,bas_vacs,HNF,LatDim)
#         agroup = []
#         for i in range(len(symm.perm.site_perm)):
#             agroup.append([symm.perm.site_perm[i],symm.perm.arrow_perm[i]])            
#         self.assertEqual(agroup,out)

# class TestGetFullHNF(ut.TestCase):
#     """ Tests of the get_full_HNF subroutine."""

#     def test1(self):
#         from phenum.grouptheory import get_full_HNF
#         from numpy import array
#         HNF = array([1,0,1,0,0,1])
#         out = [[1,0,0],[0,1,0],[0,0,1]]
#         self.assertEqual(get_full_HNF(HNF),out)

#     def test2(self):
#         from phenum.grouptheory import get_full_HNF
#         from numpy import array
#         HNF = array([2,1,2,1,0,4])
#         out = [[2,0,0],[1,2,0],[1,0,4]]
#         self.assertEqual(get_full_HNF(HNF),out)

#     def test3(self):
#         from phenum.grouptheory import get_full_HNF
#         from numpy import array
#         HNF = array([1,0,3,1,2,3])
#         out = [[1,0,0],[0,3,0],[1,2,3]]
#         self.assertEqual(get_full_HNF(HNF),out)

#     def test4(self):
#         from phenum.grouptheory import get_full_HNF
#         from numpy import array
#         HNF = array([0,0,0,0,0,0])
#         out = [[0,0,0],[0,0,0],[0,0,0]]
#         self.assertEqual(get_full_HNF(HNF),out)

#     def test5(self):
#         from phenum.grouptheory import get_full_HNF
#         from numpy import array
#         HNF = array([3,0,3,0,0,3])
#         out = [[3,0,0],[0,3,0],[0,0,3]]
#         self.assertEqual(get_full_HNF(HNF),out)

#     def test1(self):
#         from phenum.grouptheory import get_full_HNF
#         from numpy import array
#         HNF = array([1,1,2,0,2,2])
#         out = [[1,0,0],[1,2,0],[0,2,2]]
#         self.assertEqual(get_full_HNF(HNF),out)

#     def test1(self):
#         from phenum.grouptheory import get_full_HNF
#         from numpy import array
#         HNF = array([2,0,2,0,2,4])
#         out = [[2,0,0],[0,2,0],[0,2,4]]
#         self.assertEqual(get_full_HNF(HNF),out)

# class TestSmithNormalForm(ut.TestCase):
#     """ Tests of the SmithNormalForm subroutine."""

#     def test1(self):
#         from phenum.grouptheory import SmithNormalForm
#         HNF =  [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
#         out =  ([[1, 0, 0], [0, 1, 0], [0, 0, 1]], [[1, 0, 0], [0, 1, 0], [0, 0, 1]], [[1, 0, 0], [0, 1, 0], [0, 0, 1]])
#         self.assertEqual(SmithNormalForm(HNF),out)

#     def test2(self):
#         from phenum.grouptheory import SmithNormalForm
#         HNF =  [[1, 0, 0], [0, 1, 0], [0, 1, 2]]
#         out =  ([[1, 0, 0], [0, 1, 0], [0, 0, 2]], [[1, 0, 0], [0, 1, 0], [0, -1, 1]], [[1, 0, 0], [0, 1, 0], [0, 0, 1]])
#         self.assertEqual(SmithNormalForm(HNF),out)

#     def test3(self):
#         from phenum.grouptheory import SmithNormalForm
#         HNF =  [[1, 0, 0], [0, 1, 0], [0, 0, 3]]
#         out =  ([[1, 0, 0], [0, 1, 0], [0, 0, 3]], [[1, 0, 0], [0, 1, 0], [0, 0, 1]], [[1, 0, 0], [0, 1, 0], [0, 0, 1]])
#         self.assertEqual(SmithNormalForm(HNF),out)

#     def test4(self):
#         from phenum.grouptheory import SmithNormalForm
#         HNF =  [[1, 0, 0], [0, 2, 0], [0, 0, 2]]
#         out =  ([[1, 0, 0], [0, 2, 0], [0, 0, 2]], [[1, 0, 0], [0, 1, 0], [0, 0, 1]], [[1, 0, 0], [0, 1, 0], [0, 0, 1]])
#         self.assertEqual(SmithNormalForm(HNF),out)

#     def test5(self):
#         from phenum.grouptheory import SmithNormalForm
#         HNF =  [[1, 0, 0], [0, 1, 0], [1, 2, 5]]
#         out =  ([[1, 0, 0], [0, 1, 0], [0, 0, 5]], [[1, 0, 0], [0, 1, 0], [-1, -2, 1]], [[1, 0, 0], [0, 1, 0], [0, 0, 1]])
#         self.assertEqual(SmithNormalForm(HNF),out)

#     def test6(self):
#         from phenum.grouptheory import SmithNormalForm
#         HNF =  [[1, 0, 0], [0, 1, 0], [2, 3, 6]]
#         out =  ([[1, 0, 0], [0, 1, 0], [0, 0, 6]], [[1, 0, 0], [0, 1, 0], [-2, -3, 1]], [[1, 0, 0], [0, 1, 0], [0, 0, 1]])
#         self.assertEqual(SmithNormalForm(HNF),out)

#     def test7(self):
#         from phenum.grouptheory import SmithNormalForm
#         HNF =  [[1, 0, 0], [0, 1, 0], [0, 6, 7]]
#         out =  ([[1, 0, 0], [0, 1, 0], [0, 0, 7]], [[1, 0, 0], [0, 1, 0], [0, -6, 1]], [[1, 0, 0], [0, 1, 0], [0, 0, 1]])
#         self.assertEqual(SmithNormalForm(HNF),out)

#     def test8(self):
#         from phenum.grouptheory import SmithNormalForm
#         HNF =  [[1, 0, 0], [1, 2, 0], [1, 0, 4]]
#         out =  ([[1, 0, 0], [0, 2, 0], [0, 0, 4]], [[1, 0, 0], [-1, 1, 0], [-1, 0, 1]], [[1, 0, 0], [0, 1, 0], [0, 0, 1]])
#         self.assertEqual(SmithNormalForm(HNF),out)

#     def test9(self):
#         from phenum.grouptheory import SmithNormalForm
#         HNF =  [[2, 0, 0], [0, 2, 0], [0, 0, 2]]
#         out =  ([[2, 0, 0], [0, 2, 0], [0, 0, 2]], [[1, 0, 0], [0, 1, 0], [0, 0, 1]], [[1, 0, 0], [0, 1, 0], [0, 0, 1]])
#         self.assertEqual(SmithNormalForm(HNF),out)

#     def test10(self):
#         from phenum.grouptheory import SmithNormalForm
#         HNF =  [[1, 0, 0], [0, 1, 0], [1, 5, 10]]
#         out =  ([[1, 0, 0], [0, 1, 0], [0, 0, 10]], [[1, 0, 0], [0, 1, 0], [-1, -5, 1]], [[1, 0, 0], [0, 1, 0], [0, 0, 1]])
#         self.assertEqual(SmithNormalForm(HNF),out)

# class TestAGroup(ut.TestCase):
#     """ Tests of the a_group subroutine."""

#     def test1(self):
#         from phenum.grouptheory import a_group
#         trans = [[0,1],[1,0]]
#         rots = [[[0,1],[0,1,2,3,4,5]],[[1,0],[2,3,0,1,5,4]],[[1,0],[2,1,0,3,5,4]],[[0,1],[0,3,2,1,5,4]]]
#         out = _read_output("agroup.out.1")
#         self.assertEqual(a_group(trans,rots),out)

#     def test2(self):
#         from phenum.grouptheory import a_group
#         trans = [[j-1 for j in i] for i in [[1, 2, 3, 4], [2, 1, 4, 3], [3, 4, 1, 2], [4, 3, 2, 1]]]
#         rots = [[[j-1 for j in i] for i in t] for t in [[[1, 2, 3, 4], [1, 2, 3, 4, 5, 6]], [[1, 4, 3, 2], [1, 3, 2, 4, 6, 5]], [[1, 2, 3, 4], [4, 2, 3, 1, 5, 6]], [[1, 4, 3, 2], [4, 3, 2, 1, 6, 5]], [[1, 2, 3, 4], [1, 5, 3, 4, 2, 6]], [[1, 4, 3, 2], [1, 3, 5, 4, 6, 2]], [[1, 2, 3, 4], [4, 5, 3, 1, 2, 6]], [[1, 4, 3, 2], [4, 3, 5, 1, 6, 2]], [[1, 2, 3, 4], [1, 2, 6, 4, 5, 3]], [[1, 4, 3, 2], [1, 6, 2, 4, 3, 5]], [[1, 2, 3, 4], [4, 2, 6, 1, 5, 3]], [[1, 4, 3, 2], [4, 6, 2, 1, 3, 5]], [[1, 2, 3, 4], [1, 5, 6, 4, 2, 3]], [[1, 4, 3, 2], [1, 6, 5, 4, 3, 2]], [[1, 2, 3, 4], [4, 5, 6, 1, 2, 3]], [[1, 4, 3, 2], [4, 6, 5, 1, 3, 2]]]]
#         out = _read_output("agroup.out.2")
#         self.assertEqual(a_group(trans,rots),out)

#     def test3(self):
#         from phenum.grouptheory import a_group
#         trans = [[j-1 for j in i] for i in [[1, 2, 3, 4, 5, 6, 7, 8], [2, 1, 4, 3, 6, 5, 8, 7], [3, 4, 5, 6, 7, 8, 1, 2], [4, 3, 6, 5, 8, 7, 2, 1], [5, 6, 7, 8, 1, 2, 3, 4], [6, 5, 8, 7, 2, 1, 4, 3], [7, 8, 1, 2, 3, 4, 5, 6], [8, 7, 2, 1, 4, 3, 6, 5]]]
#         rots = [[[j-1 for j in i] for i in t] for t in [[[1, 2, 3, 4, 5, 6, 7, 8], [1, 2, 3, 4, 5, 6]], [[1, 2, 3, 4, 5, 6, 7, 8], [4, 2, 3, 1, 5, 6]], [[1, 2, 3, 4, 5, 6, 7, 8], [1, 5, 3, 4, 2, 6]], [[1, 2, 3, 4, 5, 6, 7, 8], [4, 5, 3, 1, 2, 6]], [[1, 2, 7, 8, 5, 6, 3, 4], [1, 2, 6, 4, 5, 3]], [[1, 2, 7, 8, 5, 6, 3, 4], [4, 2, 6, 1, 5, 3]], [[1, 2, 7, 8, 5, 6, 3, 4], [1, 5, 6, 4, 2, 3]], [[1, 2, 7, 8, 5, 6, 3, 4], [4, 5, 6, 1, 2, 3]]]]
#         out = _read_output("agroup.out.3")
#         self.assertEqual(a_group(trans,rots),out)

#     def test4(self):
#         from phenum.grouptheory import a_group
#         trans =[[j - 1 for j in i] for i in[[1,2,3,4,5,6,7,8], [2,1,4,3,6,5,8,7], [3,4,5,6,7,8,1,2], [4,3,6,5,8,7,2,1], [5,6,7,8,1,2,3,4], [6,5,8,7,2,1,4,3], [7,8,1,2,3,4,5,6], [8,7,2,1,4,3,6,5]]]
#         rots = [[[0,1,2,3,4,5,6,7],[0,1,2,3]],[[0,1,2,3,4,5,6,7],[2,1,0,3]],[[0,1,6,7,4,5,2,3],[0,3,2,1]],[[0,1,6,7,4,5,2,3],[2,3,0,1]]]
#         out = _read_output("agroup.out.4")
#         self.assertEqual(a_group(trans,rots),out)

#     def test5(self):
#         from phenum.grouptheory import a_group
#         trans =[[j - 1 for j in i] for i in[[1, 2, 3, 4, 5, 6, 7, 8], [2, 1, 4, 3, 6, 5, 8, 7], [3, 4, 5, 6, 7, 8, 1, 2], [4, 3, 6, 5, 8, 7, 2, 1], [5, 6, 7, 8, 1, 2, 3, 4], [6, 5, 8, 7, 2, 1, 4, 3], [7, 8, 1, 2, 3, 4, 5, 6], [8, 7, 2, 1, 4, 3, 6, 5]]]
#         rots = [[[0,1,2,3,4,5,6,7],[0,1,2,3]],[[0,1,2,3,4,5,6,7],[2,1,0,3]],[[0,1,6,7,4,5,2,3],[0,3,2,1]],[[0,1,6,7,4,5,2,3],[2,3,0,1]]]
#         out = _read_output("agroup.out.5")
#         self.assertEqual(a_group(trans,rots),out)

# class TestAGroupGen(ut.TestCase):
#     """ Tests of the a_group subroutine."""

#     def test1(self):
#         from phenum.grouptheory import a_group_gen
#         trans = [[0,1],[1,0]]
#         rots = [[[0,1],[0,1,2,3,4,5]],[[1,0],[2,3,0,1,5,4]],[[1,0],[2,1,0,3,5,4]],[[0,1],[0,3,2,1,5,4]]]
#         out = _read_output("agroupgen.out.1")
#         self.assertEqual(a_group_gen(trans,rots),out)

#     def test2(self):
#         from phenum.grouptheory import a_group_gen
#         trans = [[j-1 for j in i] for i in [[1, 2, 3, 4], [2, 1, 4, 3], [3, 4, 1, 2], [4, 3, 2, 1]]]
#         rots = [[[j-1 for j in i] for i in t] for t in [[[1, 2, 3, 4], [1, 2, 3, 4, 5, 6]], [[1, 4, 3, 2], [1, 3, 2, 4, 6, 5]], [[1, 2, 3, 4], [4, 2, 3, 1, 5, 6]], [[1, 4, 3, 2], [4, 3, 2, 1, 6, 5]], [[1, 2, 3, 4], [1, 5, 3, 4, 2, 6]], [[1, 4, 3, 2], [1, 3, 5, 4, 6, 2]], [[1, 2, 3, 4], [4, 5, 3, 1, 2, 6]], [[1, 4, 3, 2], [4, 3, 5, 1, 6, 2]], [[1, 2, 3, 4], [1, 2, 6, 4, 5, 3]], [[1, 4, 3, 2], [1, 6, 2, 4, 3, 5]], [[1, 2, 3, 4], [4, 2, 6, 1, 5, 3]], [[1, 4, 3, 2], [4, 6, 2, 1, 3, 5]], [[1, 2, 3, 4], [1, 5, 6, 4, 2, 3]], [[1, 4, 3, 2], [1, 6, 5, 4, 3, 2]], [[1, 2, 3, 4], [4, 5, 6, 1, 2, 3]], [[1, 4, 3, 2], [4, 6, 5, 1, 3, 2]]]]
#         out = _read_output("agroupgen.out.2")
#         self.assertEqual(a_group_gen(trans,rots),out)

#     def test3(self):
#         from phenum.grouptheory import a_group_gen
#         trans = [[j-1 for j in i] for i in [[1, 2, 3, 4, 5, 6, 7, 8], [2, 1, 4, 3, 6, 5, 8, 7], [3, 4, 5, 6, 7, 8, 1, 2], [4, 3, 6, 5, 8, 7, 2, 1], [5, 6, 7, 8, 1, 2, 3, 4], [6, 5, 8, 7, 2, 1, 4, 3], [7, 8, 1, 2, 3, 4, 5, 6], [8, 7, 2, 1, 4, 3, 6, 5]]]
#         rots = [[[j-1 for j in i] for i in t] for t in [[[1, 2, 3, 4, 5, 6, 7, 8], [1, 2, 3, 4, 5, 6]], [[1, 2, 3, 4, 5, 6, 7, 8], [4, 2, 3, 1, 5, 6]], [[1, 2, 3, 4, 5, 6, 7, 8], [1, 5, 3, 4, 2, 6]], [[1, 2, 3, 4, 5, 6, 7, 8], [4, 5, 3, 1, 2, 6]], [[1, 2, 7, 8, 5, 6, 3, 4], [1, 2, 6, 4, 5, 3]], [[1, 2, 7, 8, 5, 6, 3, 4], [4, 2, 6, 1, 5, 3]], [[1, 2, 7, 8, 5, 6, 3, 4], [1, 5, 6, 4, 2, 3]], [[1, 2, 7, 8, 5, 6, 3, 4], [4, 5, 6, 1, 2, 3]]]]
#         out = _read_output("agroupgen.out.3")
#         self.assertEqual(a_group_gen(trans,rots),out)

#     def test4(self):
#         from phenum.grouptheory import a_group_gen
#         trans =[[j - 1 for j in i] for i in[[1,2,3,4,5,6,7,8], [2,1,4,3,6,5,8,7], [3,4,5,6,7,8,1,2], [4,3,6,5,8,7,2,1], [5,6,7,8,1,2,3,4], [6,5,8,7,2,1,4,3], [7,8,1,2,3,4,5,6], [8,7,2,1,4,3,6,5]]]
#         rots = [[[0,1,2,3,4,5,6,7],[0,1,2,3]],[[0,1,2,3,4,5,6,7],[2,1,0,3]],[[0,1,6,7,4,5,2,3],[0,3,2,1]],[[0,1,6,7,4,5,2,3],[2,3,0,1]]]
#         out = _read_output("agroupgen.out.4")
#         self.assertEqual(a_group_gen(trans,rots),out)

#     def test5(self):
#         from phenum.grouptheory import a_group_gen
#         trans =[[j - 1 for j in i] for i in[[1, 2, 3, 4, 5, 6, 7, 8], [2, 1, 4, 3, 6, 5, 8, 7], [3, 4, 5, 6, 7, 8, 1, 2], [4, 3, 6, 5, 8, 7, 2, 1], [5, 6, 7, 8, 1, 2, 3, 4], [6, 5, 8, 7, 2, 1, 4, 3], [7, 8, 1, 2, 3, 4, 5, 6], [8, 7, 2, 1, 4, 3, 6, 5]]]
#         rots = [[[0,1,2,3,4,5,6,7],[0,1,2,3]],[[0,1,2,3,4,5,6,7],[2,1,0,3]],[[0,1,6,7,4,5,2,3],[0,3,2,1]],[[0,1,6,7,4,5,2,3],[2,3,0,1]]]
#         out = _read_output("agroupgen.out.5")
#         self.assertEqual(a_group_gen(trans,rots),out)
