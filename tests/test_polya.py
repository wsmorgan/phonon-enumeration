"""Methods for testing the subroutines in the polyaburnside module."""
import unittest as ut
from phenum.polyaburnside import polya
import pytest

class Testpolya(ut.TestCase):
    """Tests of the polya subroutine."""

    def test_toy_2(self):
        from phenum.grouptheory import a_group_gen
        from phenum.phonons import how_many_arrows, how_many_arrows
        col = [[-1,1],[1,2]]
        trans = [[0,1],[1,0]]
        rots = [[[0,1],[0,1,2,3,4,5]],[[1,0],[2,3,0,1,5,4]],[[1,0],[2,1,0,3,5,4]],
                [[0,1],[0,3,2,1,5,4]]]
        dim = 6
        (narrows,arrow_types,Concs) = how_many_arrows(col)
        agroup = a_group_gen(trans,rots)
        out = 3
        self.assertEqual(polya(Concs,agroup,arrowings=arrow_types),out)
        
    def test_toy_3(self):
        from phenum.grouptheory import a_group
        from phenum.phonons import how_many_arrows, how_many_arrows
        col = [[-1,1],[1,2]]
        trans = [[0,1],[1,0]]
        rots = [[[0,1],[0,1,2,3,4,5]],[[1,0],[2,3,0,1,5,4]],[[1,0],[2,1,0,3,5,4]],
                [[0,1],[0,3,2,1,5,4]],[[0,1],[0,4,2,5,3,1]],[[0,1],[0,5,2,4,1,3]]]
        dim = 6
        (narrows,arrow_types,Concs) = how_many_arrows(col)
        agroup = a_group(trans,rots)
        out = 2
        self.assertEqual(polya(Concs,agroup,arrowings=arrow_types),out)
        
    def test_toy_1(self):
        from phenum.grouptheory import a_group
        from phenum.phonons import how_many_arrows, how_many_arrows
        col = [[-1,1],[-1,2]]
        trans = [[0,1],[1,0]]
        rots = [[[0,1],[0,1,2,3,4,5]],[[1,0],[2,3,0,1,5,4]],[[1,0],[2,1,0,3,5,4]],
                [[0,1],[0,3,2,1,5,4]]]
        dim = 6
        (narrows,arrow_types,Concs) = how_many_arrows(col)
        agroup = a_group(trans,rots)
        out = 1
        self.assertEqual(polya(Concs,agroup,arrowings=arrow_types),out)

    def test_3g_1(self):
        from phenum.grouptheory import a_group
        from phenum.phonons import how_many_arrows, how_many_arrows
        dim = 6
        col = [[1, 2], [1, 3], [1, 4], [1, 4]]
        trans = [[j-1 for j in i] for i in [[1, 2, 3, 4], [2, 1, 4, 3], [3, 4, 1, 2], [4, 3, 2, 1]]]
        rots = [[[j-1 for j in i] for i in t] for t in [[[1, 2, 3, 4], [1, 2, 3, 4, 5, 6]], [[1, 4, 3, 2], [1, 3, 2, 4, 6, 5]], [[1, 2, 3, 4], [4, 2, 3, 1, 5, 6]], [[1, 4, 3, 2], [4, 3, 2, 1, 6, 5]], [[1, 2, 3, 4], [1, 5, 3, 4, 2, 6]], [[1, 4, 3, 2], [1, 3, 5, 4, 6, 2]], [[1, 2, 3, 4], [4, 5, 3, 1, 2, 6]], [[1, 4, 3, 2], [4, 3, 5, 1, 6, 2]], [[1, 2, 3, 4], [1, 2, 6, 4, 5, 3]], [[1, 4, 3, 2], [1, 6, 2, 4, 3, 5]], [[1, 2, 3, 4], [4, 2, 6, 1, 5, 3]], [[1, 4, 3, 2], [4, 6, 2, 1, 3, 5]], [[1, 2, 3, 4], [1, 5, 6, 4, 2, 3]], [[1, 4, 3, 2], [1, 6, 5, 4, 3, 2]], [[1, 2, 3, 4], [4, 5, 6, 1, 2, 3]], [[1, 4, 3, 2], [4, 6, 5, 1, 3, 2]]]]
        (narrows,arrow_types,Concs) = how_many_arrows(col)
        agroup = a_group(trans,rots)
        out = 400
        self.assertEqual(polya(Concs,agroup,arrowings=arrow_types),out)

    def test_3g_3(self):
        from phenum.grouptheory import a_group
        from phenum.phonons import how_many_arrows, how_many_arrows
        dim = 6
        col = [[1,1],[1,2],[1,3],[1,4]]
        trans = [[j-1 for j in i] for i in [[1, 2, 3, 4], [2, 1, 4, 3], [3, 4, 1, 2], [4, 3, 2, 1]]]
        rots = [[[j-1 for j in i] for i in t] for t in [[[1, 2, 3, 4], [1, 2, 3, 4, 5, 6]], [[1, 4, 3, 2], [1, 3, 2, 4, 6, 5]], [[1, 2, 3, 4], [4, 2, 3, 1, 5, 6]], [[1, 4, 3, 2], [4, 3, 2, 1, 6, 5]], [[1, 2, 3, 4], [1, 5, 3, 4, 2, 6]], [[1, 4, 3, 2], [1, 3, 5, 4, 6, 2]], [[1, 2, 3, 4], [4, 5, 3, 1, 2, 6]], [[1, 4, 3, 2], [4, 3, 5, 1, 6, 2]], [[1, 2, 3, 4], [1, 2, 6, 4, 5, 3]], [[1, 4, 3, 2], [1, 6, 2, 4, 3, 5]], [[1, 2, 3, 4], [4, 2, 6, 1, 5, 3]], [[1, 4, 3, 2], [4, 6, 2, 1, 3, 5]], [[1, 2, 3, 4], [1, 5, 6, 4, 2, 3]], [[1, 4, 3, 2], [1, 6, 5, 4, 3, 2]], [[1, 2, 3, 4], [4, 5, 6, 1, 2, 3]], [[1, 4, 3, 2], [4, 6, 5, 1, 3, 2]]]]
        (narrows,arrow_types,Concs) = how_many_arrows(col)
        agroup = a_group(trans,rots)
        out = 792
        self.assertEqual(polya(Concs,agroup,arrowings=arrow_types),out)

    def test_3g_2(self):
        from phenum.grouptheory import a_group
        from phenum.phonons import how_many_arrows, how_many_arrows
        dim = 6
        col = [[1,1],[1,1],[1,1],[1,1]]
        trans = [[j-1 for j in i] for i in [[1, 2, 3, 4], [2, 1, 4, 3], [3, 4, 1, 2], [4, 3, 2, 1]]]
        rots = [[[j-1 for j in i] for i in t] for t in [[[1, 2, 3, 4], [1, 2, 3, 4, 5, 6]], [[1, 4, 3, 2], [1, 3, 2, 4, 6, 5]], [[1, 2, 3, 4], [4, 2, 3, 1, 5, 6]], [[1, 4, 3, 2], [4, 3, 2, 1, 6, 5]], [[1, 2, 3, 4], [1, 5, 3, 4, 2, 6]], [[1, 4, 3, 2], [1, 3, 5, 4, 6, 2]], [[1, 2, 3, 4], [4, 5, 3, 1, 2, 6]], [[1, 4, 3, 2], [4, 3, 5, 1, 6, 2]], [[1, 2, 3, 4], [1, 2, 6, 4, 5, 3]], [[1, 4, 3, 2], [1, 6, 2, 4, 3, 5]], [[1, 2, 3, 4], [4, 2, 6, 1, 5, 3]], [[1, 4, 3, 2], [4, 6, 2, 1, 3, 5]], [[1, 2, 3, 4], [1, 5, 6, 4, 2, 3]], [[1, 4, 3, 2], [1, 6, 5, 4, 3, 2]], [[1, 2, 3, 4], [4, 5, 6, 1, 2, 3]], [[1, 4, 3, 2], [4, 6, 5, 1, 3, 2]]]]
        (narrows,arrow_types,Concs) = how_many_arrows(col)
        agroup = a_group(trans,rots)
        out = 50
        self.assertEqual(polya(Concs,agroup,arrowings=arrow_types),out)

    def test_3f_1(self):
        from phenum.grouptheory import a_group
        from phenum.phonons import how_many_arrows, how_many_arrows
        dim = 6
        col = [[-1, 1], [1, 2], [1, 2], [1, 3], [1, 3], [1, 3]]
        trans = [[j-1 for j in i] for i in [[1, 2, 3, 4, 5, 6], [2, 3, 4, 5, 6, 1], [3, 4, 5, 6, 1, 2], [4, 5, 6, 1, 2, 3], [5, 6, 1, 2, 3, 4], [6, 1, 2, 3, 4, 5]]]
        rots = [[[j-1 for j in i] for i in t] for t in [[[1, 2, 3, 4, 5, 6], [1, 2, 3, 4, 5, 6]], [[1, 2, 3, 4, 5, 6], [4, 2, 3, 1, 5, 6]], [[1, 2, 3, 4, 5, 6], [1, 5, 3, 4, 2, 6]], [[1, 2, 3, 4, 5, 6], [4, 5, 3, 1, 2, 6]], [[1, 6, 5, 4, 3, 2], [1, 2, 6, 4, 5, 3]], [[1, 6, 5, 4, 3, 2], [4, 2, 6, 1, 5, 3]], [[1, 6, 5, 4, 3, 2], [1, 5, 6, 4, 2, 3]], [[1, 6, 5, 4, 3, 2], [4, 5, 6, 1, 2, 3]]]]
        (narrows,arrow_types,Concs) = how_many_arrows(col)
        agroup = a_group(trans,rots)
        out = 12392
        self.assertEqual(polya(Concs,agroup,arrowings=arrow_types),out)

    def test_3e_1(self):
        from phenum.grouptheory import a_group
        from phenum.phonons import how_many_arrows, how_many_arrows
        dim = 6
        col = [[-1, 1], [-1, 2], [1, 3], [1, 3], [1, 4], [1, 4]]
        trans = [[j-1 for j in i] for i in [[1, 2, 3, 4, 5, 6], [2, 3, 4, 5, 6, 1], [3, 4, 5, 6, 1, 2], [4, 5, 6, 1, 2, 3], [5, 6, 1, 2, 3, 4], [6, 1, 2, 3, 4, 5]]]
        rots = [[[j-1 for j in i] for i in t] for t in [[[1, 2, 3, 4, 5, 6], [1, 2, 3, 4, 5, 6]], [[1, 2, 3, 4, 5, 6], [4, 2, 3, 1, 5, 6]], [[1, 2, 3, 4, 5, 6], [1, 5, 3, 4, 2, 6]], [[1, 2, 3, 4, 5, 6], [4, 5, 3, 1, 2, 6]], [[1, 6, 5, 4, 3, 2], [1, 2, 6, 4, 5, 3]], [[1, 6, 5, 4, 3, 2], [4, 2, 6, 1, 5, 3]], [[1, 6, 5, 4, 3, 2], [1, 5, 6, 4, 2, 3]], [[1, 6, 5, 4, 3, 2], [4, 5, 6, 1, 2, 3]]]]
        (narrows,arrow_types,Concs) = how_many_arrows(col)
        agroup = a_group(trans,rots)
        out = 6876
        self.assertEqual(polya(Concs,agroup,arrowings=arrow_types),out)

    def test_3h_1(self):
        from phenum.grouptheory import a_group
        from phenum.phonons import how_many_arrows, how_many_arrows
        dim = 6
        col=[[1,2],[1,2],[1,3],[1,3],[1,3],[1,3]]
        trans = [[j-1 for j in i] for i in[[1, 2, 3, 4, 5, 6],[2, 3, 4, 5, 6, 1], [3, 4, 5, 6, 1, 2], [4, 5, 6, 1, 2, 3], [5, 6, 1, 2, 3, 4], [6, 1, 2, 3, 4, 5]]]
        rots = [[[j -1 for j in i] for i in t] for t in[[[1, 2, 3, 4, 5, 6], [1, 2, 3, 4, 5, 6]], [[1, 2, 3, 4, 5, 6], [4, 2, 3, 1, 5, 6]], [[1, 2, 3, 4, 5, 6], [1, 5, 3, 4, 2, 6]], [[1, 2, 3, 4, 5, 6], [4, 5, 3, 1, 2, 6]], [[1, 6, 5, 4, 3, 2], [1, 2, 6, 4, 5, 3]], [[1, 6, 5, 4, 3, 2], [4, 2, 6, 1, 5, 3]], [[1, 6, 5, 4, 3, 2], [1, 5, 6, 4, 2, 3]], [[1, 6, 5, 4, 3, 2], [4, 5, 6, 1, 2, 3]]]]
        (narrows,arrow_types,Concs) = how_many_arrows(col)
        agroup = a_group(trans,rots)
        out = 17538
        self.assertEqual(polya(Concs,agroup,arrowings=arrow_types),out)

    def test_pg(self):
        from phenum.grouptheory import a_group
        from phenum.phonons import how_many_arrows, how_many_arrows
        dim = 6
        col = [[-1,1],[-1,1],[-1,1],[-1,2],[-1,2],[-1,2],[-1,2],[1,1],[1,1]]
        trans = [[j - 1 for j in i] for i in [[1,2,3,4,5,6,7,8,9], [3,1,2,6,4,5,9,7,8], [2,3,1,5,6,4,8,9,7],[7,8,9,1,2,3,4,5,6], [9,7,8,3,1,2,6,4,5], [8,9,7,2,3,1,5,6,4], [4,5,6,7,8,9,1,2,3], [6,4,5,9,7,8,3,1,2], [5,6,4,8,9,7,2,3,1]]]
        rots = [[[j - 1 for j in i] for i in t] for t in [[[1,2,3,4,5,6,7,8,9],[1,2,3,4,5,6]], [[3,6,9,2,5,8,1,4,7],[3,4,2,1,5,6]], [[9,8,7,6,5,4,3,2,1],[2,1,4,3,5,6]], [[7,4,1,8,5,2,9,6,3],[4,3,1,2,5,6]], [[1,4,7,2,5,8,3,6,9],[3,4,1,2,6,5]], [[9,6,3,8,5,2,7,4,1],[4,3,2,1,6,5]], [[7,8,9,4,5,6,1,2,3],[1,2,4,3,6,5]], [[3,2,1,6,5,4,9,8,7],[2,1,3,4,6,5]]]]
        (narrows,arrow_types,Concs) = how_many_arrows(col)
        agroup = a_group(trans,rots)
        out = 663
        self.assertEqual(polya(Concs,agroup,arrowings=arrow_types),out)

    def test_3d_1(self):
        from phenum.grouptheory import a_group
        from phenum.phonons import how_many_arrows, how_many_arrows
        dim = 6
        col = [[-1, 1], [-1, 1], [-1, 2], [-1, 2], [-1, 2], [0, 3], [0, 4], [0, 4]]
        trans =[[j - 1 for j in i] for i in[[1, 2, 3, 4, 5, 6, 7, 8], [2, 1, 4, 3, 6, 5, 8, 7], [3, 4, 5, 6, 7, 8, 1, 2], [4, 3, 6, 5, 8, 7, 2, 1], [5, 6, 7, 8, 1, 2, 3, 4], [6, 5, 8, 7, 2, 1, 4, 3], [7, 8, 1, 2, 3, 4, 5, 6], [8, 7, 2, 1, 4, 3, 6, 5]]]
        rots = [[[j - 1 for j in i] for i in t] for t in [[[1, 2, 3, 4, 5, 6, 7, 8], [1, 2, 3, 4, 5, 6]], [[1, 2, 3, 4, 5, 6, 7, 8], [4, 2, 3, 1, 5, 6]], [[1, 2, 3, 4, 5, 6, 7, 8], [1, 5, 3, 4, 2, 6]], [[1, 2, 3, 4, 5, 6, 7, 8], [4, 5, 3, 1, 2, 6]], [[1, 2, 7, 8, 5, 6, 3, 4], [1, 2, 6, 4, 5, 3]], [[1, 2, 7, 8, 5, 6, 3, 4], [4, 2, 6, 1, 5, 3]], [[1, 2, 7, 8, 5, 6, 3, 4], [1, 5, 6, 4, 2, 3]], [[1, 2, 7, 8, 5, 6, 3, 4], [4, 5, 6, 1, 2, 3]]]]
        (narrows,arrow_types,Concs) = how_many_arrows(col)
        agroup = a_group(trans,rots)
        out = 9348
        self.assertEqual(polya(Concs,agroup,arrowings=arrow_types),out)

    def test_p10(self):
        from phenum.grouptheory import a_group
        from phenum.phonons import how_many_arrows, how_many_arrows
        dim = 4
        col = [[-1,1],[-1,1],[-1,2],[-1,2],[0,4],[0,4],[0,5],[0,5]]
        trans =[[j - 1 for j in i] for i in[[1,2,3,4,5,6,7,8], [2,1,4,3,6,5,8,7], [3,4,5,6,7,8,1,2], [4,3,6,5,8,7,2,1], [5,6,7,8,1,2,3,4], [6,5,8,7,2,1,4,3], [7,8,1,2,3,4,5,6], [8,7,2,1,4,3,6,5]]]
        rots = [[[0,1,2,3,4,5,6,7],[0,1,2,3]],[[0,1,2,3,4,5,6,7],[2,1,0,3]],[[0,1,6,7,4,5,2,3],[0,3,2,1]],[[0,1,6,7,4,5,2,3],[2,3,0,1]]]
        (narrows,arrow_types,Concs) = how_many_arrows(col)
        agroup = a_group(trans,rots)
        out = 21720
        self.assertEqual(polya(Concs,agroup,arrowings=arrow_types),out)

    def test_p9(self):
        from phenum.grouptheory import a_group
        from phenum.phonons import how_many_arrows, how_many_arrows
        dim = 4
        col = [[-1,1],[-1,1],[-1,1],[0,1],[0,2],[0,2],[0,3],[0,3]]
        trans =[[j - 1 for j in i] for i in[[1,2,3,4,5,6,7,8], [2,1,4,3,6,5,8,7], [3,4,5,6,7,8,1,2], [4,3,6,5,8,7,2,1], [5,6,7,8,1,2,3,4], [6,5,8,7,2,1,4,3], [7,8,1,2,3,4,5,6], [8,7,2,1,4,3,6,5]]]
        rots = [[[0,1,2,3,4,5,6,7],[0,1,2,3]],[[0,1,2,3,4,5,6,7],[2,1,0,3]],[[0,1,6,7,4,5,2,3],[0,3,2,1]],[[0,1,6,7,4,5,2,3],[2,3,0,1]]]
        (narrows,arrow_types,Concs) = how_many_arrows(col)
        agroup = a_group(trans,rots)
        out = 55552
        self.assertEqual(polya(Concs,agroup,arrowings=arrow_types),out)

    # to long for testing purposes
    # def test_p8(self):
    #     from phenum.grouptheory import a_group
    #     from phenum.phonons import how_many_arrows, how_many_arrows
    #     dim = 4
    #     col = [[-1,1],[-1,1],[-1,1],[-1,1],[-1,1],[-1,1],[-1,1],[-1,1],[-1,1],[-1,2],[-1,2],[-1,2],[-1,2],[-1,2],[-1,2],[-1,2],[-1,2],[-1,2]]
    #     trans =[[j - 1 for j in i] for i in[[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18], [2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,1], [3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,1,2], [4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,1,2,3], [5,6,7,8,9,10,11,12,13,14,15,16,17,18,1,2,3,4], [6,7,8,9,10,11,12,13,14,15,16,17,18,1,2,3,4,5], [7,8,9,10,11,12,13,14,15,16,17,18,1,2,3,4,5,6], [8,9,10,11,12,13,14,15,16,17,18,1,2,3,4,5,6,7], [9,10,11,12,13,14,15,16,17,18,1,2,3,4,5,6,7,8], [10,11,12,13,14,15,16,17,18,1,2,3,4,5,6,7,8,9], [11,12,13,14,15,16,17,18,1,2,3,4,5,6,7,8,9,10], [12,13,14,15,16,17,18,1,2,3,4,5,6,7,8,9,10,11], [13,14,15,16,17,18,1,2,3,4,5,6,7,8,9,10,11,12], [14,15,16,17,18,1,2,3,4,5,6,7,8,9,10,11,12,13], [15,16,17,18,1,2,3,4,5,6,7,8,9,10,11,12,13,14], [16,17,18,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15], [17,18,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16], [18,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17]]]
    #     rots = [[[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17],[0,1,2,3]],[[0,17,16,15,14,13,12,11,10,9,8,7,6,5,4,3,2,1],[2,3,0,2]]]
    #     (narrows,arrow_types,Concs) = how_many_arrows(col)
    #     agroup = a_group(trans,rots)
    #     out = 1387
    #     self.assertEqual(polya(Concs,agroup,arrowings=arrow_types),out)

    def test_o22(self):
        from phenum.grouptheory import a_group
        from phenum.phonons import how_many_arrows, how_many_arrows
        dim = 4
        col = [[-1,1],[-1,1],[-1,1],[-1,2],[-1,2],[-1,2],[3,3],[3,4],[3,4]]
        trans =[[j - 1 for j in i] for i in [[1, 2, 3, 4, 5, 6, 7, 8, 9], [2, 3, 1, 5, 6, 4, 8, 9, 7], [3, 1, 2, 6, 4, 5, 9, 7, 8], [4, 5, 6, 7, 8, 9, 1, 2, 3], [5, 6, 4, 8, 9, 7, 2, 3, 1], [6, 4, 5, 9, 7, 8, 3, 1, 2], [7, 8, 9, 1, 2, 3, 4, 5, 6], [8, 9, 7, 2, 3, 1, 5, 6, 4], [9, 7, 8, 3, 1, 2, 6, 4, 5]]]
        rots = [[[0,1,2,3,4,5,6,7,8],[0,1,2,3]],[[0,1,2,8,6,7,4,5,3],[2,1,0,3]],[[0,2,1,4,3,5,8,7,6],[0,3,2,1]],[[0,2,1,6,8,7,3,5,4],[2,3,0,1]],[[0,7,5,6,4,2,3,1,8],[1,0,3,2]],[[0,7,5,8,3,1,4,2,6],[3,0,1,2]],[[0,5,7,4,6,2,8,1,3],[1,2,3,0]],[[0,5,7,3,8,1,6,2,4],[3,2,1,0]]]
        (narrows,arrow_types,Concs) = how_many_arrows(col)
        agroup = a_group(trans,rots)
        out = 4504
        self.assertEqual(polya(Concs,agroup,arrowings=arrow_types),out)
        
    # To long for testing purposes
    # def test_p7(self):
    #     from phenum.grouptheory import a_group
    #     from phenum.phonons import how_many_arrows, how_many_arrows
    #     dim = 4
    #     col = [[-1,1],[-1,1],[-1,1],[-1,1],[-1,1],[-1,1],[-1,1],[-1,1],[-1,2],[-1,2],[-1,2],[-1,2],[-1,2],[-1,2],[-1,2],[-1,2],[0,3],[0,4]]
    #     trans =[[j - 1 for j in i] for i in[[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18], [2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,1], [3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,1,2], [4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,1,2,3], [5,6,7,8,9,10,11,12,13,14,15,16,17,18,1,2,3,4], [6,7,8,9,10,11,12,13,14,15,16,17,18,1,2,3,4,5], [7,8,9,10,11,12,13,14,15,16,17,18,1,2,3,4,5,6], [8,9,10,11,12,13,14,15,16,17,18,1,2,3,4,5,6,7], [9,10,11,12,13,14,15,16,17,18,1,2,3,4,5,6,7,8], [10,11,12,13,14,15,16,17,18,1,2,3,4,5,6,7,8,9], [11,12,13,14,15,16,17,18,1,2,3,4,5,6,7,8,9,10], [12,13,14,15,16,17,18,1,2,3,4,5,6,7,8,9,10,11], [13,14,15,16,17,18,1,2,3,4,5,6,7,8,9,10,11,12], [14,15,16,17,18,1,2,3,4,5,6,7,8,9,10,11,12,13], [15,16,17,18,1,2,3,4,5,6,7,8,9,10,11,12,13,14], [16,17,18,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15], [17,18,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16], [18,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17]]]
    #     rots = [[[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17],[0,1,2,3]],[[0,17,16,15,14,13,12,11,10,9,8,7,6,5,4,3,2,1],[2,3,0,2]]]
    #     (narrows,arrow_types,Concs) = how_many_arrows(col)
    #     agroup = a_group(trans,rots)
    #     out = _read_output("p7")
    #     self.assertEqual(tree.brancher(Concs,agroup,col,dim),out)
        
    def test_o10(self):
        from phenum.grouptheory import a_group
        from phenum.phonons import how_many_arrows, how_many_arrows
        dim = 4
        col = [[-1,1],[-1,1],[-1,1],[-1,1],[-1,1],[-1,2],[-1,2],[-1,2],[-1,2],[-1,2],[3,3],[3,4]]
        trans =[[j - 1 for j in i] for i in[[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12], [2, 1, 4, 3, 6, 5, 8, 7, 10, 9, 12, 11], [3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 1, 2], [4, 3, 6, 5, 8, 7, 10, 9, 12, 11, 2, 1], [5, 6, 7, 8, 9, 10, 11, 12, 1, 2, 3, 4], [6, 5, 8, 7, 10, 9, 12, 11, 2, 1, 4, 3], [7, 8, 9, 10, 11, 12, 1, 2, 3, 4, 5, 6], [8, 7, 10, 9, 12, 11, 2, 1, 4, 3, 6, 5], [9, 10, 11, 12, 1, 2, 3, 4, 5, 6, 7, 8], [10, 9, 12, 11, 2, 1, 4, 3, 6, 5, 8, 7], [11, 12, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10], [12, 11, 2, 1, 4, 3, 6, 5, 8, 7, 10, 9]]]
        rots = [[[0,1,2,3,4,5,6,7,8,9,10,11],[0,1,2,3]],[[0,1,2,3,4,5,6,7,8,9,10,11],[2,1,0,3]],[[0,1,10,11,8,9,6,7,4,5,2,3],[0,3,2,1]],[[0,1,10,11,8,9,6,7,4,5,2,3],[2,3,0,1]]]
        (narrows,arrow_types,Concs) = how_many_arrows(col)
        agroup = a_group(trans,rots)
        out = 13896
        self.assertEqual(polya(Concs,agroup,arrowings=arrow_types),out)
        
    def test_o25(self):
        from phenum.grouptheory import a_group
        from phenum.phonons import how_many_arrows, how_many_arrows
        dim = 4
        col = [[-1,1],[-1,2],[-1,2],[-1,2],[3,3],[3,3],[3,4],[3,4]]
        trans =[[j - 1 for j in i] for i in[[1, 2, 3, 4, 5, 6, 7, 8], [2, 1, 4, 3, 6, 5, 8, 7], [3, 4, 5, 6, 7, 8, 1, 2], [4, 3, 6, 5, 8, 7, 2, 1], [5, 6, 7, 8, 1, 2, 3, 4], [6, 5, 8, 7, 2, 1, 4, 3], [7, 8, 1, 2, 3, 4, 5, 6], [8, 7, 2, 1, 4, 3, 6, 5]]]
        rots = [[[0,1,2,3,4,5,6,7],[0,1,2,3]],[[0,1,2,3,4,5,6,7],[2,1,0,3]],[[0,1,6,7,4,5,2,3],[0,3,2,1]],[[0,1,6,7,4,5,2,3],[2,3,0,1]]]
        (narrows,arrow_types,Concs) = how_many_arrows(col)
        agroup = a_group(trans,rots)
        out = 14344
        self.assertEqual(polya(Concs,agroup,arrowings=arrow_types),out)
        
    def test_o21(self):
        from phenum.grouptheory import a_group
        from phenum.phonons import how_many_arrows, how_many_arrows
        dim = 4
        col = [[-1,1],[-1,1],[-1,2],[-1,2],[-1,2],[3,3],[3,4],[3,4]]
        trans =[[j - 1 for j in i] for i in[[1, 2, 3, 4, 5, 6, 7, 8], [2, 1, 4, 3, 6, 5, 8, 7], [3, 4, 5, 6, 7, 8, 1, 2], [4, 3, 6, 5, 8, 7, 2, 1], [5, 6, 7, 8, 1, 2, 3, 4], [6, 5, 8, 7, 2, 1, 4, 3], [7, 8, 1, 2, 3, 4, 5, 6], [8, 7, 2, 1, 4, 3, 6, 5]]]
        rots = [[[0,1,2,3,4,5,6,7],[0,1,2,3]],[[0,1,2,3,4,5,6,7],[2,1,0,3]],[[0,1,6,7,4,5,2,3],[0,3,2,1]],[[0,1,6,7,4,5,2,3],[2,3,0,1]]]
        (narrows,arrow_types,Concs) = how_many_arrows(col)
        agroup = a_group(trans,rots)
        out = 3808
        self.assertEqual(polya(Concs,agroup,arrowings=arrow_types),out)

    # to long for testing purposes
    # def test_p6(self):
    #     from phenum.grouptheory import a_group
    #     from phenum.phonons import how_many_arrows, how_many_arrows
    #     dim = 4
    #     col = [[-3,1],[-3,1],[-3,1],[-3,1],[-3,1],[-3,1],[-3,1],[-3,1],[-3,1],[-3,1],[-3,2],[-3,2],[-3,2],[-3,2],[-3,2],[-3,2],[-3,2],[-3,2]]
    #     trans =[[j - 1 for j in i] for i in[[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18], [2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,1], [3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,1,2], [4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,1,2,3], [5,6,7,8,9,10,11,12,13,14,15,16,17,18,1,2,3,4], [6,7,8,9,10,11,12,13,14,15,16,17,18,1,2,3,4,5], [7,8,9,10,11,12,13,14,15,16,17,18,1,2,3,4,5,6], [8,9,10,11,12,13,14,15,16,17,18,1,2,3,4,5,6,7], [9,10,11,12,13,14,15,16,17,18,1,2,3,4,5,6,7,8], [10,11,12,13,14,15,16,17,18,1,2,3,4,5,6,7,8,9], [11,12,13,14,15,16,17,18,1,2,3,4,5,6,7,8,9,10], [12,13,14,15,16,17,18,1,2,3,4,5,6,7,8,9,10,11], [13,14,15,16,17,18,1,2,3,4,5,6,7,8,9,10,11,12], [14,15,16,17,18,1,2,3,4,5,6,7,8,9,10,11,12,13], [15,16,17,18,1,2,3,4,5,6,7,8,9,10,11,12,13,14], [16,17,18,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15], [17,18,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16], [18,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17]]]
    #     rots = [[[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17],[0,1,2,3]],[[0,17,16,15,14,13,12,11,10,9,8,7,6,5,4,3,2,1],[2,3,0,2]]]
    #     (narrows,arrow_types,Concs) = how_many_arrows(col)
    #     agroup = a_group(trans,rots)
    #     out = 1282
    #     self.assertEqual(polya(Concs,agroup,arrowings=arrow_types),out)

    # to large for testing purposes
    # def test_p5(self):
    #     from phenum.grouptheory import a_group
    #     from phenum.phonons import how_many_arrows, how_many_arrows
    #     dim = 4
    # col = [[-3,1],[-3,1],[-3,1],[-3,1],[-3,1],[-3,1],[-3,1],[-3,1],[-3,2],[-3,2],[-3,2],[-3,2],[-3,2],[-3,2],[3,3],[3,3],[3,4],[3,4]]
    #     trans =[[j - 1 for j in i] for i in[[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18], [2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,1], [3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,1,2], [4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,1,2,3], [5,6,7,8,9,10,11,12,13,14,15,16,17,18,1,2,3,4], [6,7,8,9,10,11,12,13,14,15,16,17,18,1,2,3,4,5], [7,8,9,10,11,12,13,14,15,16,17,18,1,2,3,4,5,6], [8,9,10,11,12,13,14,15,16,17,18,1,2,3,4,5,6,7], [9,10,11,12,13,14,15,16,17,18,1,2,3,4,5,6,7,8], [10,11,12,13,14,15,16,17,18,1,2,3,4,5,6,7,8,9], [11,12,13,14,15,16,17,18,1,2,3,4,5,6,7,8,9,10], [12,13,14,15,16,17,18,1,2,3,4,5,6,7,8,9,10,11], [13,14,15,16,17,18,1,2,3,4,5,6,7,8,9,10,11,12], [14,15,16,17,18,1,2,3,4,5,6,7,8,9,10,11,12,13], [15,16,17,18,1,2,3,4,5,6,7,8,9,10,11,12,13,14], [16,17,18,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15], [17,18,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16], [18,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17]]]
    #     rots = [[[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17],[0,1,2,3]],[[0,17,16,15,14,13,12,11,10,9,8,7,6,5,4,3,2,1],[2,3,0,2]]]
    #     (narrows,arrow_types,Concs) = how_many_arrows(col)
    #     agroup = a_group(trans,rots)
    #     self.assertEqual(tree.brancher(Concs,agroup,col,dim),out)
        
    def test_p4(self):
        from phenum.grouptheory import a_group
        from phenum.phonons import how_many_arrows, how_many_arrows
        dim = 4
        col = [[-1,3],[-1,2],[-1,3],[1,1]]
        trans = [[0,1,2,3],[1,2,3,0],[2,3,0,1],[3,0,1,2]]
        rots = [[[0,1,2,3],[0,1,2,3]],[[0,1,2,3],[2,1,0,3]],[[0,3,2,1],[0,3,2,1]],[[0,3,2,1],[2,3,0,1]]]
        (narrows,arrow_types,Concs) = how_many_arrows(col)
        agroup = a_group(trans,rots)
        out = 5
        self.assertEqual(polya(Concs,agroup,arrowings=arrow_types),out)
        
    def test_p3(self):
        from phenum.grouptheory import a_group
        from phenum.phonons import how_many_arrows, how_many_arrows
        dim = 4
        col = [[-1,1],[1,1],[1,2],[1,3],[1,1],[1,2]]
        trans = [[0,1,2,3,4,5],[1,2,3,4,5,0],[2,3,4,5,0,1],[3,4,5,0,1,2],[4,5,0,1,2,3],[5,0,1,2,3,4]]
        rots = [[[0,1,2,3,4,5],[0,1,2,3]],[[0,1,2,3,4,5],[2,1,0,3]],[[0,5,4,3,2,1],[0,3,2,1]],[[0,5,4,3,2,1],[2,3,0,1]]]
        (narrows,arrow_types,Concs) = how_many_arrows(col)
        agroup = a_group(trans,rots)
        out = 7936
        self.assertEqual(polya(Concs,agroup,arrowings=arrow_types),out)
        
    def test_p2(self):
        from phenum.grouptheory import a_group
        from phenum.phonons import how_many_arrows, how_many_arrows
        dim = 4
        col = [[-3,1],[-3,1],[-3,1],[-3,2],[3,3],[3,3],[3,3],[3,4]]
        trans = [[0,1,2,3,4,5,6,7],[1,0,3,2,5,4,7,6],[2,3,4,5,6,7,0,1],[3,2,5,4,7,6,1,0],[4,5,6,7,0,1,2,3],[5,4,7,6,1,0,3,2],[6,7,0,1,2,3,4,5],[7,6,1,0,3,2,5,4]]
        rots = [[[0,1,2,3,4,5,6,7],[0,1,2,3]],[[0,1,2,3,4,5,6,7],[2,1,0,3]],[[0,1,6,7,4,5,2,3],[0,3,2,1]],[[0,1,6,7,4,5,2,3],[2,3,0,1]]]
        (narrows,arrow_types,Concs) = how_many_arrows(col)
        agroup = a_group(trans,rots)
        out = 9568
        self.assertEqual(polya(Concs,agroup,arrowings=arrow_types),out)
        
    def test_p1(self):
        from phenum.grouptheory import a_group
        from phenum.phonons import how_many_arrows, how_many_arrows
        dim = 4
        col = [[-1,1],[-1,2],[0,3],[0,4]]
        trans = [[0,1,2,3],[1,0,3,2],[2,3,0,1],[3,2,1,0]]
        rots = [[[0,1,2,3],[0,1,2,3]],[[0,1,2,3],[2,1,0,3]],[[0,1,2,3],[0,3,2,1]],[[0,1,2,3],[2,3,0,1]],[[0,3,2,1],[1,0,3,2]],[[0,3,2,1],[3,0,1,2]],[[0,3,2,1],[1,2,3,0]],[[0,3,2,1],[3,2,1,0]]]
        (narrows,arrow_types,Concs) = how_many_arrows(col)
        agroup = a_group(trans,rots)
        out = 18
        self.assertEqual(polya(Concs,agroup,arrowings=arrow_types),out)
        
    def test_r1(self):
        from phenum.grouptheory import a_group
        from phenum.phonons import how_many_arrows, how_many_arrows
        dim = 4
        col = [[-1,1],[-1,2],[0,3],[0,4]]
        trans = [[0,1,2,3],[1,0,3,2],[2,3,0,1],[3,2,1,0]]
        rots = [[[0,1,2,3],[0,1,2,3]],[[0,1,2,3],[2,1,0,3]],[[0,1,2,3],[0,3,2,1]],[[0,1,2,3],[2,3,0,1]],[[0,3,2,1],[1,0,3,2]],[[0,3,2,1],[3,0,1,2]],[[0,3,2,1],[1,2,3,0]],[[0,3,2,1],[3,2,1,0]]]
        (narrows,arrow_types,Concs) = how_many_arrows(col)
        temp_agroup = a_group(trans,rots)
        agroup = []
        print("ta",temp_agroup)
        for i in temp_agroup:
            sites = [j+1 for j in i[0]]
            ars = [j+1 for j in i[1]]
            agroup.append([sites,ars])
        print("a",agroup)
        out = 18
        self.assertEqual(polya(Concs,agroup,arrowings=arrow_types),out)
        
    def test_r2(self):
        from phenum.grouptheory import a_group
        from phenum.phonons import how_many_arrows, how_many_arrows
        dim = 4
        col = [[-1,1],[-1,2],[0,3],[0,4]]
        trans = [[0,1,2,3],[1,0,3,2],[2,3,0,1],[3,2,1,0]]
        rots = [[[0,1,2,3],[0,1,2,3]],[[0,1,2,3],[2,1,0,3]],[[0,1,2,3],[0,3,2,1]],[[0,1,2,3],[2,3,0,1]],[[0,3,2,1],[1,0,3,2]],[[0,3,2,1],[3,0,1,2]],[[0,3,2,1],[1,2,3,0]],[[0,3,2,1],[3,2,1,0]]]
        (narrows,arrow_types,Concs) = how_many_arrows(col)
        Concs = [1]
        temp_agroup = a_group(trans,rots)
        agroup = []
        for i in temp_agroup:
            sites = [j+1 for j in i[0]]
            ars = [j+1 for j in i[1]]
            agroup.append([sites,ars])

        out = 18
        with pytest.raises(ValueError):
            polya(Concs,agroup,arrowings=arrow_types)

if __name__ == '__main__':
    untittest.main()
