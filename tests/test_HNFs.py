"""Methods for testing the subroutines in the HNFs module."""
import unittest as ut
import numpy as np

class TestGetHNFDiagonals(ut.TestCase):
    """Tests the _get_HNF_diagonals subroutine."""

    def test_1(self):
        from phenum.HNFs import _get_HNF_diagonals
        n = 4
        out = [[1,1,4],[1,2,2],[1,4,1],[2,1,2],[2,2,1],[4,1,1]]
        test = _get_HNF_diagonals(n)
        self.assertEqual(test,out)

    def test_2(self):
        from phenum.HNFs import _get_HNF_diagonals
        n = 15
        out = [[1,1,15],[1,3,5],[1,5,3],[1,15,1],[3,1,5],[3,5,1],[5,1,3],[5,3,1],[15,1,1]]
        test = _get_HNF_diagonals(n)
        self.assertEqual(test,out)

class TestGetHNFs(ut.TestCase):
    """Tests the get_HNFs subroutine."""
    
    def test_get_HNFs1(self):
        from phenum.HNFs import get_HNFs
        n = 10
        pLV = [[1,0,0],[0,1,0],[0,0,1]]
        base_vecs = [[0,0,0]]
        LatDim = 3
        out = [[[1, 0, 0], [0, 1, 0], [0, 0, 10]], [[1, 0, 0], [0, 1, 0], [0, 1, 10]], [[1, 0, 0], [0, 1, 0], [0, 2, 10]], [[1, 0, 0], [0, 1, 0], [0, 3, 10]], [[1, 0, 0], [0, 1, 0], [0, 4, 10]], [[1, 0, 0], [0, 1, 0], [0, 5, 10]], [[1, 0, 0], [0, 1, 0], [1, 1, 10]], [[1, 0, 0], [0, 1, 0], [1, 2, 10]], [[1, 0, 0], [0, 1, 0], [1, 3, 10]], [[1, 0, 0], [0, 1, 0], [1, 4, 10]], [[1, 0, 0], [0, 1, 0], [1, 5, 10]], [[1, 0, 0], [0, 1, 0], [2, 2, 10]], [[1, 0, 0], [0, 1, 0], [2, 3, 10]], [[1, 0, 0], [0, 1, 0], [2, 4, 10]], [[1, 0, 0], [0, 1, 0], [2, 5, 10]], [[1, 0, 0], [0, 1, 0], [3, 5, 10]], [[1, 0, 0], [0, 1, 0], [4, 4, 10]], [[1, 0, 0], [0, 1, 0], [4, 5, 10]], [[1, 0, 0], [0, 1, 0], [5, 5, 10]], [[1, 0, 0], [0, 2, 0], [0, 0, 5]], [[1, 0, 0], [0, 2, 0], [1, 0, 5]], [[1, 0, 0], [0, 2, 0], [2, 0, 5]], [[1, 0, 0], [1, 2, 0], [0, 0, 5]]]
        test = get_HNFs(n,pLV,base_vecs,LatDim,eps_=1E-7)
        self.assertEqual(test,out)
        
