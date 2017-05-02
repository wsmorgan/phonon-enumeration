"""Methods for testing the subroutines in the HNFs module."""
import unittest as ut
import numpy as np

class TestGetHNFDiagonals(ut.TestCase):
    """Tests the _get_HNF_diagonals subroutine."""

    def test1(self):
        from phenum.HNFs import _get_HNF_diagonals
        n = 4
        out = [[1,1,4],[1,2,2],[1,4,1],[2,1,2],[2,2,1],[4,1,1]]
        test = _get_HNF_diagonals(n)
        self.assertEqual(test,out)

    def test2(self):
        from phenum.HNFs import _get_HNF_diagonals
        n = 15
        out = [[1,1,15],[1,3,5],[1,5,3],[1,15,1],[3,1,5],[3,5,1],[5,1,3],[5,3,1],[15,1,1]]
        test = _get_HNF_diagonals(n)
        self.assertEqual(test,out)
