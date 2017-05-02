"""Methods for testing the subroutines in the structures module."""
import unittest as ut
import pytest

class TestDistribution(ut.TestCase):
    """Tests the _distribution subroutine."""

    def test1(self):
        from phenum.structures import _distribution
        with pytest.raises(ValueError):
            _distribution("shape",None,10,cellsizes=[1,2],res_type="stuff")

    def test2(self):
        from phenum.structures import _distribution
        with pytest.raises(ValueError):
            _distribution("shape",None,10)
        
    def test3(self):
        from phenum.structures import _distribution
        with pytest.raises(ValueError):
            _distribution("clean",[1,2,3],10)

class TestMakeEnumIn(ut.TestCase):
    """Tests the _distribution_summary subroutine."""

    def test1(self):
        from phenum.structures import make_enum_in
        with pytest.raises(ValueError):
            make_enum_in('all','enum.in',sizes=[1,2,3])
        
    def test2(self):
        from phenum.structures import make_enum_in
        with pytest.raises(ValueError):
            make_enum_in('all','enum.in')
        
