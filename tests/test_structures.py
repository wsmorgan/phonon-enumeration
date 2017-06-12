"""Methods for testing the subroutines in the structures module."""
import unittest as ut
import pytest

class TestDistribution(ut.TestCase):
    """Tests the _distribution subroutine."""

    def test_1(self):
        from phenum.structures import _distribution
        with pytest.raises(ValueError):
            _distribution('polya.out',"shape",None,10,cellsizes=[1,2],res_type="stuff")

    # def test_2(self):
    #     from phenum.structures import _distribution
    #     with pytest.raises(ValueError):
    #         _distribution('polya.out',"shape",None,10)
        
    def test_3(self):
        from phenum.structures import _distribution
        with pytest.raises(ValueError):
            _distribution('polya.out',"clean",[1,2,3],10)

class TestMakeEnumIn(ut.TestCase):
    """Tests the _distribution_summary subroutine."""

    def test_1(self):
        from phenum.structures import make_enum_in
        with pytest.raises(ValueError):
            make_enum_in('polya.out','all','enum.in',sizes=[1,2,3])
        
    def test_2(self):
        from phenum.structures import make_enum_in
        with pytest.raises(ValueError):
            make_enum_in('polya.out','all','enum.in')
        
