"""Tests the subroutines in element_data.
"""

import pytest

def test_get_lattice_parameter():
    """Tests the get_lattice_parameter subroutine.
    """
    from phenum.element_data import get_lattice_parameter
    
    (a,title) = get_lattice_parameter(["Al","Ni"],[1,1],"test")
    assert a == 3.785
    assert title == "test  Al  Ni \n"

    (a,title) = get_lattice_parameter(None,[2,2],"test")
    assert a == 1.0
    assert title == "test"

    with pytest.raises(ValueError):
        get_lattice_parameter(["Al","Ni"],[1,2,3],"test")
