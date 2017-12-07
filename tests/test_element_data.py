"""Tests the subroutines in element_data.
"""

import pytest
import numpy as np

def test_get_lattice_parameter():
    """Tests the get_lattice_parameter subroutine.
    """
    from phenum.element_data import get_lattice_parameter
    
    (a,title) = get_lattice_parameter(["Al","Ni"],[1,1],[[0.5,0.5,0],[0.5,0,0.5],[0,0.5,0.5]],
                                      1,"test")
    assert np.allclose(a,3.785)
    assert title == " Al  Ni  test\n"

    (a,title) = get_lattice_parameter(None,[2,2],[[0.5,0.5,0],[0.5,0,0.5],[0,0.5,0.5]],
                                      1,"test")
    assert a == 1.0
    assert title == "test"

    with pytest.raises(ValueError):
        get_lattice_parameter(["Al","Ni"],[1,2,3],[[0.5,0.5,0],[0.5,0,0.5],[0,0.5,0.5]],
                                      1,"test")
        
    (a,title) = get_lattice_parameter(["C","Cr"],[1,1],[[0.5,0.5,0],[0.5,0,0.5],[0,0.5,0.5]],
                                      2,"test")
    assert np.allclose(a,4.07086)
    assert title == " C  Cr  test\n"

