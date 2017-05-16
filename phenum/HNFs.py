"""Subtroutines needed to find the HNFs."""
import numpy as np

def _get_HNF_diagonals(n):
    """Finds the diagonals of the HNF that reach the target n value.
    
    Args:
        n (int): The target determinant for the HNF.
        
    Retruns:
        diags (list of lists): The allowed values of the determinant.
    """
    
    diags = []
    for i in range(1,n+1):
        if not n%i == 0:
            continue
        else:
            q = int(n//i)
            for j in range(1,q+1):
                if not q%j == 0:
                    continue
                else:
                    diags.append([i,j,q//j])
                    
    return diags

def find_all_HNFs(n):
    """Finds the HNFs for a given volume factor.
    
    Args:
        n (int): The determinant of the HNFs.

    Retruns:
        HNFs (array-like): A list of the HNFs.
    """

    diags = _get_HNF_diagonals(n)

    HNFs = []
    for diag in diags:
        a = int(diag[0])
        c = int(diag[1])
        f = int(diag[2])
        
        for b in range(c):
            for d in range(f):
                for e in range(f):
                    HNF = [[a,0,0],[b,c,0],[d,e,f]]
                    HNFs.append(HNF)

    return HNFs

def remove_duplicate_lattices(HNFs,pLV,base_vecs,LatDim,base_perms,eps_=None):
    """Removes symmetrically equivalent lattices from the list of possible lattices.

    Args:
        HNFs (array-like): A list of the possible HNFs.
        pLV (array-like): The parent lattice vectors as a column matrix.
        base_vecs (array-like): The atomic basis vectors.
        LatDim (int): 3 if bulk, 2 if surf.
        base_perms (RotPermList): Permutations of the basis vectors.
        eps_ (float optional): Floating point tollerance.

    Returns:
        uHNFs (array-like): List of unique HNFs.
        latts (array-like): List of unique lattices.
        fix_Ops (opList): List of operations that preserve the symmetry group.
        RPList (RotPermList): Permutations of basis vectors.
        degeneracy (list): The degeneracy of each HNF.
    """

    from .symmetry import get_spaceGroup
    from .grouptheory import _is_equiv_lattice, _get_sLV_fixing_operations

    if eps_:
        eps = eps_
    else: #pragma: no cover
        eps = 1E-7

    aType = [1]*len(base_vecs)

    (rots, shifts) = get_spaceGroup(pLV,aType,base_vecs,eps=eps)
    A = np.transpose(pLV)
    
    uHNFs = []
    if len(HNFs) == 1:
        uHNFs = HNFs
    else:
        for i in range(len(HNFs)):
            duplicate = False
            for j in range(len(uHNFs)):
                for r in range(len(rots)):
                    test_latticei = np.matmul(rots[r],np.matmul(A,HNFs[i]))
                    test_latticej = np.matmul(A,uHNFs[j])
                    
                    if _is_equiv_lattice(test_latticei,test_latticej,eps):
                        duplicate = True
                        break
                if duplicate == True:
                    break

            if not duplicate:
                uHNFs.append(HNFs[i])

    latts = []
    fix_Ops = []
    RPList = []
    degens = []
    for HNF in uHNFs:
        B = np.matmul(pLV,HNF)
        latts.append(B)
        (fix_op, RPlisti,degeneracy) = _get_sLV_fixing_operations(HNF,pLV,len(base_vecs),
                                                                 rots,shifts,base_perms,eps)
        fix_Ops.append(fix_op)
        RPList.append(RPlisti)
        degens.append(degeneracy+1)

    return (uHNFs,latts,fix_Ops,RPList,degens)

def get_HNFs(n,pLV,base_vecs,LatDim,eps_=None):
    """Finds the unique HNFs for the system.

    Args:
        n (array-like): The determinant of the HNFs.
        pLV (array-like): The parent lattice vectors as a column matrix.
        base_vecs (array-like): The atomic basis vectors.
        LatDim (int): 3 if bulk, 2 if surf.
        eps_ (float optional): Floating point tollerance.

    Returns:
        HNFs (array-like): A list of the unique HNFS
    """

    from .grouptheory import _get_dvector_permutations
    from .symmetry import bring_into_cell
    from copy import deepcopy

    if eps_:
        eps = eps_
    else:
        eps = 1E-7

    HNFs = find_all_HNFs(n)
    temp_basis = deepcopy(base_vecs)
    pLVinv = np.linalg.inv(pLV)
    for i in range(len(base_vecs)):
        # par_lat_inv = linalg.inv(np.transpose(par_lat))
        base_vecs[i] = bring_into_cell(base_vecs[i],pLVinv,pLV,eps)
        if not np.allclose(base_vecs[i], temp_basis[i], rtol=0, atol=eps):
            from phenum.msg import warn
            warn("An atomic basis vector was not inside the unit cell. It has been "
                 "remapped.\n Original vector: {} \n Remapped vector: {}"
                 .format(" ".join([str(p) for p in temp_basis[i]]),
                         " ".join([str(p) for p in base_vecs[i]])))
            
    base_perms = _get_dvector_permutations(pLV,base_vecs,LatDim,eps)

    (HNFs,latts,fix_Ops,RPList,degens) = remove_duplicate_lattices(HNFs,pLV,base_vecs,
                                                                   LatDim,base_perms,eps_=eps)
    return HNFs
