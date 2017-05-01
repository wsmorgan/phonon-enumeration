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
            q = n/i
            for j in range(1,q+1):
                if not q%j == 0:
                    continue
                else:
                    diags.append([i,j,q/j])
                    
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
        a = diag[0]
        c = diag[1]
        f = diag[2]
        
        for b in range(c):
            for e in range(f):
                for d in range(f):
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
    from copy import deepcopy

    if eps_:
        eps = eps_
    else:
        eps = 1E-7

    aType = [1]*len(base_vecs)
    
    (rots, shifts) = get_spaceGroup(pLV,aType,base_vecs,eps=eps)

    temp_hnf = deepcopy(HNFs)

    uHNFs = []
    if len(HNFs) == 1:
        uHNFs = HNFs
    else:
        for i in range(1,len(HNFs)):
            duplicate = False
            for j in range(len(uHNFs)):
                for r in range(len(rots)):
                    test_latticei = np.matmul(rots[r],np.matmul(pLV,HNFs[i]))
                    test_latticej = np.matmul(pLV,uHNFs[j])
                    
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
    if eps_:
        eps = eps_
    else:
        eps = 1E-7

    HNFs = find_all_HNFs(n)
    base_perms = _get_dvector_permutations(pLV,base_vecs,LatDim,eps)

    (HNFs,latts,fix_Ops,RPList,degens) = remove_duplicate_lattices(HNFs,pLV,base_vecs,
                                                                   LatDim,base_perms,eps_=eps)

    return HNFs

def get_SNFs(HNFs,fix_ops,RPList):
    """Finds the unique SNFs, and right and left transforms. It also reorders all the inputs.

    Args:
        HNFs (array-like): A list of the unique HNFs.
        fix_ops (opList): The lattice fixing operations.
        RPList (RotPermList): The permutations of the basis vectors.

    Returns:
        HNFs (array-like): A list of unique HNFs.
        SNFs (array-like): A list of unique SNFs.
        Ls (array-like): A list of the left Transforms.
        SNF_labels (list): The integer SNF labels.
        fix_ops (opList): The lattice fixing operations.
        RPList (opList): The permutations of the basis set.
    """

    from .grouptheory import SmithNormalForm
    
    nfound = 0

    uSNFs = []
    SNFs = []
    Ls = []
    SNF_labels = np.zeros(len(HNFs))
    for ihnf in range(len(HNFs)):
        (SNF,L,R) = SmithNormalForm(HNF[ihnf])
        SNFs.append(SNF)
        Ls.append(L)
        duplicate = False
        for ifound in range(nfound):
            if SNF == uSNFs[ifound]:
                duplicate = True
                SNF_labels[ihnf] = ifound
                break

        if not duplicate:
            nfound += 1
            uSNFs.append(SNF)
            SNF_labels[ihnf] = nfound

    # Now we want to reorder the groups so that they each come in the
    # same order and in actual groups instead of scattered.
    idx = []
    for i in range(nfound):
        for j in range(len(SNF_labels)):
            if SNF_labels[j] == i:
                idx.append(j)

    HNFs = np.array(HNFs)
    SNFs = np.array(SNFs)
    Ls = np.array(Ls)
    SNF_labels = np.array(SNF_labels)
    fix_ops = np.array(fix_ops)
    RPList = np.array(RPlist)

    HNFs = HNFs[idx]
    SNFs = SNFs[idx]
    Ls = Ls[idx]
    SNF_labels = SNF_labels[idx]
    fix_ops = fix_ops[idx]
    RPList = RPList[idx]

    return(HNFs,SNFs,Ls,SNF_labels,fix_ops,RPList)
