"""Methods for taking an HNF and finding the site permutations and the
arrow permutations. All of these have been modified for python from
their original fortran implementations which can be found at:
https://github.com/msg-byu/enumlib/. Of these only get_rotation
get_rotation_perms_lists() has been modified in order to produce the
arrow group. The class ArrowPerm has also been added and implemented
inside the RotPermList class.
"""

import itertools
import numpy as np
import operator

class ArrowPerm(object):
    """ArrowPerm pairs a site permutation with an arrow permutation.

    Attributes:
        site_perm (list): A 1D array containing a site permutation.
        arrow_perm (list): A 1D array containing an arrow permutation.
    """
    def __init__(self,site_perm = None,arrow_perm = None):
        """Initializes the ArrowPerm.
        
        Args:
            site_perm (list): A 1D array containing a site permutation.
            arrow_perm (list): A 1D array containing an arrow permutation.
        """
        self.site_perm = site_perm
        self.arrow_perm = arrow_perm

class RotPermList(object):
    """RotPermList maintains the data structure for the rotations and
    permutations from the fortran code.

    Attributes:
        nL (int): An integer indicating the number of operations.
        v (array-like): A 2D array containing the lattice vectors.
        perm (ArrowPerm): An ArrowPerm object.
        RotIndx (list): List of the rotation indices.
    """
    def __init__(self,nL = None,v =None,perm = None,RotIndx = None,arrows=None):
        """Initializes the RotPermList.

        Args:
            nL (int): An integer indicating the number of operations.
            v (array-like): A 2D array containing the lattice vectors.
            perm (ArrowPerm): An ArrowPerm object.
            RotIndx (list):  List of the rotation indices.
        """
        self.nL = nL
        self.v = v
        self.perm = ArrowPerm(site_perm = perm, arrow_perm=arrows)
        self.RotIndx = RotIndx

class OpList(object):
    """OpList maintains the data structure of the fixing operations form
    the fortran code.

    Attributes:
        rot (array-like): A 3D array containing the rotations.
        shift (array-like): A 2D array containing the translations of the lattice.
    """
    def __init__(self,rot = None, shift = None):
        """Initializes the OpList.
        
        Args:
            rot (array-like): A 3D array containing the rotations.
            shift (array-like): An 2D array containing the translations of th lattice.
        """
        self.rot = rot
        self.shift = shift

def _make_member_list(n):
    """Takes the length of three cyclic group and constructs the member
    list so that the permutations can be determined. Each member has
    three components, corresponding to the entries for each of the
    three cyclic groups.
    
    Args:
        n (int): array containing the diagonal elements of the SNF.

    Returns:
        p (list): The member list for the symmetry group.
    """
    from operator import mul
    from functools import reduce
    
    depth = int(round(reduce(mul,n,1)))
    p = np.zeros([depth,3])
    for im in range(1,depth):  # Loop over the members of the translation group
        p[im] = list(p[im-1]) # Start with the same digits as in the previous increment        
        p[im][2] = (p[im-1][2]+1)%n[2]  # Increment the first cyclic group
        if (p[im][2]==0):             # If it rolled over then
            p[im][1] = (p[im-1][1] +1) % n[1] # increment the next cyclic group
            if (p[im][1]==0):          # If this one rolled over too
                p[im][0] = (p[im-1][0]+1)%n[0] # Then increment the third one
    return p

def _find_permutation_of_group(g,gp):
    """Constructs the permutation group from the point group.
    
    Args:
        g (list): The unpermuted symmetry group.
        gp (list): The permuted symmetry group.

    Returns:
        perm (list): Permutation of g to get gp.
    """
    n = len(g)
    perm = []
    skip = [False]*n
    for im in range(n):
        for jm in range(n):
            if skip[jm]:
                continue # This is just for efficiency
            if np.allclose(gp[jm],g[im]):
                perm.append(jm)
                skip[jm] = True
                break # don't keep looking if you already found the match
            
    return perm

def _is_equiv_lattice(lat1,lat2,eps):
    """This function determines whether two lattices are equivalent. That
    is, if one is an equal-volume derivative lattice of the other. It
    uses the idea of Santoro and Mighell (Acta. Cryst. 1972) that, for
    equivalent lattices, the transformation matrix that takes the
    vectors of str2 to str1 has determinant of 1 and all integer
    elements. [LV2]=[LV1][S]=&gt;S=LV1^-1*LV2. But I had to use abs on
    the determinant. This is not mentioned by Santoro and Mighell
    (issue of opposite handedness---which doesn't matter to me).

    Args:
        lat1 (array-like): 2D array containing the first lattice.
        lat2 (array-like): 2D array containing the second lattice.
        eps (float): The finite precision tolerance.

    Returns:
        is_equiv_lattice (bool): True if lattice are equivalent.
    """
    from numpy import linalg, allclose, matmul
    is_equiv_lattice = False
    lat1inv = linalg.inv(lat1)
    S = matmul(lat1inv,lat2).tolist()
    if (allclose(abs(linalg.det(S)),1.0,rtol=0,atol=eps) and
        allclose(S,[[int(round(S[i][j])) for j in range(3)] for i in range(3)],rtol=0,atol=eps)):
        is_equiv_lattice = True
    return is_equiv_lattice
        
def _get_sLV_fixing_operations(HNF,pLV,rot,shift,dPerm,eps):
    """This subroutine deteremines which operations of the parent cell
    don't change the supercell. These operations are the symmetry
    operations of the supercell.

    Args:
        HNF (list): A 2D integer list containing the hermite normal form matrix.
        pLV (array-like): A 2D array containing the parent lattice vectors.
        rot (array-like): A 3D array containing all the rotation matrices.
        shift (array-like): A 2D array containing the translations. 
        dPerm (list): List of basis vector permutations.
        eps (float): Finite precision tolerance.

    Returns:
        fix_op (OpList): A list of the lattice fixing symmetry operations.
        rot_perm (list): The rotation permutations.
        degeneracy (int): Degeneracy of the supercell.
    """

    n_rot = len(rot)
    degen_lattices = np.zeros([n_rot,3,3])
        
    c_degen = 0
    ic = 0 # Counter for the fixing operations
    tv = 0
    tmp_op_rot = []
    tmp_op_shift = []
    tv = []
    t_index = [] # temp variables
    for iRot in range(n_rot):  # Loop over each rotation
        this_rot = rot[iRot] # Store the rotation
        orig_lat = np.matmul(np.transpose(pLV),HNF)  # Compute the superlattice
        rot_lat = np.matmul(this_rot,orig_lat) # Compute the rotated superlattice
        if _is_equiv_lattice(rot_lat,orig_lat,eps):
            # this operation fixes the lattice and should be recorded
            ic += 1
            tmp_op_rot.append(this_rot)
            tmp_op_shift.append(shift[iRot])
            tv.append(dPerm.v[iRot])
            t_index.append(iRot)
            # Added by LN from here
        else:
            in_list = False
            for iDegen in range(c_degen):
                if _is_equiv_lattice(degen_lattices[iDegen],rot_lat,eps):
                    in_list = True
                    break
          
            if not in_list:
                degen_lattices[c_degen] = rot_lat
                c_degen += 1
        # End of LN additions for degeneracy
        # Now we know which rotations fix the lattice and how many
        # there are so store them
    degeneracy = c_degen
    # Allocate the storage for them
    fix_op = OpList(tmp_op_rot,tmp_op_shift) # Stuff the rotations into the permanent array

    rot_perm = RotPermList(v=tv,RotIndx=t_index)
    
    return(fix_op,rot_perm,degeneracy)

def _map_dvector_permutation(rd,d,eps):
    """Maps the basis vectors to a permutation.
    
    Args:
        rd (array-like): 2D array of the rotated basis vectors.
        d (array-like): 2D array of the original basis vectors.
        eps (float): Finite precision tolerance.

    Returns:
       RP (list): The permutation of the basis vectors.
    """
    n_d = len(rd) # of d-vectors
    found = [False]*n_d
    RP = []
    for iD in range(n_d):
        for jD in range(n_d):
            if found[jD]:
                continue
            if np.allclose(rd[iD],d[jD],atol=eps,rtol=0):
                RP.append(jD)
                found[jD] = True
                break

    if len(RP) != len(d): #pragma: no cover
        print("d-vector didn't permute in map_dvector_permutation "
              "This usually means that the d-set from the input structure and the d-set"
              " from the struct_enum.out have a different origin or don't live in the same"
              " unit cell. This probably isn't your fault---the code should overcome this."
              ,RP)
        exit()
    return(RP)

        
def _find_minmax_indices(invec):
    """Finds the indices corresponding to the minimum and maximum values
    in an integer vector.

    Args:
        invec (list): The input array.

    Returns:
        this_min (int): The minimum index.
        this_max (int): The maximum index.
    """

    vec = [abs(i) for i in invec]
    
    this_min = vec.index(min([i for i in vec if i > 0]))

    rvec = [i for i in reversed(vec)]
    
    this_max = 2 - rvec.index(max(rvec))

    return (this_min,this_max)

def SmithNormalForm(HNF):
    """Finds the Smith Normal Form (SNF) of an HNF as well as the left and
    right transforms.
      
    Args:
        HNF (list of list): An integer matrix for which the SNF is to be found.

    Returns:

        (M, A, B) (matrix, matrix, matrix): The M is the Smith Normal Form
        matrix of the input matrix, A is the left transform matrix,
        and B is the right transform matrix. The matrices are such
        that np.dot(np.dot(A,HNF),B) = M.

    """
    from numpy import dot
    from copy import deepcopy
    
    if np.linalg.det(HNF) < 1:
        raise ValueError("SmithNormalForm routine failed because the input matrix had a "
                         "determinant less than 1.")

    A = [[0,0,0],[0,0,0],[0,0,0]]
    Mpre = deepcopy(list(HNF))
    M = [[int(x) for x in i] for i in Mpre]
    B = [[0,0,0],[0,0,0],[0,0,0]]

    if not np.allclose(M,Mpre):
        raise ValueError("The input matrix to SmithNormalForm was not integer.")
    
    for i in range(3):
        A[i][i] = 1
        B[i][i] = 1

    j = 0
    it_cnt = 0

    stop_loop = False
    while not stop_loop:
        it_cnt += 1
        if (it_cnt >=100): #pragma: no cover
            raise RuntimeError("Bad programming in SmithNormalForm")
        
        while (3-[M[0][j],M[1][j],M[2][j]].count(0)) > 1:
            (minidx,maxidx) = _find_minmax_indices([M[0][j],M[1][j],M[2][j]])
            minm = M[minidx][j]
            mult = M[maxidx][j]//minm

            M[maxidx] = list(map(operator.sub,M[maxidx],[mult*i for i in M[minidx]]))
            A[maxidx] = list(map(operator.sub,A[maxidx],[mult*i for i in A[minidx]]))
            if dot(dot(A,HNF),B).any() != np.array(M).any(): #pragma: no cover
                print("ROW: Transformation matrices didn't work")
                exit()
        if M[j][j] == 0:
            maxidx = [abs(M[0][j]),abs(M[1][j]),abs(M[2][j])].index(max([abs(M[0][j]),abs(M[1][j]),abs(M[2][j])]))
            tmprow = list(A[j])
            A[j] = list(A[maxidx])
            A[maxidx] = tmprow

            tmprow = list(M[j])
            M[j] = list(M[maxidx])
            M[maxidx] = tmprow    
        if np.dot(np.dot(A,HNF),B).any() != np.array(M).any(): #pragma: no cover
            print("ROWSWAP: Transformation matrices didn't work")
            exit()

        while (3-M[j].count(0)) >1:
            (minidx,maxidx) = _find_minmax_indices(M[j])

            minm = M[j][minidx]
            mult = M[j][maxidx]//M[j][minidx]
            for i in range(3):
                M[i][maxidx] = M[i][maxidx]-mult * M[i][minidx]
                B[i][maxidx] = B[i][maxidx]-mult * B[i][minidx]

            if np.dot(np.dot(A,HNF),B).any() != np.array(M).any(): #pragma: no cover
                print("COLS: Transformation matrices didn't work")
                exit()

        if M[j][j] < 0:
            for i in range(3):
                M[i][j] = -M[i][j]
                B[i][j] = -B[i][j]
        elif M[j][j] == 0:
            maxidx = [abs(i) for i in M[j]].index(max([abs(i) for i in M[j]]))
            for i in range(3):
                tmp = B[i][j]
                B[i][j] = B[i][maxidx]
                B[i][maxidx] = tmp 
                tmp = M[i][j]
                M[i][j] = M[i][maxidx]
                M[i][maxidx] = tmp

        if ((3-M[j].count(0)) >1) or ((3-[M[0][j],M[1][j],M[2][j]].count(0)) >1):
            continue
        if np.dot(np.dot(A,HNF),B).any() != np.array(M).any(): #pragma: no cover
            print("COLSWAP: Transformation matrices didn't work")
            exit()

        l_div = []
        for i in range(1,3):
            for k in range(1,3):
                l_div.append((M[i][k]%M[0][0] == 0))

        if j == 0 and False in l_div:
            local = [[0,0],[0,0]] 
            for i in range(1,3):
                for k in range(1,3):
                    local[i-1][k-1] = abs(M[i][k]%M[0][0])
            nondividx = local.index(max(local))
            M[0] = list(map(operator.add,M[0],M[nondividx+1]))
            A[0] = list(map(operator.add,A[0],A[nondividx+1]))
            continue

        if j == 1 and M[2][2]%M[1][1] != 0:
            M[1] = list(map(operator.add,M[1],M[2]))
            A[1] = list(map(operator.add,A[1],A[2]))
            continue
        elif j != 1:
            j = 1
            continue
        if j == 1 and (M[2][1] != 0 or M[1][2] != 0): #pragma: no cover
            # I literally could not trigger this section of code in my testse.
            continue
        stop_loop = True

    if M[2][2] < 0:
        for i in range(3):
            M[i][2]=-M[i][2]
            B[i][2]=-B[i][2]

    if np.dot(np.dot(A,HNF),B).any() != np.array(M).any(): #pragma: no cover
        print("END: Transformation matrices didn't work.")
        exit()
    check = [M[0][1],M[0][2],M[1][0],M[1][2],M[2][0],M[2][1]]
    if any(check) != 0: #pragma: no cover
        print("Not diagonal")
        exit()
    if M[1][1] % M[0][0] != 0 or M[2][2] % M[1][1] != 0: #pragma: no cover
        print("SNF conditions not met")
        exit()
        
    return(M,A,B)

def _get_dvector_permutations(par_lat,bas_vecs,LatDim,eps):
    """This routine applies the symmetry operations of the parent lattice
    to the interior points (i.e., the d-set) to see which ones are
    equivalent. Labelings of the lattice points that are contain
    permutations only of labels on equivalant sites are physically
    equivalent and therefore redundant. We use these permutations to
    eliminate those duplicate labelings

    Args:
        par_lat (list): A 2D array containing the parent lattice vectors.
        bas_vacs (list): A 2D array containing the basis vectors for the cell.
        LatDim (int): 2 if a 2D case 3 if 3D.
        eps (float): Finite precesion tolerance.

    Returns:
        drp_list (list): The basis vector permutations.
    """

    from phenum.symmetry import get_spaceGroup, bring_into_cell
    from copy import deepcopy
    
    n_d = len(bas_vecs)
    a_type = [1]*n_d

    bv_copy = deepcopy(bas_vecs)
    (rot,shift) = get_spaceGroup(par_lat,a_type,bv_copy,eps = eps)

    if LatDim==2:
        (rot,shift) = _rm_3D_operations(par_lat,rot,shift,eps)

    n_op = len(rot)
    inv_par_lat = np.linalg.inv(np.transpose(par_lat))
    nl_temp= n_op  # Number of operations that fix the parent

    # lattice (but may permute the d-vectors)
    v_temp = []
    perm_temp = []
    for iOp in range(n_op): # Try each operation in turn and see how the
        # d-vectors are permuted for each
        # Rotate each d and add the shift
        if len(bas_vecs) == 1:
            temp_b = bas_vecs[0]
        else:
            temp_b = bas_vecs
        rd_rot = np.matmul(temp_b,rot[iOp]).tolist()

        if n_d > 1:
            rd_shift = [shift[iOp]]*n_d
            rd =  [[rd_rot[i][j] + rd_shift[i][j] for j in range(len(rd_rot[i]))] for i in range(len(rd_rot))]
            
        else:
            rd = [rd_rot[i] + shift[iOp][i] for i in range(len(rd_rot))]
            
        trd = list(rd)
        if n_d > 1:
            for iD in range(n_d):
                rd[iD] = bring_into_cell(rd[iD],inv_par_lat,np.transpose(par_lat),eps)
        else:
            rd = bring_into_cell(rd,inv_par_lat,np.transpose(par_lat),eps)
            
        # The v vector is the vector that must be added (it's a lattice
        # vector) to move a rotated d-vector back into the parent cell.

        if n_d > 1:
            v_temp.append([[rd[i][j] - trd[i][j] for j in range(len(rd[i]))] for i in range(len(rd))])
        else:
            v_temp.append([[rd[i] - trd[i] for i in range(len(rd))]])
        if n_d > 1:
            perm_temp.append(_map_dvector_permutation(rd,bas_vecs,eps))
        else:
            perm_temp.append(_map_dvector_permutation([rd],[bas_vecs],eps))

    drp_list = RotPermList(nL=nl_temp,v=v_temp,perm=perm_temp)

    # I don't think we should reduce this list to a unique one. Some
    # rotations that don't permute the d's could still permute the
    # g's. So we have to keep all the d permutations, even if they
    # look redundant here.
    return(drp_list)
    
def _get_rotation_perms_lists(A,HNF,L,SNF,Op,RPlist,dperms,eps, arrows=False):
    """For each HNF, we have a list of the operations (rotations + shifts,
    if present) that leave the superlattice fixed. Given this set of
    fixing operations, make a list of the permutations on the
    d-vectors (interior points of the multilattice) effected by the
    rotations. Then sort the HNFs into blocks that have the same
    rotation permutations.

    Args:
        A (numpy ndarray): Lattice vectors of the parent lattice as rows of a matrix.
        HNF (numpy ndarray): The HNF matrix.
        L (numpy ndarray): Left transforms matrix.
        SNF (numpy ndarray): The SNF matrix.
        Op (OpList): A list of symmetry ops (rots and shifts) for the parent multilattice.
        RPlist (list of RotPermList): A list of permutations effected by the Ops.
        dperms (list): The permutation of the basis vectors.
        eps (float): Finite precision tolerance.
        arrows (bool, optional): True if arrow group is to be found as well.

    Returns:
        RPlist (list of RotPermList): The symmetry operations represented as permutations
          of a list.
    """

    # Index of the superlattices; Number of d-vectors in d set    
    n = int(round(np.linalg.det(HNF[0])))
    n_h = len(HNF)    
    n_d = np.array(dperms.perm.site_perm).shape[-1]

    skip = []
    gp = []

    arrowg = [[1,0,0],[-1,0,0],[0,1,0],[0,-1,0],[0,0,1],[0,0,-1]]
    
    tg = []
    perm = []
    ident = []
    tperms_perm = []
    temp_rperms_perm = []
    ident = np.reshape([range(n*n_d)],(n_d,n))

    len_rp = len(RPlist)
    for i in range(len_rp):
        RPlist[i].nL=0  # initialize the number
    # Make the group member list for the first SNF
    diag = [SNF[0][0][0],SNF[0][1][1],SNF[0][2][2]]
    g = np.transpose(_make_member_list(diag))

    a_t = np.transpose(A)
    a_tinv = np.linalg.inv(a_t)

    for iH in range(n_h):
        if iH > 0 and not (SNF[iH][0][0] == SNF[iH-1][0][0] and SNF[iH][1][1] == SNF[iH-1][1][1] and SNF[iH][2][2] == SNF[iH-1][2][2]):
            # Make the transform matrices for taking the g's and rotating them
            diag = [SNF[iH][0][0],SNF[iH][1][1],SNF[iH][2][2]]
            g = np.transpose(_make_member_list(diag))
    
        t_inv = np.matmul(L[iH],a_tinv)
        T = np.linalg.inv(t_inv)
        n_op = len(Op[iH].rot)
        na_op = len(arrowg)
        
        temp_rperms_perm = []
        temp_arrow_perm = []
        for iOp in range(n_op): # For each rotation, find the permutation
            dap = np.zeros(na_op)
            dgp = []
            for i in range(n):
                tt = []
                tt = np.zeros(n_d).tolist()
                dgp.append(tt)

            for iD in range(n_d): # Loop over each row in the (d,g) table
                temp1 = np.transpose([RPlist[iH].v[iOp][iD]]*n)
                temp2 = np.matmul(np.matmul(Op[iH].rot[iOp],T),g)
                rag = np.transpose(np.matmul(arrowg,Op[iH].rot[iOp]))
                rgp = np.matmul(t_inv,(-temp1+temp2))
                temp_gp = np.round(rgp).astype(int)
                temp_ag = np.round(rag).astype(int)
                if not np.allclose(rgp,temp_gp,rtol=0,atol=eps): #pragma: no cover
                    print("Transform left big fractional parts")
                    exit()

                gp = temp_gp
                ag = temp_ag
                temp_diag = np.transpose([diag]*n)
                gp = gp%temp_diag # Mod by each entry of
                # the SNF to bring into group Now that the rotated group
                # is known, find the mapping of the elements between the
                # original group and the permuted group. This is the
                # permutation.
                
                skip = [False]*n
                for im in range(n):
                    for jm in range(n):
                        if skip[jm]:
                            continue # Skip elements whose mapping is already known
                        if np.allclose(gp[:,jm],g[:,im]): # these elements
                            # map to each other The list of operations that fix
                            # the superlattice are a subset of those that fix
                            # the parent lattice. RotIndx stores the indicies
                            # of the parent lattice operations in a list with
                            # as many entries as supercell fixing operations.

                            # dperms%perm stores a list of d-vector
                            # permutations, one permutation (an n_d list) for
                            # each operation in the parent lattice symmetries
                            op_indx_in_supercell_list = RPlist[iH].RotIndx[iOp]
                            row_indx_gtable = dperms.perm.site_perm[op_indx_in_supercell_list][iD]
                            dgp[im][row_indx_gtable] = jm+iD*n
                            skip[jm] = True
                            break
                        
                # do the same thing for the arrows
                skip = [False]*na_op
                for im in range(na_op):
                    for jm in range(na_op):
                        if skip[jm]:
                            continue
                        if (ag[:,jm] == arrowg[im]).all():
                            dap[im] = jm
                            skip[jm] = True
                            break

                if np.count_nonzero(dgp == 0) > 1 or np.count_nonzero(dap == 0) > 1: #pragma: no cover
                    print("(d,g)-->(d',g') mapping failed in get_rotation_perm_lists")
                    exit()

            # Now we have the (d',g') table for this rotation. Now
            # record the permutation
            # permutation in the "long form"
            temp_rperms_perm.append(np.transpose(dgp).reshape(n_d*n).tolist()) # store
            temp_arrow_perm.append(dap)

        # nomenclature:
        # N+t = rotation (N) + fractional translation (t)  (me bethinks....)
        # r = lattice translation

        # Now that we have the permutations that are effected by N+t
        # type of rotations (for each N in the point group), we need to
        # compose them with the permutations effected by lattice
        # translations, (parent) lattice translations contained inside the
        # supercell (N+t+r, in Rod's nomenclature). Only when we have these
        # two lists, and compose them, do we have a guarantee that the
        # resulting list is actually a group. For efficiency, we reduce the
        # N+t list (remove duplicates). Rod claims that the compositions
        # (with the translations) will not have duplicates.
        if len(temp_rperms_perm) > 1 and arrows == False:
            temp_rperms_perm.sort()
            temp_rperms_perm = list(temp_rperms_perm for temp_rperms_perm, _ in itertools.groupby(temp_rperms_perm))
        elif len(temp_rperms_perm) > 1 and arrows == True:
            nr = len(temp_rperms_perm[0])
            na = len(temp_arrow_perm[0])
            len_temp_rp = len(temp_rperms_perm)
            perms = []
            for i in range(len_temp_rp):
                perms.append(list(temp_rperms_perm[i])+list(temp_arrow_perm[i]))

            perms = list(perms for perms, _ in itertools.groupby(perms))
            perms.sort()
            temp_rperms_perm = []
            temp_arrow_perm = []
            len_perms = len(perms)
            for i in range(len_perms):
                temp_rperms_perm =[p[:nr] for p in perms]
                temp_arrow_perm = [p[-na:] for p in perms]
        # The rotations permutations list is now in "alphabetical"
        # order and contains no duplicates
    
        # To get the permutations effected by the r's, we don't need
        # the r's. We can merely take the member list, the group, and
        # add each element of the group to the group itself and see
        # what permutation happens. (This part is somewhat redundant
        # with make_translation_group in labeling_related module.)
        tperms_perm = []
        for ig in range(n): # The number of r's inside the superlattice (the
            # translation perms) is the same as the index n
            temp_g = np.transpose([g[:,ig]]*n)
            tg = g+temp_g#[[g[i][j]+temp_g[i][j] for j in range(len(g[i]))] for i in range(len(g))] # Add the element to the group
            temp_diag = np.transpose([diag]*n)
            tg = tg%temp_diag#[[tg[i][j]%temp_diag[i][j] for j in range(len(tg[i]))] for i in range(len(tg))] # mod by the SNF entries to
            # bring it back to the "primitive" representation
            perm = _find_permutation_of_group(np.transpose(g),np.transpose(tg))
            temp_ident = []

            len_ident = len(ident)
            for il in range(len_ident):
                temp_ident.append([ident[il][i] for i in perm])
            tperms_perm.append(np.reshape(temp_ident,(n*n_d)).tolist())
            
        rp_list_nl = len(temp_rperms_perm)*n
        rp_list_perm_sites = [0]*rp_list_nl
        rp_list_perm_arrows = [0]*rp_list_nl
        len_t_perm = len(temp_rperms_perm)
        for it in range(n): # Loop over unique rotation
            for iOp in range(len_t_perm): # Loop over translation perms (r type)
                # perms (N+t type) Form the permutation effected by
                # composing the iOp-th one with the it-th one
                rp_list_temp = [tperms_perm[it][i] for i in temp_rperms_perm[iOp]]
                rp_list_perm_sites[iOp*n + it] = rp_list_temp
                rp_list_perm_arrows[iOp*n+it] = temp_arrow_perm[iOp]
                # ^--- Having gotten both the rotations and the
                # translation in the sections above (sort_permutations_list
                # etc...), the "operators" of the rotation (one of them is
                # R_i) and the translations (one of them is T_j) are now
                # "scrambled", i.e., (T_j) o (R_i).  Since the *first*
                # Rotation R_1 = Id, the *first* entries in the
                # RPlist(iH)%perm are equivalent to pure translations only.
                # The following entries are a combination of both.

        RPlist[iH].perm.site_perm = rp_list_perm_sites
        if arrows: 
            RPlist[iH].perm.arrow_perm = rp_list_perm_arrows
        else:
            RPlist[iH].perm.arrow_perm = [range(na_op)]*len(rp_list_perm_arrows)
        RPlist[iH].n = rp_list_nl
        
    return (RPlist)
    
def get_sym_group(par_lat,bas_vecs,HNF,LatDim,arrows=True):
    """Generates the symmetry group for a given lattice and HNF

    Args:

        par_lat (numpy ndarray): a 2D array that contains the parrent
          lattice vectors

        bas_vecs (list): A list of the atomic basis vectors for the lattice.
        HNF (numpy ndarray): The HNF matrix that defines the supercell.
        LatDim (int): 2 if a surface case case 3 if bulk.
        arrows (bool, optional): True if the arrow group is to be found as well.

    Returns:
        sym_group (RotPermList): The symmetry operations represented as permutations
          of a list.
    """

    from numpy import linalg, allclose
    from phenum.symmetry import bring_into_cell, get_spaceGroup
    from copy import deepcopy

    eps = 1E-10
    # map any atoms in the basis that aren't within the cell to be in
    # the cell
    par_lat_inv = np.transpose(linalg.inv(par_lat))
    temp_basis = deepcopy(bas_vecs)
    for i, b_vecs in enumerate(bas_vecs):
        # par_lat_inv = linalg.inv(np.transpose(par_lat))
        bas_vecs[i] = bring_into_cell(b_vecs,par_lat_inv,np.transpose(par_lat),eps)
        if not allclose(b_vecs, temp_basis[i], rtol=0, atol=eps): # pragma: no cover
            from phenum.msg import warn
            warn("An atomic basis vector was not inside the unit cell. It has been "
                 "remapped.\n Original vector: {} \n Remapped vector: {}"
                 .format(" ".join([str(p) for p in temp_basis[i]]),
                         " ".join([str(p) for p in b_vecs])))

    par_rp_list = _get_dvector_permutations(par_lat,bas_vecs,LatDim,eps)

    a_type = [1]*len(bas_vecs)
    
    (sgrots,sgshifts) = get_spaceGroup(par_lat,a_type,bas_vecs,eps)
    (fixing_ops,RPList,degeneracy) = _get_sLV_fixing_operations(HNF,par_lat,sgrots,sgshifts,
                                                                par_rp_list,eps)

    (SNF,L,R) = SmithNormalForm(deepcopy(HNF))
    sym_group = _get_rotation_perms_lists(par_lat,[HNF],[L],[SNF],[fixing_ops],[RPList],
                                          par_rp_list,eps,arrows=arrows)

    return(sym_group[0])

#a_group take the set of translations of the lattice and the set of
#rotations of the lattice paired with their effect on the arrows of
#the lattice. It makes a direct product of the transformations and the
#rotations and pairs them with the combination of th arraws. If the
#resultant operation is not in the group it then adds it to the group.
def a_group(trans,rots):
    """This subroutine combines that translations of the lattice with the
    rotatians of the lattice.

    Args:
        trans (list): A 2D array where each row is an translation
          of the lattice.

        rots (list): A 3D array. Each row's first entry is the
          site permutations and the second entry is the arrow
          permutations.

    Returns:
        groupi (list): The full symmetry group.
    """
    groupi = []
    for i in trans:
        for j in rots:
            c = []            
            c.append([j[0][i[l]] for l in range(0,len(i))])
            c.append(j[1])
            if c not in groupi:
                groupi.append(c)

    return(groupi)

def a_group_gen(trans,rots):
    """This subroutine combines that translations of the lattice with the
    rotatians of the lattice for the generators of the system.
    
    Args:
        trans (list): A 2D array where each row is an translation
          of the lattice.

        rots (list): A 3D array. Each row's first entry is the
          site permutations and the second entry is the arrow
          permutations.

    Returns:
        groupi (list): The full symmetry group.
    """
    groupi = []
    for i in trans:
        for j in rots:
            c = []            
            c.append([j[0][i[l]] for l in range(0,len(i))])
            c.append(j[1])
            if c not in groupi:
                groupi.append(c)

    groupi = _group(groupi)
    
    return(groupi)

def _group(gen): # pragma: no cover
    """This subroutine takes the generators of the group then uses them
    to form the entire group.

    Args:
        gen (list): A 2D array of the generators of the group
          where each row is a seperate generator.

    Returns:
        groupi (list): The full symmetry group.
    """
    # q is a counter to track how many group elements we've found
    q = 1
    # the final output group
    groupi = []
    # loop over the generators
    for i in gen:
        # a second set of generators used to ensure that a group
        # action doesn't act an its self
        gen2 = []      
        # if k is not the ith element of the generators add it to the
        # set gen2
        for k in gen:
            if k != i:
                gen2.append(k)
        # loop over all except the ith generator
        for j in gen2:
            # handle the site group and arrow group independently
            trans = [j[0][i[0][l]] for l in range(0,len(i[0]))]
            rots = [j[1][i[1][l]] for l in range(0,len(i[1]))]
            # c is the combination of the trans and rots, if it's not
            # in groupi then it should be added.
            c = [trans,rots]
            if c not in groupi:
                groupi.append(c)

            # for each new group element, c, we should also see if
            # the ith group element acting on it will make another
            # new group element.
            trans = [c[0][i[0][l]] for l in range(0,len(i[0]))]
            rots = [c[1][i[1][l]] for l in range(0,len(i[1]))]
            # d is the group action after the ith element has acted on
            # c, if it's not in groupi then add it.
            d = [trans, rots]
            if d not in groupi:
                groupi.append(d)
            # until d and c are the same again we want to apply the
            # ith group element to d so that we've found all possible
            # group elements for the combination
            while d != c:
                trans = [d[0][i[0][l]] for l in range(0,len(i[0]))]
                rots = [d[1][i[1][l]] for l in range(0,len(i[1]))]
                d = [trans,rots]
                # increment q because we found a new group element
                q += 1
                # every time d is not in group i add it
                if d not in groupi:
                    groupi.append(d)
    # new is a tracker that is 0 while we are still finding new group
    # elements and 1 when no new ones are found
    new = 0

    #Here each generator is applied to the group over and over until no new group
    #elements are found
    while new == 0:
        # a second set to store any additional group elements
        group2 = []
        for i in gen:
            for h in groupi:
                # apply the ith generator to the hth element of groupi
                # from above to search for new group elements
                trans = [h[0][i[0][l]] for l in range(0,len(i[0]))]
                rots = [h[1][i[1][l]] for l in range(0,len(i[1]))]
                
                d = [trans,rots]
                # if do isn't in groupi or in group2 then add it to
                # group2
                if d not in groupi and d not in group2:
                    group2.append(d)
        # add anything that is in group2 to groupi
        if len(group2) > 0:
            groupi.extend(group2)
        else:
        # if no new group elements were found then we'v found the
        # whole group and new is set to 1
            new = 1
        
    return(groupi)

def get_full_HNF(HNF):
    """Turns the reduced HNF to a full HNF matrix.

    Args:
        HNF (list): The 1D array that contains all the lower
          diagonal entries of the HNF matrix.

    Returns:
        full_hnf (list of list): A 2D list containing the HNF matrix.
    """

    if type(HNF).__module__ == np.__name__:
        temp_hnf = HNF.tolist()
    else:
        temp_hnf = HNF
    full_hnf = [[temp_hnf[0],0,0],[temp_hnf[1],temp_hnf[2],0],temp_hnf[3:]]

    return full_hnf

def _rm_3D_operations(aVecs,sgrots,sgshifts,eps):
    """This subroutine removes operations that are 3 dimmensional.

    Args:
        aVecs (list): The 2D array of the primitive real space lattice vectors.
        sgrots (list): The 3D array containing the space group rotations.
        sgshifts (list): The 2D array containing the space group shifts.
        eps (float): Finite precisions tolerance.

    Returns:
        t_sg_rots (list): The 3D array of the surface rotations.
        t_sg_shifts (list): The 2D array of the surface translations.
    
    Raises:
        ValueError: if the input vectors aren't primitive.
    """
  
    if not np.allclose(
            aVecs[0][1:3],0.0,rtol=0, atol=eps) or not np.allclose(aVecs[2][0:2]
                                                                   ,0.0,rtol=0,atol=eps):
        raise ValueError("Error in rm_3d_operations: only allowed for primitive vectors x00,0xx,0xx")
    
    n_rot = len(sgrots)
    irot = 0
    t_sg_rots = []
    t_sg_shifts = []
    for i in range(n_rot):
        if (np.allclose(sgrots[i][0][1:3],0.0,rtol=0,atol=eps) and
            np.allclose([sgrots[i][1][0],sgrots[i][2][0]],0.0,rtol=0,atol=eps)
            and np.allclose(abs(sgrots[i][0][0]),1.0,rtol=0,atol=eps)):
            # this operation is "2D"         
            irot += 1
            t_sg_rots.append(sgrots[i])
            t_sg_shifts.append(sgshifts[i])

    return(t_sg_rots,t_sg_shifts)
