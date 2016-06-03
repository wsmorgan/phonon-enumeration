"""Methods for taking an HNF and finding the site permutations and the
arrow permutations. All of these have been modified for python from
their original fortran implementations which can be found at:
https://github.com/msg-byu/enumlib/. Of these only get_rotation
get_rotation_perms_lists() has been modified in order to produce the
arrow group. The class ArrowPerm has also been added and implemented
inside the RotPermList class.
"""

import itertools
import numpy
from copy import deepcopy
import operator
import symmetry as sym 

class ArrowPerm(object):
    """ArrowPerm pairs a site permutation with an arrow permutation."""
    def __init__(self,site_perm = None,arrow_perm = None):
        """Initializes the ArrowPerm.
        :args site_perm: A 1D integer array containing a site permutation.
        :args arrow_perm: A 1D integer array containing an arrow permutation.
        """
        self.site_perm = site_perm
        self.arrow_perm = arrow_perm

class RotPermList(object):
    """RotPermList maintains the data structure for the rotations and
      permutations from the fortran code.
    """
    def __init__(self,nL = None,v =None,perm = None,RotIndx = None,arrows=None):
        """Initializes the RotPermList.
          :args nL: An integer indicating the number of operations
          :args v: A 2D integer array containing the lattice vectors
          :args perm: 
        """
        self.nL = nL
        self.v = v
        self.perm = ArrowPerm(site_perm = perm, arrow_perm=arrows)
        self.RotIndx = RotIndx

class opList(object):
    """opList maintains the data structure of the fixing operations form
      the fortran code.
    """
    def __init__(self,rot = None, shift = None):
        """Initializes the opList.
          :args v: A 3D array containing the rotations
          :args perm: An 2D array containing the shifts
        """
        self.rot = rot
        self.shift = shift

def _make_member_list(n):
    """Takes the length of three cyclic group and constructs the member
    list so that the permutations can be determined. Each member has
    three components, corresponding to the entries for each of the
    three cyclic groups.
    
    :args n: Integer array containing the diagonal elements of the SNF.
    """
    from operator import mul
    depth = reduce(mul,n,1)
    p = []
    for i in range(depth):
        p.append([0,0,0])
    for im in range(1,reduce(mul,n,1)):  # Loop over the members of the translation group
        p[im] = deepcopy(p[im-1]) # Start with the same digits as in the previous increment        
        p[im][2] = (p[im-1][2]+1)%n[2]  # Increment the first cyclic group
        if (p[im][2]==0):             # If it rolled over then
            p[im][1] = (p[im-1][1] +1) % n[1] # increment the next cyclic group
            if (p[im][1]==0):          # If this one rolled over too
                p[im][0] = (p[im-1][0]+1)%n[0] # Then increment the third one
    return p

def _find_permutation_of_group(g,gp):
    """
      :args g: Unpermuted groups
      :args gp: Permuted groups.
      :args perm: Permuted of gp.
    """
    n = len(g)
    perm = []
    skip = []
    for i in range(n):
        skip.append(False) # This is just for efficiency
    for im in range(n):
        for jm in range(n):
            if skip[jm]:
                continue # This is just for efficiency
            if gp[jm]==g[im]:
                perm.append(jm)
                skip[jm] = True
                break # don't keep looking if you already found the match
            
    return perm

def _is_equiv_lattice(lat1,lat2,eps):

    """This function determines whether two lattices are equivalent. That
      is, if one is an equal-volume derivative lattice of the
      other. It uses the idea of Santoro and Mighell
      (Acta. Cryst. 1972) that, for equivalent lattices, the
      transformation matrix that takes the vectors of str2 to str1 has
      determinant of 1 and all integer
      elements. [LV2]=[LV1][S]=&gt;S=LV1^-1*LV2. But I had to use abs
      on the determinant. This is not mentioned by Santoro and Mighell
      (issue of opposite handedness---which doesn't matter to me).

      :args lat1: 2D integer array containing the first lattice.
      :args lat2: 2D integer array containing the second lattice.
      :args eps: The finite precision tolerance.
    """
    is_equiv_lattice = False
    lat1inv = numpy.linalg.inv(numpy.array(lat1)).tolist()
    S = numpy.matmul(numpy.array(lat2),lat1inv).tolist()
    if (numpy.allclose(abs(numpy.linalg.det(numpy.array(S))),1.0,rtol=eps,atol=eps) and
        numpy.allclose(numpy.array(S),numpy.array([[int(round(S[i][j])) for j in range(3)] for i in range(3)]),rtol=eps,atol=eps)):
        is_equiv_lattice = True
    return is_equiv_lattice
        
def _get_sLV_fixing_operations(HNF,pLV,nD,rot,shift,dPerm,eps):
    """
      :args HNF:
      :args pLV:
      :args nD:
      :args rot:
      :args shift:
      :args dPerm:
      :args eps: Finite precision tolerance.
    """
    
    nRot = len(rot)

    degen_lattices = []
    for i in range(nRot):
        degen_lattices.append([[0,0,0],[0,0,0],[0,0,0]])
        
    cDegen = 0
    ic = 0 # Counter for the fixing operations
    tv = 0
    tmpOp_rot = []
    tmpOp_shift = []
    tv = []
    tIndex = [] # temp variables
    for iRot in range(nRot):  # Loop over each rotation
        thisRot = rot[iRot] # Store the rotation
        origLat = numpy.matmul(HNF,pLV).tolist()  # Compute the superlattice
        rotLat = numpy.matmul(numpy.array(origLat),numpy.array(thisRot)).tolist() # Compute the rotated superlattice
        if _is_equiv_lattice(rotLat,origLat,eps):
            # this operation fixes the lattice and should be recorded
            ic += 1
            tmpOp_rot.append(thisRot)
            tmpOp_shift.append(shift[iRot])
            tv.append(dPerm.v[iRot])
            tIndex.append(iRot)
            # Added by LN from here
        else:
            inList = False
            for iDegen in range(cDegen):
                if _is_equiv_lattice(degen_lattices[iDegen],rotLat,eps):
                    inList = True
                    break
          
            if not inList:
                degen_lattices[cDegen] = rotLat
                cDegen += 1
        # End of LN additions for degeneracy
        # Now we know which rotations fix the lattice and how many
        # there are so store them
    degeneracy = cDegen
    # Allocate the storage for them
    fixOp = opList(tmpOp_rot,tmpOp_shift) # Stuff the rotations into the permanent array
    if nD > 1:
        rotPerm = RotPermList(v=tv,RotIndx=tIndex)
    else:
        rotPerm = RotPermList(v=[tv],RotIndx=tIndex)

    return(fixOp,rotPerm,degeneracy)

def _map_dvector_permutation(rd,d,eps,n):
    """
      :args rd: 2D array of rotated basis vectors
      :args d: 2D array of original basis vectors
      :args eps: Finite precision tolerance.
      :args n: number of basis vectors
    """

    found = []
    for i in range(len(rd)):
        found.append(False)
    nD = len(rd) # of d-vectors
    RP = []
    for iD in range(nD):
        for jD in range(nD):
            if found[jD]:
                continue
            if numpy.allclose(numpy.array(rd[iD]),numpy.array(d[jD]),atol=eps,rtol=eps):
                RP.append(jD)
                found[jD] = True
                break

    if len(RP) != len(d):
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

      :args invec: The input integer array.
    """

    vec = [abs(i) for i in invec]
    
    this_min = invec.index(min([i for i in vec if i > 0]))

    rvec = [i for i in reversed(vec)]
    
    this_max = 2 - rvec.index(max(rvec))

    return (this_min,this_max)

def SmithNormalForm(HNF):
    """Finds the Smith Normal Form (SNF) of an HNF as well as the left and
      right transforms.
      
      :args HNF: The integer matrix HNF for which the SNFs are to be
            found
    """

    if numpy.linalg.det(numpy.array(HNF)) < 1:
        print("SmithNormalForm routine failed because the input matrix had a negative determinant")
        exit()
    
    A = [[0,0,0],[0,0,0],[0,0,0]]
    M = deepcopy(HNF)
    B = [[0,0,0],[0,0,0],[0,0,0]]

    for i in range(3):
        A[i][i] = 1
        B[i][i] = 1

    j = 0
    itCnt = 0

    stop_loop = False
    while not stop_loop:
        itCnt += 1
        if (itCnt >=100):
            print("Bad progromming in SmithNormalForm")
            exit()

        while (3-[M[0][j],M[1][j],M[2][j]].count(0)) > 1:
            (minidx,maxidx) = _find_minmax_indices([M[0][j],M[1][j],M[2][j]])
            minm = M[minidx][j]
            mult = M[maxidx][j]/minm

            M[maxidx] = map(operator.sub,M[maxidx],[mult*i for i in M[minidx]])
            A[maxidx] = map(operator.sub,A[maxidx],[mult*i for i in A[minidx]])

            if numpy.dot(numpy.dot(numpy.array(A),numpy.array(HNF)),numpy.array(B)).any() != numpy.array(M).any():
                print("ROW: Transformation matrices didn't work")
                exit()

        if M[j][j] == 0:
            maxidx = [abs(M[0][j]),abs(M[1][j]),abs(M[2][j])].index(max([abs(M[0][j]),abs(M[1][j]),abs(M[2][j])]))
            tmprow = deepcopy(A[j])
            A[j] = deepcopy(A[maxidx])
            A[maxidx] = tmprow

            tmprow = deepcopy(M[j])
            M[j] = deepcopy(M[maxidx])
            M[maxidx] = tmprow    
        if numpy.dot(numpy.dot(numpy.array(A),numpy.array(HNF)),numpy.array(B)).any() != numpy.array(M).any():
            print("ROWSWAP: Transformation matrices didn't work")
            exit()

        while (3-M[j].count(0)) >1:
            (minidx,maxidx) = _find_minmax_indices(M[j])
            minm = M[j][minidx]
            mult = M[j][maxidx]/M[j][minidx]

            for i in range(3):
                M[i][maxidx] = M[i][maxidx]-mult * M[i][minidx]
                B[i][maxidx] = B[i][maxidx]-mult * B[i][minidx]

            if numpy.dot(numpy.dot(numpy.array(A),numpy.array(HNF)),numpy.array(B)).any() != numpy.array(M).any():
                print("COLS: Transformation matrices didn't work")
                exit()

            if M[j][j] < 0:
                for i in range(3):
                    M[i][j] = -M[i][j]
                    B[i][j] = -B[i][j]
            else:
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
        if numpy.dot(numpy.dot(numpy.array(A),numpy.array(HNF)),numpy.array(B)).any() != numpy.array(M).any():
            print("COLSWAP: Transformation matrices didn't work")
            exit()

        Ldiv = [[True,True],[True,True]]
        for i in range(1,3):
            for k in range(1,3):
                Ldiv[i-1][k-1] = (M[i][k]%M[0][0] == 0)

        if j == 0 and any(Ldiv) == False:
            local = [[0,0],[0,0]]
            for i in range(1,3):
                for k in range(1,3):
                    local[i][k] = abs(M[i][k]%M[0][0])
            nondividx = local.index(max(local))
            M[0] = map(operator.add,M[0],M[nondividx+1])
            A[0] = map(operator.add,A[0],A[nondividx+1])
            continue

        if j == 1:
            if M[2][2]%M[1][1] != 0:
                M[1] = map(operator.add,M[1],M[2])
                A[1] = map(operator.add,A[1],A[2])
                continue
        else:
            j = 1
            continue
        if j == 1 and (M[2][1] != 0 or M[1][2] != 0):
            continue
        stop_loop = True
         
    if M[2][2] < 0:
        for i in range(3):
            M[i][2]=-M[i][2]
            B[i][2]=-B[i][2]
    
    if numpy.dot(numpy.dot(numpy.array(A),numpy.array(HNF)),numpy.array(B)).any() != numpy.array(M).any():
        print("END: Transformation matrices didn't work.")
        exit()
    check = [M[0][1],M[0][2],M[1][0],M[1][2],M[2][0],M[2][1]]
    if any(check) != 0:
        print("Not diagonal")
        exit()
    if M[1][1] % M[0][0] != 0 or M[2][2] % M[1][1] != 0:
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

      :args pLV: A 2D integer array containing the parent lattice vectors
      :args bas_vacs: A 2D integer array containing the basis vectors for the cell
      :args LatDim: 2 if a 2D case 3 if 3D
      :args eps: finite precesion tolerance
    """

    nD = len(bas_vecs)
    aTyp = []
    for i in range(nD):
        aTyp.append(1)

    (rot,shift)= sym.get_spaceGroup(par_lat,aTyp,bas_vecs)

    if LatDim==2:
        (rot,shift) = _rm_3D_operations(par_lat,rot,shift,eps)

    nOp = len(rot)
    
    inv_par_lat = numpy.linalg.inv(numpy.array(par_lat)).tolist()
    
    nL_temp= nOp  # Number of operations that fix the parent

    # lattice (but may permute the d-vectors)
    v_temp = []
    perm_temp = []
    for iOp in range(nOp): # Try each operation in turn and see how the
        # d-vectors are permuted for each
        # Rotate each d and add the shift
        if len(bas_vecs) == 1:
            temp_b = bas_vecs[0]
        else:
            temp_b = bas_vecs
        rd_rot = numpy.matmul(temp_b,rot[iOp]).tolist()
        if nD > 1:
            rd_shift = [shift[iOp]]*nD
            rd =  [[rd_rot[i][j] + rd_shift[i][j] for j in range(len(rd_rot[i]))] for i in range(len(rd_rot))]
            
        else:
            rd = [rd_rot[i] + shift[iOp][i] for i in range(len(rd_rot))]
        tRD = deepcopy(rd)
        if nD > 1:
            for iD in range(nD):
                rd[iD] = sym.bring_into_cell(rd[iD],inv_par_lat,par_lat,eps)
        else:
            rd = sym.bring_into_cell(rd,inv_par_lat,par_lat,eps)

        # The v vector is the vector that must be added (it's a lattice
        # vector) to move a rotated d-vector back into the parent cell.
        if nD > 1:
            v_temp.append([[rd[i][j] - tRD[i][j] for j in range(len(rd[i]))] for i in range(len(rd))])
        else:
            v_temp.append([rd[i] - tRD[i] for i in range(len(rd))])
        if nD > 1:
            perm_temp.append(_map_dvector_permutation(rd,bas_vecs,eps,nD))
        else:
            perm_temp.append(_map_dvector_permutation([rd],[bas_vecs],eps,nD))

    dRPList = RotPermList(nL=nL_temp,v=v_temp,perm=perm_temp)

    # I don't think we should reduce this list to a unique one. Some
    # rotations that don't permute the d's could still permute the
    # g's. So we have to keep all the d permutations, even if they
    # look redundant here.

    return(dRPList)
    
def _get_rotation_perms_lists(A,HNF,L,SNF,Op,RPlist,dperms,eps):
    """For each HNF, we have a list of the operations (rotations + shifts,
      if present) that leave the superlattice fixed. Given this set of
      fixing operations, make a list of the permutations on the
      d-vectors (interior points of the multilattice) effected by the
      rotations. Then sort the HNFs into blocks that have the same
      rotation permutations.

      :args A: Lattice vectors of the primary lattice (parent lattice).
      :args HNF: HNF matrix.
      :args L: Left transforms matrix.
      :args SNF: SNF matrix.
      :args Op: A list of symmetry ops (rots and shifts) for the parent multilattice.
      :args RlList:A list of permutations effected by the Ops.
      :args dperms:
      :args eps: Finite precision tolerance.d
    """

    # Index of the superlattices; Number of d-vectors in d set
    n = int(round(numpy.linalg.det(numpy.array(HNF))))

    nD = numpy.array(dperms.perm.site_perm).shape[-1]
    # nD = len(RPlist.v[0])

    skip = []
    gp = []

    arrowg = [[0,0,-1],[0,-1,0],[-1,0,0],[1,0,0],[0,1,0],[0,0,1]]
    
    tg = []
    perm = []
    ident = []
    identT = []
    tperms_perm = []
    temp_rperms_perm = []
    identT = [list(i) for i in zip(*[iter(range(0,n*nD))]*n)]  # we could combine
    # those two lines, but gfortran
    ident  = numpy.transpose(identT).tolist()            # does some strange things then.
    RPlist.nL=0  # initialize the number

    # Make the group member list for the first SNF
    diag = [SNF[0][0],SNF[1][1],SNF[2][2]]
    g = _make_member_list(diag)
    # arrowg = make_arrow_list(diag)
    Ainv = numpy.linalg.inv(A)
    L = numpy.transpose(L).tolist()
    # Make the transform matrices for taking the g's and rotating them
    Tinv = numpy.matmul(L,Ainv).tolist()
    T = numpy.linalg.inv(Tinv).tolist()
       
    nOp = len(Op.rot)

    naOp = len(arrowg)
    
    temp_rperm_perm = []
    temp_arrow_perm = []

    for iOp in range(nOp): # For each rotation, find the permutation
        dap = []
        dgp = []
        for i in range(n):
            tt = []
            for j in range(nD):
                tt.append(0)
            dgp.append(tt)

        for i in range(naOp):
            dap.append(0)

        for iD in range(nD): # Loop over each row in the (d,g) table
            if nD == 1:
                temp1 = [RPlist.v[iD][iOp]]*n
            else:
                temp1 = [RPlist.v[iOp][iD]]*n

            temp2 = numpy.transpose(numpy.matmul(numpy.matmul(Op.rot[iOp],T),numpy.transpose(g))).tolist()
            # temp2a = numpy.transpose(numpy.matmul(numpy.matmul(Op.rot[iOp],T),numpy.transpose(arrowg))).tolist()
            rag = numpy.transpose(numpy.matmul(Op.rot[iOp],numpy.transpose(arrowg))).tolist()
            rgp = numpy.transpose(numpy.matmul(Tinv,numpy.transpose([[temp1[i][j]+temp2[i][j] for j in range(len(temp1[i]))] for i in range(len(temp1))]))).tolist()
            # rag = numpy.transpose(numpy.matmul(Tinv,numpy.transpose([[temp1a[i][j]+temp2a[i][j] for j in range(len(temp1a[i]))] for i in range(len(temp1a))]))).tolist()
            temp_gp = [[int(round(rgp[i][j])) for j in range(len(rgp[i]))] for i in range(len(rgp))] # Move the rotated group into an integer array
            temp_ag = [[int(round(rag[i][j])) for j in range(len(rag[i]))] for i in range(len(rag))]
            if not numpy.allclose(rgp,temp_gp,rtol=eps,atol=eps):
                print("Transform left big fractional parts")
                exit()

            gp = temp_gp
            ag = temp_ag
            temp_diag = [diag]*n
            gp = [[gp[i][j] % temp_diag[i][j] for j in range(len(gp[i]))] for i in range(len(gp))] # Mod by each entry of
            # the SNF to bring into group Now that the rotated group
            # is known, find the mapping of the elements between the
            # original group and the permuted group. This is the
            # permutation.
            skip = []
            for mm in range(n):
                skip.append(False) # This is just for efficiency
            for im in range(n):
                for jm in range(n):
                    if skip[jm]:
                        continue # Skip elements whose mapping is already known
                    if gp[jm] == g[im]: # these elements
                        # map to each other The list of operations that fix
                        # the superlattice are a subset of those that fix
                        # the parent lattice. RotIndx stores the indicies
                        # of the parent lattice operations in a list with
                        # as many entries as supercell fixing operations.

                        # dperms%perm stores a list of d-vector
                        # permutations, one permutation (an nD list) for
                        # each operation in the parent lattice symmetries
                        OpIndxInSuperCellList = RPlist.RotIndx[iOp]
                        RowInDxGTable = numpy.transpose(dperms.perm.site_perm).tolist()[iD][OpIndxInSuperCellList]
                        dgp[im][RowInDxGTable] = jm+iD*n
                        skip[jm] = True
                        break
            # do the some thing for the arrows
            skip = []
            for mm in range(naOp):
                skip.append(False)
            for im in range(naOp):
                for jm in range(naOp):
                    if skip[jm]:
                        continue
                    if ag[jm] == arrowg[im]:
                        dap[im] = jm
                        skip[jm] = True
                        break

            if dgp.count(0) > 1 or dap.count(0) > 1:
                print("(d,g)-->(d',g') mapping failed in get_rotation_perm_lists")
                exit()

        # Now we have the (d',g') table for this rotation. Now
        # record the permutation
        temp_rperms_perm.append(numpy.transpose(dgp).reshape(nD*n).tolist()) # store
        temp_arrow_perm.append(dap)
        # permutation in the "long form"
        temp_rperms_nL = nOp
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
    if len(temp_rperms_perm) > 1 and len(temp_arrow_perm) == 0:
        temp_rperms_perm.sort()
        temp_rperms_perm = list(temp_rperms_perm for temp_rperms_perm, _ in itertools.groupby(temp_rperms_perm))
    # The rotations permutations list is now in "alphabetical"
    # order and contains no duplicates
    
    # To get the permutations effected by the r's, we don't need
    # the r's. We can merely take the member list, the group, and
    # add each element of the group to the group itself and see
    # what permutation happens. (This part is somewhat redundant
    # with make_translation_group in labeling_related module.)
    for ig in range(n): # The number of r's inside the superlattice (the
        # translation perms) is the same as the index n
        temp_g = [g[ig]]*n
        tg = [[g[i][j]+temp_g[i][j] for j in range(len(g[i]))] for i in range(len(g))] # Add the element to the group
        temp_diag = [diag]*n
        tg = [[tg[i][j]%temp_diag[i][j] for j in range(len(tg[i]))] for i in range(len(tg))] # mod by the SNF entries to
        # bring it back to the "primitive" representation
        perm = _find_permutation_of_group(g,tg)
        temp_ident = []
        trans_ident = numpy.transpose(ident).tolist()
        for il in range(len(numpy.transpose(ident).tolist())):
            temp_ident.append([trans_ident[il][i] for i in perm])
        tperms_perm.append(numpy.reshape(temp_ident,(n*nD)).tolist())

    RPlist_perm_sites = []
    RPlist_perm_arrows = []
    RPlist_nL = len(temp_rperms_perm)*n
    for it in range(n): # Loop over translation perms (r type)
        for iOp in range(len(temp_rperms_perm)): # Loop over unique rotation
            # perms (N+t type) Form the permutation effected by
            # composing the iOp-th one with the it-th one
            RPlist_temp = [tperms_perm[it][i] for i in temp_rperms_perm[iOp]]
            RPlist_perm_sites.append(RPlist_temp)
            RPlist_perm_arrows.append(temp_arrow_perm[iOp])
            # ^--- Having gotten both the rotations and the
            # translation in the sections above (sort_permutations_list
            # etc...), the "operators" of the rotation (one of them is
            # R_i) and the translations (one of them is T_j) are now
            # "scrambled", i.e., (T_j) o (R_i).  Since the *first*
            # Rotation R_1 = Id, the *first* entries in the
            # RPlist(iH)%perm are equivalent to pure translations only.
            # The following entries are a combination of both.

    RPlist.perm.site_perm = RPlist_perm_sites
    RPlist.perm.arrow_perm = RPlist_perm_arrows
    RPlist.n = RPlist_nL
    return (RPlist)
    
    
def get_sym_group(par_lat,bas_vecs,HNF,LatDim):
    """Generates the symmetry group for a given lattice and HNF

    :args par_lat: a 2D integer array that contains the parrent
    lattice vectors
    :args bas_vecs: The atomic basis for the lattice
    :args HNF: The HNF that defines the supercell
    :args LatDim: 2 if a 2D case 3 if 3D
    """

    eps = 1E-10
    ParRPList = _get_dvector_permutations(par_lat,bas_vecs,LatDim,eps)

    aTyp = []
    for i in range(len(bas_vecs)):
        aTyp.append(1)
        
    (sgrots,sgshifts) = sym.get_spaceGroup(par_lat,aTyp,bas_vecs,eps)
    (fixing_ops,RPList,degeneracy) = _get_sLV_fixing_operations(HNF,par_lat,len(bas_vecs),sgrots,
                                                    sgshifts,ParRPList,eps)
    
    (SNF,L,R) = SmithNormalForm(HNF)

    symm = _get_rotation_perms_lists(par_lat,HNF,L,SNF,fixing_ops,RPList,ParRPList,eps)

    return(symm)

#a_group take the set of translations of the lattice and the set of
#rotations of the lattice paired with their effect on the arrows of
#the lattice. It makes a direct product of the transformations and the
#rotations and pairs them with the combination of th arraws. If the
#resultant operation is not in the group it then adds it to the group.
def a_group(trans,rots):
    """This subroutine combines that translations of the lattice with the
    rotatians of the lattice.

    :arg trans: a 2D integer array where each row is an translation
    of the lattice
    :arg rots: a 3D integer array. Each row's first entry is the
    site permutations and the second entry is the arrow
    permutations.
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

    :arg trans: a 2D integer array where each row is an translation
    of the lattice
    :arg rots: a 3D integer array. Each row's first entry is the
    site permutations and the second entry is the arrow
    permutations.
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

def _group(gen):
    """This subroutine takes the generators of the group then uses them
    to form the entire group.

    :arg gen: a 2D integer array of the generators of the group
    where each row is a seperate generator
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
                if d not in groupi:
                    if d not in group2:
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
    
    :args HNF: The 1D integer numpy array that contains all the lower
    diagonal entries of the HNF matrix.
    """

    temp_HNF = HNF.tolist()
    full_HNF = [[temp_HNF[0],0,0],[temp_HNF[1],temp_HNF[2],0],temp_HNF[3:]]

    return full_HNF


def _rm_3D_operations(aVecs,sgrots,sgshifts,eps):
    """This subroutine removes operations that are 3 dimmensional.

    :args aVecs: 2D integer array of the primitive real space lattice vectors
    :args sgrots: 3D integer array containing the space group rotations.
    :args sgshifts: 2D integer array containing the space group shifts.
    :args eps: Finite precisions tolerance.
    """
  
    if not numpy.allclose(numpy.array(aVecs[0][1:3]),0.0,rtol=eps,
                          atol=eps) or not numpy.allclose(numpy.array(aVecs[2][0:2])
                                                          ,0.0,rtol=eps,atol=eps):
        print("Error in rm_3d_operations: only allowed for primitive vectors x00,0xx,0xx")
        exit()
    
    nRot = len(sgrots)
    irot = 0
    tSGrots = []
    tSGshifts = []
    for i in range(nRot):
        if (numpy.allclose(numpy.array(sgrots[i][0][1:3]),0.0,rtol=eps,atol=eps) and
            numpy.allclose(numpy.array([sgrots[i][1][0],sgrots[i][2][0]]),0.0,rtol=eps,atol=eps)
            and numpy.allclose(abs(sgrots[i][0][0]),1.0,rtol=eps,atol=eps)):
            # this operation is "2D"         
            irot += 1
            tSGrots.append(sgrots[i])
            tSGshifts.append(sgshifts[i])

    nRot = irot

    return(tSGrots,tSGshifts)
