"""Methods used for generating the symmetry group. All of these
methods were adapted to python from their original fortran
implementations which can be found at:
https://github.com/msg-byu/symlib
"""

from copy import deepcopy
import numpy
import math

def get_concs_for_size(size,nspecies,res_concs,nB,concs):
    """Gets the concentration ranges for the atoms within the cells of
    certain sizes given the constraints provided such as
    concentration restrictions and the presence of arrows. Code
    rewritten from the get_conetration_list subroutine of:
    https://github.com/msg-byu/enumlib/blob/master/src/derivative_structure_generator.f90  

    :arg size: the cell size in integer form
    :arg nspecies: the number of atomic species in the system
    :arg res_concs: a logical that indicates of the concentrations
    are being restricted
    :arg nB: the number of basis vectors being used
    :arg concs: an 2D integer array that contains the concentration
    ranges for each atom in the system
    """
    eps = 1E-10
    
    from itertools import product
    if res_concs == True:
        denom = concs[0][2]
        volTable = []
        for atom in concs:
            minc = min(atom[0:2])
            maxc = max(atom[0:2])
            volTable.append([int(math.floor(float(minc)/denom*size*nB)),int(math.ceil(float(maxc)/denom*size*nB)),size*nB])

        n = volTable[0][2]
        digit = [volTable[i][1]-volTable[i][0] for i in range(len(volTable))]
        digCnt = [0*i for i in range(len(volTable))]
        k = len(volTable)

        label = []
        minv = []
        maxv = []
        for i in range(k):
            label.append(list(range(volTable[i][0],volTable[i][1]+1)))
            minv.append(float(min([concs[i][0],concs[i][1]]))/concs[i][2])
            maxv.append(float(max([concs[i][0],concs[i][1]]))/concs[i][2])

        a = [label[i][0] for i in range(len(label))]
        
        cc = 0
        cList = []
        done = False
        while done == False:
            if sum(a) == n:
                conc = []
                for i in range(len(a)):
                    conc.append(a[i]/float(n))
                if not ((any(conc[i] < (minv[i]-eps) for i in range(len(minv)))) or (any(conc[i] > (maxv[i]+eps) for i in range(len(maxv))))):
                    cList.append(deepcopy(a))
            j = k-1
            done2 = False
            while done2 == False:
                if digCnt[j] != digit[j]:
                    done2 = True
                    break
                a[j] = label[j][0]
                digCnt[j] = 0
                j -= 1
                if j < 0:
                    done2 = True
                    break
            if j < 0:
                done = True
                break
            digCnt[j] += 1
            a[j] = label[j][digCnt[j]]
            
    else:
        cList = []
        crange = list(range(0,size+1))
        aranges = []
        for i in range(nspecies):
            aranges.append(crange)
        for p in product(*aranges):
            if sum(p) == size:
                # if not any([list(c) in cList for c in it.permutations(p)]) == True:
                cList.append(list(p))
    return(cList)

def _does_mapping_exist(v,this_type,atom_pos,atomType,eps):
    """Checks to see if a mapping exists between the vector v and the
      position of any of the atoms of type "this_type".  If a mapping
      exists, then the logical "mapped" is returned .true., otherwise
      .false.

      :args v: Array of the position to check mapping for
      :args this_type: Integer that indicates which type of atom that is
            being checked.
      :args atom_pos: 2D array of the positions of the basis
            atoms.
      :args atomType: Array of integers of the types of atoms in the
            basis.
      :args eps: Epsilon for checking equivalence.
    """
    mapped = False
    for i in range(len(atomType)):
        if atomType[i] == this_type:
            # if the coordinates are the same, 
            # their difference will be zero for every component
            this_position = atom_pos[i]
            if(numpy.allclose(numpy.array(v), numpy.array(this_position), rtol=0,atol=eps)):
                mapped = True
                break
    return mapped


def _get_transformations(par_lat):

    """This routine generates the matrices for converting vectors from
      lattice coordinates to cartesion coordinates and vice versa.

      :args par_lat: A 2D integer array that contains the parent lattice vectors
    """

    prim_to_cart = deepcopy(par_lat)
    cart_to_prim = numpy.linalg.inv(numpy.array(prim_to_cart)).tolist()

    return(prim_to_cart,cart_to_prim)

def bring_into_cell(vec,cart_to_latt,latt_to_cart,eps):
    """This subroutine translates a point back into the unit cell in
      lattice coordinates, the coefficients of the point must all be
      less than one and at least zero.
    
      :args vec: A 1D integer array of the atom position vector
      :args cart_to_latt: The matrix that transforms from cartesian to lattice coordinates
      :args latt_to_cart: The matrix that transforms form lattice to cartesian coordinates
      :args eps: Finite precision tolerance.
    """

    from numpy import matmul
    # Put the representation of the point into lattice coordinates
    vec = matmul(vec,cart_to_latt).tolist()
    # counter to catch compiler bug
    c = 0
    maxc = max(math.ceil(abs(max(vec))),math.ceil(abs(min(vec))))*2

    # If a component >= 1, translate by subtracting a lattice vector
    # If a component < 0, translate by adding a lattice vector
    while any([i > 1.0-eps for i in vec]) or any([i < 0.0-eps for i in vec]):
        c = c +1
        for i in range(len(vec)):
            if vec[i] >= 1-eps:
                vec[i] -= 1
            elif vec[i] < 0-eps:
                vec[i] += 1
        if (c>maxc): #pragma: no cover
            print("ERROR: loop does not end in bring_into_cell. Probably compiler bug.")
            exit()

    # Put the point back into cartesion coordinate representation
    vec = matmul(vec,latt_to_cart).tolist()
    return vec

def _get_lattice_pointGroup(aVecs, eps=1E-10):
    """This routine returns only the point group of the rather than the
      space group of the given crystal structure.

      :args aVecs: The 2D integer array that contains the parent lattice vectors
      :args eps: (Optional) Finite precision tolerance
    """

    inverse_aVecs = numpy.linalg.inv(numpy.array(aVecs)).tolist()
    # Store the norms of the three lattice vectors
    norm_avecs = []
    for i in range(3):
        norm_avecs.append(numpy.linalg.norm(aVecs[i]).tolist())

    # Decide how many lattice points to look in each direction to get all the 
    # points in a sphere that contains all of the longest _primitive_ vectors
    cell_volume = abs(numpy.dot(aVecs[0],numpy.cross(aVecs[1],aVecs[2])))
    max_norm = max([numpy.linalg.norm(i) for i in aVecs])
    n1 = math.ceil(max_norm*numpy.linalg.norm(numpy.cross(aVecs[1],aVecs[2])/cell_volume+eps))
    n2 = math.ceil(max_norm*numpy.linalg.norm(numpy.cross(aVecs[2],aVecs[0])/cell_volume+eps))
    n3 = math.ceil(max_norm*numpy.linalg.norm(numpy.cross(aVecs[0],aVecs[1])/cell_volume+eps))

    Rvecs = []
    Rlengths = []

    aVecs = numpy.array(aVecs)
    # Store the R vectors that lie within the sphere
    num_Rs = 0
    for i in range(-int(round(n1)), int(round(n1))+1):
        for j in range(-int(round(n2)), int(round(n2))+1):
            for k in range(-int(round(n3)), int(round(n3))+1):
                this_vector = i*aVecs[0] + j*aVecs[1] + k*aVecs[2]
                length = numpy.linalg.norm(this_vector)
                if (length > max_norm + eps):
                    continue # This vector is outside sphere
                num_Rs += 1
                Rvecs.append(this_vector.tolist())
                Rlengths.append(length)
    # Try all R vector triplets in the sphere and see which ones are valid 
    # rotations of the original basis vectors.
    # 
    # The length of all vectors must be preserved under a unitary
    # transformation so skip any trial vectors that aren't the same
    # length as the original. We also skip any set of vectors that
    # have the right lengths but do not form a parallelpiped that has
    # the same volume as the original set. Also, note that the we skip
    # sets of vectors that contain the same vector more than once
    # (i.e. the indices i, j, k must be unique).
    num_ops = 0
    lattpg_op = []
    for i in range(num_Rs):
        if (abs(Rlengths[i] - norm_avecs[0]) > eps):
            continue
        for j in range(num_Rs):
            if (abs(Rlengths[j] - norm_avecs[1]) > eps):
                continue
            if (j == i):
                continue
            for k in range(num_Rs):
                if (abs(Rlengths[k] - norm_avecs[2]) > eps):
                    continue
                if (k == i or k == j):
                    continue
                if (abs(cell_volume - abs(numpy.linalg.det([Rvecs[i],Rvecs[j],Rvecs[k]]))) > eps):
                    continue
                # Form the new set of "rotated" basis vectors
                new_vectors = [Rvecs[i],Rvecs[j],Rvecs[k]]
                # If the transformation matrix that takes the original set to the new set is
                # an orthogonal matrix then this rotation is a point symmetry of the lattice.
                rotation_matrix = numpy.matmul(inverse_aVecs,new_vectors).tolist()
                # Check orthogonality of rotation matrix by [R][R]^T = [1]
                test_matrix = numpy.matmul(rotation_matrix,numpy.transpose(rotation_matrix)).tolist()
                if (numpy.allclose(test_matrix, [[1,0,0],[0,1,0],[0,0,1]], rtol=0,atol=eps)): # Found valid rotation
                    num_ops +=  1 # Count number of rotations
                    lattpg_op.append(rotation_matrix)

    return(lattpg_op)
    
def get_spaceGroup(par_lat,atomType,bas_vecs,eps=1E-10,lattcoords = False):

    """This routine takes a crystal structure (basis vectors and basis
      atoms and positions) and returns the point operators and
      fractional translations of the space group. The routine assumes
      that the given crystal structure is already primitive.

      :args par_lat: A 2D integere array that contains the parent lattice vectors
      :args atomType: Integer array representing the type of each basis atom
      :args bas_vecs: A 2D integere array that contains the basis vectors for the cell
      :args eps: (Optional) Finite precisions tolerance
      :args lattcoords: (Optional) True if vectors are in lattice
            coordinates rather than cartesian
    """
    # Get number of atoms in the basis    
    nAtoms = len(atomType)

    # save original atomic input positions
    atom_pos = deepcopy(bas_vecs)

    # A vector can be represented as either a set of cartesian coordi-
    # nates or as a linear combination of primitive lattice vectors
    # Get transformation matrices to take us back and forth
    (latt_to_cart,cart_to_latt) = _get_transformations(par_lat)
    
    # If we're in lattice coordinates Convert the position of the
    # basis atoms from lattice coordinates.
    if lattcoords:
        for i in range(nAtoms):
            atom_pos[i] = numpy.matmul(latt_to_cart,atom_pos[i]).tolist()

    # bring all the basis atoms into the unit cell
    for i in range(len(atom_pos)):
        atom_pos[i] = bring_into_cell(atom_pos[i],cart_to_latt,latt_to_cart,eps)

    # Now find the point group
    lattpg_op = _get_lattice_pointGroup(par_lat,eps=eps)

    # **** Find the elements of the space group ****
    # Count the elements
    sgop_count = 0
    sg_ops = []
    sg_fracts = []
    # Apply each of the point operators in the point group to the crystal
    for iop in range(len(lattpg_op)):
        # rotate atom 1 and store its position in the vector v
        v = numpy.matmul(lattpg_op[iop],atom_pos[0]).tolist()
        # Loop over all possible fractional translations
        for jAtom in range(nAtoms):
            if (atomType[jAtom] != atomType[0]):
                continue #pragma: no cover
            fract = [atom_pos[jAtom][i] - v[i] for i in range(3)]
            fract = bring_into_cell(fract, cart_to_latt, latt_to_cart, eps)
            # Is each atom of every type mapped by this rotation + translation?
            for kAtom in range(nAtoms):
                this_type = atomType[kAtom]
                # Rotate and translate each atom        
                v2 = numpy.matmul(lattpg_op[iop],atom_pos[kAtom]).tolist()
                v2 = [v2[i] + fract[i] for i in range(3)]
                v2 = bring_into_cell(v2, cart_to_latt, latt_to_cart, eps)
                
                # Try to map this rotated atom onto another the same type
                mapped = _does_mapping_exist(v2, this_type, atom_pos, atomType, eps)
                if not mapped:
                    break # no mapping for this atom

            # if all atoms mapped successfully then count this
            # element (rotation + translation)
            if mapped:
                sgop_count += 1 # Count the number of elements
                sg_fracts.append(fract) # Save the translational part
                sg_ops.append(lattpg_op[iop]) # Store the rotational part
                # loop over fractional translations and try next op
                # By removing the preceding exit, we include fractional translations
                # for non-primitive lattices. (GLWH 10/26/2009)

    return(sg_ops,sg_fracts)
