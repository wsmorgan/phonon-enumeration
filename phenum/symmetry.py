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

    Args:
        size (int): The cell size in integer form
        nspecies (int): the number of atomic species in the system
        nB (int): the number of basis vectors being used

        res_concs (bool): a logical that indicates of the concentrations
          are being restricted

        concs (list of int): A 2D integer array that contains the concentration
          ranges for each atom in the system

    Returns:
        c_list (array-like): The allowed concentration ranges.
    """
    eps = 1E-10
    
    from itertools import product
    if res_concs == True:
        denom = concs[0][2]
        vol_table = []
        for atom in concs:
            minc = min(atom[0:2])
            maxc = max(atom[0:2])
            vol_table.append([int(math.floor(float(minc)/denom*size*nB)),int(math.ceil(float(maxc)/denom*size*nB)),size*nB])

        n = vol_table[0][2]
        digit = [vol_table[i][1]-vol_table[i][0] for i in range(len(vol_table))]
        dig_cnt = [0*i for i in range(len(vol_table))]
        k = len(vol_table)

        label = []
        minv = []
        maxv = []
        for i in range(k):
            label.append(list(range(vol_table[i][0],vol_table[i][1]+1)))
            minv.append(float(min([concs[i][0],concs[i][1]]))/concs[i][2])
            maxv.append(float(max([concs[i][0],concs[i][1]]))/concs[i][2])

        a = [label[i][0] for i in range(len(label))]
        
        cc = 0
        c_list = []
        done = False
        while done == False:
            if sum(a) == n:
                conc = []
                len_a = len(a)
                for i in range(len_a):
                    conc.append(a[i]/float(n))
                if not ((any(conc[i] < (minv[i]-eps) for i in range(len(minv)))) or (any(conc[i] > (maxv[i]+eps) for i in range(len(maxv))))):
                    c_list.append(deepcopy(a))
            j = k-1
            done2 = False
            while done2 == False:
                if dig_cnt[j] != digit[j]:
                    done2 = True
                    break
                a[j] = label[j][0]
                dig_cnt[j] = 0
                j -= 1
                if j < 0:
                    done2 = True
                    break
            if j < 0:
                done = True
                break
            dig_cnt[j] += 1
            a[j] = label[j][dig_cnt[j]]
            
    else:
        c_list = []
        crange = list(range(0,(size*nB)+1))
        aranges = []
        for i in range(nspecies):
            aranges.append(crange)
        for p in product(*aranges):
            if sum(p) == size*nB:
                c_list.append(list(p))
    return(c_list)

def _does_mapping_exist(v,this_type,atom_pos,atomType,eps):
    """Checks to see if a mapping exists between the vector v and the
    position of any of the atoms of type "this_type".  If a mapping
    exists, then the logical "mapped" is returned .true., otherwise
    .false.

    Args:
        v (list of float): Array of the position to check mapping for
        this_type (int): Integer that indicates which type of atom that is
          being checked.

        atom_pos (array-like): 2D array of the positions of the basis
          atoms.

        atomType (list of int): Array of integers of the types of atoms in the
          basis.

        eps (float): Epsilon for checking equivalence.

    Returns:
        mapped (bool): True if mapping exists.
    """
    mapped = False
    for i, a_type in enumerate(atomType):
        if a_type == this_type:
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

    Args:
        par_lat (array-like): A 2D integer array that contains the parent 
          lattice vectors

    Returns:
        prim_to_cart (numpy ndarray): The matrix that transforms from lattice to 
          cartesian coordinates.
        
        cart_to_prim (numpy ndarray): The matrix that tranforms form cartesian to
          lattice coordinates.
    """

    prim_to_cart = deepcopy(par_lat)
    cart_to_prim = numpy.linalg.inv(numpy.array(prim_to_cart)).tolist()

    return(prim_to_cart,cart_to_prim)

def bring_into_cell(vec,cart_to_latt,latt_to_cart,eps):
    """This subroutine translates a point back into the unit cell in
    lattice coordinates, the coefficients of the point must all be
    less than one and at least zero.
    
    Args:
        vec (list of int): A 1D integer array of the atom position vector

        cart_to_latt (numpy ndarray): The matrix that transforms from cartesian 
          to lattice coordinates

        latt_to_cart (numpy ndarray): The matrix that transforms form lattice to 
          cartesian coordinates

        eps (float): Finite precision tolerance.

    Returns:
        vec (list of float): The vector brought back into the cell.
    """

    from numpy import matmul, transpose
    # Put the representation of the point into lattice coordinates
    vec = matmul(transpose(cart_to_latt),vec).tolist()
    # counter to catch compiler bug
    c = 0
    maxc = max(math.ceil(abs(max(vec))),math.ceil(abs(min(vec))))*2

    # If a component >= 1, translate by subtracting a lattice vector
    # If a component < 0, translate by adding a lattice vector
    while any(i > 1.0-eps for i in vec) or any(i < 0.0-eps for i in vec):
        c = c +1
        for i, v in enumerate(vec):
            if v >= 1-eps:
                vec[i] -= 1
            elif v < 0-eps:
                vec[i] += 1
        if (c>maxc): #pragma: no cover
            print("ERROR: loop does not end in bring_into_cell. Probably compiler bug.")
            exit()

    # Put the point back into cartesion coordinate representation
    vec = matmul(transpose(latt_to_cart),vec).tolist()
    return vec

def get_lattice_pointGroup(a_vecs, eps=1E-10):
    """This routine returns only the point group of the lattite rather
      than the space group of the given crystal structure.

    Args:
        a_vecs (array-like): The 2D array that contains the parent lattice 
            vectors as row vectors.

        eps (float, optional): Finite precision tolerance

    Returns:
        lattpg_op (array-like): The point group for the lattice in 
            cartesian coordinates.
    """

    inverse_avecs = numpy.linalg.inv(numpy.array(a_vecs)).tolist()
    # Store the norms of the three lattice vectors
    norm_avecs = []
    for i in range(3):
        norm_avecs.append(numpy.linalg.norm(a_vecs[i]).tolist())

    # Decide how many lattice points to look in each direction to get all the 
    # points in a sphere that contains all of the longest _primitive_ vectors
    cell_volume = abs(numpy.dot(a_vecs[0],numpy.cross(a_vecs[1],a_vecs[2])))
    max_norm = max([numpy.linalg.norm(i) for i in a_vecs])
    n1 = math.ceil(max_norm*numpy.linalg.norm(numpy.cross(a_vecs[1],a_vecs[2])/cell_volume+eps))
    n2 = math.ceil(max_norm*numpy.linalg.norm(numpy.cross(a_vecs[2],a_vecs[0])/cell_volume+eps))
    n3 = math.ceil(max_norm*numpy.linalg.norm(numpy.cross(a_vecs[0],a_vecs[1])/cell_volume+eps))

    r_vecs = []
    r_lengths = []

    a_vecs = numpy.array(a_vecs)
    # Store the R vectors that lie within the sphere
    num_rs = 0
    for i in range(-int(round(n1)), int(round(n1))+1):
        for j in range(-int(round(n2)), int(round(n2))+1):
            for k in range(-int(round(n3)), int(round(n3))+1):
                this_vector = i*a_vecs[0] + j*a_vecs[1] + k*a_vecs[2]
                length = numpy.linalg.norm(this_vector)
                if (length > max_norm + eps):
                    continue # This vector is outside sphere
                num_rs += 1
                r_vecs.append(this_vector.tolist())
                r_lengths.append(length)
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
    from itertools import permutations
    for i,j,k in permutations(range(num_rs),3):
        if (abs(r_lengths[i] - norm_avecs[0]) > eps) or (abs(r_lengths[j] - norm_avecs[1]) > eps) or (abs(r_lengths[k] - norm_avecs[2]) > eps) or (abs(cell_volume - abs(numpy.linalg.det([r_vecs[i],r_vecs[j],r_vecs[k]]))) > eps):
            continue
        # Form the new set of "rotated" basis vectors
        new_vectors = [r_vecs[i],r_vecs[j],r_vecs[k]]
        # If the transformation matrix that takes the original set to the new set is
        # an orthogonal matrix then this rotation is a point symmetry of the lattice.
        rotation_matrix = numpy.matmul(inverse_avecs,new_vectors).tolist()
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

    Args:
        par_lat (array-like): A 2D integere array that contains the parent lattice vectors
        atomType (list of int): Integer array representing the type of each basis atom
        bas_vecs (array-like): A 2D integere array that contains the basis vectors for the cell
        eps (float, optional): Finite precisions tolerance

        lattcoords (bool, optional): True if vectors are in lattice coordinates 
          rather than cartesian

    Returns:
        sg_ops (array-like): The rotation and mirror operations of the space group.
        sg_fracts (array-like): The translation operations of the space group.
    """
    # Get number of atoms in the basis    
    n_atoms = len(atomType)

    # save original atomic input positions
    atom_pos = deepcopy(bas_vecs)

    # A vector can be represented as either a set of cartesian coordi-
    # nates or as a linear combination of primitive lattice vectors
    # Get transformation matrices to take us back and forth
    (latt_to_cart,cart_to_latt) = _get_transformations(par_lat)
    
    # If we're in lattice coordinates Convert the position of the
    # basis atoms from lattice coordinates.
    if lattcoords:
        for i in range(n_atoms):
            atom_pos[i] = numpy.matmul(latt_to_cart,atom_pos[i]).tolist()

    # bring all the basis atoms into the unit cell
    for i, a_pos in enumerate(atom_pos):
        atom_pos[i] = bring_into_cell(a_pos,cart_to_latt,latt_to_cart,eps)

    # Now find the point group
    lattpg_op = get_lattice_pointGroup(par_lat,eps=eps)

    # **** Find the elements of the space group ****
    # Count the elements
    sgop_count = 0
    sg_ops = []
    sg_fracts = []
    # Apply each of the point operators in the point group to the crystal
    for op in lattpg_op:
        # rotate atom 1 and store its position in the vector v
        v = numpy.matmul(op,atom_pos[0])
        # Loop over all possible fractional translations
        for jAtom in range(n_atoms):
            if (atomType[jAtom] != atomType[0]):
                continue #pragma: no cover
            fract = [atom_pos[jAtom][i] - v[i] for i in range(3)]
            fract = bring_into_cell(fract, cart_to_latt, latt_to_cart, eps)
            # Is each atom of every type mapped by this rotation + translation?
            for kAtom in range(n_atoms):
                this_type = atomType[kAtom]
                # Rotate and translate each atom        
                v2 = numpy.matmul(op,atom_pos[kAtom])
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
                sg_ops.append(op) # Store the rotational part
                # loop over fractional translations and try next op
                # By removing the preceding exit, we include fractional translations
                # for non-primitive lattices. (GLWH 10/26/2009)

    return(sg_ops,sg_fracts)
