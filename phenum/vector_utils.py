"""Contains routines used to map the enumerated structures to real space."""

def map_enumStr_to_real_space(system_data,structure_data,minkowskiReduce):
    """Maps an enumerated structure back to real space. Returns a
    dictionary containing the real space data.

    Args:
        system_data (dict): A dictionary containing all the information about the sysytem with
          keys:
          "plattice": The parrent lattice vectors.
          "dvecs": The atomic basis vectors.
          "title": The name of the system.
          "bulksurf": 'bulk' or 'surface'.
          "nD": The number of atomic basis vectors.
          "k": The number of atomic species.
          "eps": Finite precision tolerance.

        sturture_data (dict): A dictionary containing the information for this structure with
          keys:
          "strN": The structure number.
          "hnfN": The HNFs number.
          "hnf_degen": The HNFs degeneracy.
          "lab_degen": The label degeneracy.
          "tot_degen": The total degeneracy for the structure.
          "sizeN": The system size.
          "n": Number of structure within this size.
          "pgOps": The number of point group operations.
          "diag": The diagonal of the SNF.
          "HNF": The HNF matrix.
          "L": The left transform matrix.
          "labeling": The atomic labeling for the structure.
          "directions": The arrow labeling for the structure.

        minkowskiReduce (bool): Logical indicating if basis should be reduced.

    Returns:
        space_data (dict): A dictionary of the structure mapped to real space with keys:
          "sLV": The super lattice vectors as rows in a matrix.
          "aBas": The attomic basis vectors.
          "spin": A list of the occupanices.
          "gIndx": The integer group index.
          "x": A list of the concentrations of each atom type.
    """
    # print("in system_data",system_data)
    # print("in structure_data",structure_data)
    # print("in minkReduce",minkowskiReduce)
    from numpy import matmul, allclose, matrix, array
    
    n_d = system_data["nD"]
    n = structure_data["n"]

    # DEFINE the non-zero elements of the HNF matrix
    a = structure_data["HNF"][0][0]
    b = structure_data["HNF"][1][0]
    c = structure_data["HNF"][1][1]
    d = structure_data["HNF"][2][0]
    e = structure_data["HNF"][2][1]
    f = structure_data["HNF"][2][2]

    p_bas = system_data["dvecs"]
    S = structure_data["diag"]
    plv = system_data["plattice"]
    HNF = structure_data["HNF"]
    eps = system_data["eps"]
    L = structure_data["L"]
    # Compute the superlattice vectors
    slv = matmul(plv,HNF).tolist()
    # Find the coordinates of the basis atoms
    g_indx = []

    if minkowskiReduce:
        slv = list(map(list,zip(*_minkowski_reduce_basis(list(map(list,zip(*slv))),eps))))
        
    # Find each atomic position from the g-space information
    a_bas = []
    ic = 0  # Keep track of the number of points mapped so far
    # Loop over the number at sites/parent cell (the d set)it
    for iD in range(0, n_d):  
        # For the limits on the loops, see the "interior_points.pdf" write-up
        for z1 in range(a):
            for z2 in range(int((b*z1)/a), int(c+(b*z1)/a)):
                for z3 in range(int(z1*(d-(e*b)/c)/a+(e*z2)/c), int(f+z1*(d-(e*b)/c)/a+(e*z2)/c)):
                    ic +=1
                    # if ic > n: #pragma: no cover
                    #     from .msg import err
                    #     err("Problem with basis atoms in map_enpStr_to_real_space...")
                    #     exit()
                    # Atomic basis vector in Cartesian coordinates
                    temp = matmul(plv,[z1,z2,z3]).tolist()
                    temp2 = [temp[i]+p_bas[iD][i] for i in range(len(p_bas[iD]))]
                    a_bas.append(temp2)
                    # Map position into the group
                    greal = matmul(L,[float(z1),float(z2),float(z3)]).tolist() 
                    g = [int(i) for i in greal] # Convert the g-vector from real to integer
                    if not allclose(greal,g,rtol=eps,atol=eps): #pragma: no cover
                        from .msg import err
                        err("map2G didn't work in map_enumStr_to_real_space")
                        exit()
                    # Bring the g-vector back into the first tile
                    g = [g[i]%S[i] for i in range(len(S))] 
                    # g_indx is the index in the configuration string that
                    # tells us which atom type is used at this position

                    g_indx.append((iD)*S[0]*S[1]*S[2]+g[0]*S[1]*S[2]+g[1]*S[2]+g[2])

    if ic != n*n_d: #pragma: no cover
        from .msg import err
        err("ERROR: map_enumStr_to_real_space: Didn't find the correct # of basis atoms")
        exit()

    k = system_data["k"]
    x = []
    for i in range(k):
        x.append(0.0)
        
    labeling = structure_data["labeling"]
    spin = []
    if k % 2 == 0:
        for iAt in range(0, n*n_d):
            i = int(labeling[g_indx[iAt]])
            digit = i-k//2 # convert 0..k-1 label to spin variable -k/2..k/2
            x[i] += 1  # Keep track of the concentration of each atom type
            if digit < 0:
                spin.append(digit)
            else:
                spin.append(digit+1) # skip 0 as a spin if k is even
    else:
        for iAt in range(0, n*n_d):
            i = int(labeling[g_indx[iAt]])
            spin.append(i-k//2)
            x[i] += 1 # Keep track of the concentration of each atom type

    x = [i/float(n*n_d) for i in x]

    space_data = {"sLV": list(map(list,zip(*slv))), "aBas": a_bas, "spin": spin, "gIndx": g_indx, "x": x}

    return space_data

def _minkowski_reduce_basis(IN,eps):
    """Performs a minkowski reduction on the basis atoms.

    Args:
        IN (list): The input basis as rows of a matrix.
        eps (float): Floating point tolerance.

    Returns:
        OUT (list): The reduced basis vectors as rows of a matrix.

    Raises:
        ValueError: if the input vectors are linearly dependent.
    """

    from numpy import allclose, linalg, array, matrix
    from copy import deepcopy

    limit = 10

    if allclose(linalg.det(IN),0.0,rtol=eps,atol=eps):
        raise ValueError("Input basis for 'minkowski_reduce_basis' was not linearly independent")

    OUT = deepcopy(IN)

    # Keep applying the greedy algorithm until the vectors come out already sorted
    for it in range(1, limit +1):
        # Sort the three vectors into ascending order
        temp = deepcopy(OUT)
        norms = linalg.norm(temp,axis=1).tolist()
        tt = list(range(3))
        tt.reverse()
        for i in tt:
            idx = norms.index(max(norms))
            temp[i] = OUT[idx]
            norms[idx] = 0

        OUT = deepcopy(temp) # Copy the sorted vectors back to OUT
        (OUT[0], OUT[1], OUT[2]) = _reduce_C_in_ABC(OUT[0],OUT[1],OUT[2],eps)
        if linalg.norm(OUT[2]) >= (linalg.norm(OUT[1])-eps):
            break

    if not _minkowski_conditions_check(OUT,eps): #pragma: no cover
        from .msg import err
        err("ERROR in minkowski_reduce_basis: Minkowski conditions not met."
            "Number of iterations: {}".format(str(limit)))
        exit()
    
    # we want to make sure that the det is positive.
    # NOTE: This *destroys* the mathematical picture of a "greedy reduced basis" (Minkowski), but
    #       from a physical point of view we don't care ;-)
    #       Either way, the basis is as orthogonal as possible.
    if linalg.det(OUT) < 0:
        temp[0] = OUT[1]
        OUT[1] = OUT[2]
        OUT[2] = temp[0]

    return OUT

def _minkowski_conditions_check(basis,eps):
    """This function checks the minkowski conditions for a 3D lattice
    basis.

    Args:
        basis (list): The atomic basis vectors.
        eps (float): Finitie precision tolerance.

    Returns:
        minkowski_check (bool): True if minkowski conditions are met.
    """
    from numpy import linalg
    from .msg import err
    
    b1 = basis[0]
    b2 = basis[1]
    b3 = basis[2]

    minkowski_check = True
    if linalg.norm(b1) > (linalg.norm(b2) + eps):
        minkowski_check = False
        err("Minkowski_condition 1 failed: b1 > b2")

    if linalg.norm(b2) > (linalg.norm(b3) + eps):
        minkowski_check = False
        err("Minkowski_condition 2 failed: b2 > b3")

    if linalg.norm(b2) > (linalg.norm([b1[i]+b2[i] for i in range(len(b1))])+eps):
        minkowski_check = False
        err("Minkowski_condition 3 failed: b2 > b1+b2")

    if linalg.norm(b2) > (linalg.norm([b1[i]-b2[i] for i in range(len(b1))])+eps):
        minkowski_check = False
        err("Minkowski_condition 4 failed: b2 > b1-b2")

    if linalg.norm(b3) > (linalg.norm([b1[i]+b3[i] for i in range(len(b1))])+eps):
        minkowski_check = False
        err("Minkowski_condition 5 failed: b3 > b1+b3")

    if linalg.norm(b3) > (linalg.norm([b3[i]-b1[i] for i in range(len(b1))])+eps):
        minkowski_check = False
        err("Minkowski_condition 6 failed: b3 > b3-b1")

    if linalg.norm(b3) > (linalg.norm([b2[i]+b3[i] for i in range(len(b2))])+eps):
        minkowski_check = False
        err("Minkowski_condition 7 failed: b3 > b2+b3")

    if linalg.norm(b3) > (linalg.norm([b3[i]-b2[i] for i in range(len(b2))])+eps):
        minkowski_check = False
        err("Minkowski_condition 8 failed: b3 > b3-b2")

    if linalg.norm(b3) > (linalg.norm([b1[i]+b2[i]+b3[i] for i in range(len(b2))])+eps):
        minkowski_check = False
        err("Minkowski_condition 9 failed: b3 > b1+b2+b3")

    if linalg.norm(b3) > (linalg.norm([b1[i]-b2[i]+b3[i] for i in range(len(b2))])+eps):
        minkowski_check = False
        err("Minkowski_condition 10 failed: b3 > b1-b2+b3")

    if linalg.norm(b3) > (linalg.norm([b1[i]+b2[i]-b3[i] for i in range(len(b2))])+eps):
        minkowski_check = False
        err("Minkowski_condition 11 failed: b3 > b1+b2-b3")

    if linalg.norm(b3) > (linalg.norm([b1[i]-b2[i]-b3[i] for i in range(len(b2))])+eps):
        minkowski_check = False
        err("Minkowski_condition 12 failed: b3 > b1-b2-b3")

    return minkowski_check

def _reduce_C_in_ABC(A,B,C,eps):
    """This routine takes three vectors, A,B,C, defining a lattice, and
    reduces the last one so that it is as close as possible to the
    origin while remaining in an affine plane, which is parallel to
    the A-B plane but which passes through the end of the C
    vector. See Lecture notes in computer science, ISSN 0302-974, ANTS
    - VI : algorithmic number theory, 2004, vol. 3076, pp. 338-357
    ISBN 3-540-22156-5

    Args:
        A (list): The first vector.
        B (list): The second vector.
        C (list): The vector to be reduced.
        eps (float): Finite precision tolerance.

    Returns:
        A (list): The vector A.
        B (list): The vector B.
        C (list): The reduced vector C.
    """
    from numpy import cross, linalg, dot, allclose, matmul, array
    from copy import deepcopy
    from math import floor
    
    old_abc = deepcopy([A,B,C])
    
    # Use Gaussian reduction to reduce the A,B 2D basis so that it is
    # itself Minkowski reduced. If this is done, then the closest
    # lattice point (in A,B plane) to the projection of C (into the
    # A,B plane) is guaranteed to be one of the corners of the unit
    # cell enclosing the projection of C
    (A,B) = _gaussian_reduce_two_vectors(A,B,eps)

    # First thing to do is find the (real, not lattice) point in the
    # affine plane A,B + C that is nearest the origin. Call this T.
    cpd_AB = [i/linalg.norm(cross(A,B)) for i in cross(A,B)]
    T = [C[i] - cpd_AB[i]*dot(C,cpd_AB) for i in range(3)]

    if not allclose(dot(T,cross(A,B)),0,atol=eps,rtol=eps): #pragma: no cover
        from .msg import err
        err("{0} Projection of C into A,B plane failed".format(str(dot(T,cross(A,B)))))

    # Now find the four points of the A,B lattice, in the affine
    # plane, that enclose the point T
    abc = [A,B,C]
    abcinv = linalg.inv(abc)

    LC = [int(floor(i +eps)) for i in matmul(T,abcinv).tolist()]
    # Compute the distance from T to each of the four corners of the cell and pick
    # the one that is the closest.
    corners = array([[0,0,0],[1,0,0],[0,1,0],[1,1,0]])
    dist = []
    for i in range(0,4):
        temp1 = corners[i] + array(LC)
        temp2 = array(T) -matmul((corners[i] + array(LC)),abc)
        dist.append(linalg.norm(array(T) -matmul((corners[i] + array(LC)),abc)))

    idx = dist.index(min(dist))

    if idx == 0:
        temp1 = [corners[0][i] + LC[i] for i in range(3)]
        temp2 = matmul(temp1,abc).tolist()
        C = [C[i] - temp2[i] for i in range(len(C))]
    elif idx == 1:
        temp1 = [corners[1][i] + LC[i] for i in range(3)]
        temp2 = matmul(temp1,abc).tolist()
        C = [C[i] - temp2[i] for i in range(len(C))]
    elif idx == 2:
        temp1 = [corners[2][i] + LC[i] for i in range(3)]
        temp2 = matmul(temp1,abc).tolist()
        C = [C[i] - temp2[i] for i in range(len(C))]
    elif idx == 3:
        temp1 = [corners[3][i] + LC[i] for i in range(3)]
        temp2 = matmul(temp1,abc).tolist()
        C = [C[i] - temp2[i] for i in range(len(C))]
    else: #pragma: no cover
        from .msg import err 
        err("Case failed in reduce_C_in_ABC Lattice coordinates in the A,B plane: {0}".format(' '.join([str(i) for i in LC])))

    abc = [A,B,C]
    abcinv = linalg.inv(abc)
    temp = matmul(list(map(list,zip(*abcinv))),list(map(list,zip(*old_abc)))).tolist()
    for i in range(3):
        for j in range(3):
            if abs(temp[i][j] - int(round(temp[i][j]))) > eps: #pragma: no cover
                from .msg import err
                err("Lattice was not preserved in reduce_C_in_ABC")
                exit()

    return A, B, C

def _gaussian_reduce_two_vectors(U,V,eps):
    """This routine takes two vectors (in three-space) and reduces them to
    form a shortest set (Minkowski reduced). The idea is to subtract
    multiples of U from V so that the new V is as close to the origin
    as any lattice point along the line that passes through U in the
    direction of V. The process is repeated until the new vector isn't
    shorter than the other. It's pretty obvious if you do an example
    by hand. Also see 3.1 of Lecture notes in computer science, ISSN
    0302-974, ANTS - VI: algorithmic number theory, 2004, vol. 3076,
    pp. 338-357 ISBN 3-540-22156-5. Fixes made Apr 2012 GLWH (not sure
    if they made a practical difference though)

    Args:
        U (list): A vector.
        V (list): Another vector.
        eps (float): Finite precision tolerance.

    Returns:
        U (list): The reduced vector U.
        V (list): The reduced vector V.
    """

    from numpy.linalg import norm
    from numpy import dot

    it = 0
    if norm(U) > (norm(V) - eps):
       # Make sure that the {U,V} are listed in ascending order; ||U||<||V||
       temp = U
       U = V
       V = temp # Keep V as the longest vector

    done = False
    it = 1
    while not done:
        if it > 10: # pragma: no cover
            from .msg import err
            err("gaussian_reduce_two_vectors failed to converge in 10 iterations")
            exit()
        R = [V[i]-int(round(dot(U,V)/dot(U,U)+1E-10))*U[i] for i in range(3)] #Shorten V as much as possible
        V = U # Swap U and V (so U remains the shortest)
        U = R
        if norm(U) >= (norm(V) - eps):
            done = True
        it += 1

    # Make sure that the {U,V} are listed in ascending order on exit; ||U||<||V||
    temp = U
    U = V
    V = temp
    return U, V

def cartesian2direct(sLV,aBas, eps):
    """This routine takes three lattice vectors and a list of atomic basis
    vector in Cartesian coordinates and converts them to direct
    ("lattice") coordinates.

    Args:
        sLV (list): The superlattice vectors in cartesian coordinates as rows of a matrix.
        aBas (list): Atomic basis vectors in cartesian coordinates.
        eps (float): Finite precision tolerance.

    Returns:
        aBas (list): The atomic basis vectors in direct coordinates.
    """
    from numpy import linalg, matmul, array
    
    n_at = len(aBas)
    slv_inv = linalg.inv(sLV)
    # Convert aBas to DIRECT COORDINATES
    for iAt in range(n_at):
        aBas[iAt] = matmul(aBas[iAt],slv_inv) # Put positions into
        # "direct" coordinates This keeps the atomic coordinates inside
        # the first unit cell--- not necessary but aesthetically
        # pleasing.

        while any(aBas[iAt] >= (1.0-eps)) or any(aBas[iAt] < (0.0-eps)):
            aBas[iAt] = array([i if i < (1.0-eps) else i-1.0 for i in aBas[iAt]])
            aBas[iAt] = array([i if i >= (0.0-eps) else i+1.0 for i in aBas[iAt]])

    return aBas
