"""Methods for reading and writing enumeration and polya-counting results."""
def read_group(fname):
    """Reads the symmetry group in from the 'rot_perms' styled group
    output by enum.x.

    Args:
        fname (str): Path to the file to read the group from.

    Returns:
        group (numpy ndarray): The symmetry group.
    """
    i=0
    groupi = []
    with open(fname) as f:
        for line in f:
            if i > 5:
                if ('Perm #:') in line:
                    groupi.append(list(map(int, line.split()[4::])))
                else:
                    groupi[-1] += list(map(int, line.split()))
            i += 1
    from numpy import array
    return(list(map(list, array(groupi)-1)))

def read_lattice(filename='lattice.in'):
    """Reads the lattice.in file.

    Args:
        filename (str, optional): The filename to be read in. Default is 'lattice.in'.

    Returns: 
        result (dict): A dictionary with the following fields:
          "sizes": the range of cell sizes,
          "lat_vecs": lattice vectors of the parent cell,
          "nspecies": the number of atomic species in the enumeration,
          "basis_vecs": basis vectors of the parent cell,
          "is_crestricted": logical that indicates if the concentrations will be restricted,
          "arrows": logical that indicates if arrows are present,
          "concs": array of the concentrations in format [1,3].
    """
    from .msg import info
    with open(filename,'r') as f:
        lines = [l for l in f.readlines() if l.strip()[0] != '#']
        
    sizes = [int(el) for el in lines[0].split(' ')]
    bulk = lines[1].strip()[0].lower() == 'b'
    lat_vecs = [list(map(float, l.split())) for l in lines[2:5]]
    nspecies = int(lines[5])
    nbas = int(lines[6])
    bas_vecs = [list(map(float, l.split())) for l in lines[7:7+nbas]]

    arrows = False
    res_conc = False
    if lines[7+nbas].strip() == 'F':
        info('No concentration restrictions given.')
        concs = []
    else:
        res_conc = True
        info('Concentrations are being restricted.')
        arrows = lines[7+nbas+1].strip()
        if arrows == 'F':
            arrows = []
            info('No displacement directions are included')
        else:
            arrows = True
        concs = [list(map(float, l.split())) for l in lines[7+nbas+2:]]

    result = {
        "bulk": bulk,
        "sizes": sizes,
        "lat_vecs": lat_vecs,
        "nspecies": nspecies,
        "basis_vecs": bas_vecs,
        "is_crestricted": res_conc,
        "arrows": arrows,
        "concs": concs
    }
    return result

# This method reads in the enum.in file.
def read_enum(filename="enum.in"):
    """Reads the list of structures and the number of random candidates to enumerate
    from the 'enum.in' styled file. It should contain only integers and have:

       HNF concs num_wanted

    where,

    HNF: the list of 6 independent entries in the HNF matrix;
    concs: list of integer species concentrations;
    num_wanted: number of random structures to draw from the enumerated list.

    Args:
        filename (str, optional): The filename and path. Default 'enum.in'.
    
    Returns:
        systems (list of lists): A list containing the lists of [HNF, concs, num_wanted].
    """    
    from numpy import loadtxt
    raw = loadtxt(filename, int, ndmin=2)
    systems = []
    for i in raw:
        HNF = i[0:6]
        sys_conc = i[6:-1]
        num_wanted = i[-1]
        systems.append((HNF, sys_conc, num_wanted))
    return systems

def write_enum(params, outfile="enum.out"):
    """Writes a 'struct_enum.out' style preamble to the specified output file
    using the enumeration parameters gleaned from a 'lattice.in' styled file.

    Args:
        params (dict): values returned from method:read_lattice().
        outfile (str, optional): path to desired output file. Default is 'enum.out'.
    """
    from datetime import datetime
    lines = []
    lines.append("Random structure enumeration: {0}".format(str(datetime.now())))
    lines.append("bulk" if params["bulk"] else "surface")
    lines.extend([" {0:.7f}       {1:.7f}       {2:.7f}".format(*l) for l in params["lat_vecs"]])
    lines.append("    {0:d}".format(len(params["basis_vecs"])))
    lines.extend(["  {0:.7f}       {1:.7f}       {2:.7f}".format(*l) for l in params["basis_vecs"]])
    lines.append(" {0:d}-nary case".format(params["nspecies"]))
    lines.append("   " + '   '.join(map(str, params["sizes"])))
    lines.append(" 1e-7")
    lines.append("Concentration check:")
    lines.append("    T" if params["is_crestricted"] else "    F")
    lines.extend(["full list of labelings (including incomplete labelings) is used ",
                  "Equivalency list: 1 #Not used in the random enumeration algorithm ",
                  "start   #tot      HNF     Hdegn   labdegn   Totdegn   #size idx    "
                  "pg    SNF             HNF                 Left "
                  "transform                          labeling                 "
                  "directions",
                  ""])
        
    with open(outfile, 'w') as f:
        f.write('\n'.join(lines))

def write_struct_enum(params):
    """Writes a 'struct_enum.in' file for the executable enum.x.

    Args:
        params (dict): values returned from method:read_lattice().
    """

    with open("struct_enum.in","w+") as struct_file:
        struct_file.write("System \n")

        if params["bulk"] == True:
            struct_file.write("bulk \n")
        else:
            struct_file.write("surface \n")

        latt = ""
        i = 0
        for vec in params["lat_vecs"]:
            for point in vec:
                latt += str(point) + " "
            if i < 2:
                latt += "\n"
            i += 1

        struct_file.write(latt + " \n")

        struct_file.write(str(params["nspecies"]) + " \n")

        occupancy = ""
        for i in range(params["nspecies"]):
            if i < params["nspecies"]-1:
                occupancy += str(i) +"/"
            else:
                occupancy += str(i)
    
        basis = ""
        for vec in params["basis_vecs"]:
            basis += "{}  {} \n".format(" ".join([str(point) for point in vec]),occupancy)
        
        struct_file.write(str(len(params["basis_vecs"])) + " \n")
        struct_file.write(basis + " \n")

        size_range = ""
        i = 0
        for point in params["sizes"]:
            size_range += str(point) + " "

        struct_file.write(size_range + " \n")
        struct_file.write("0.00001 \n")
        struct_file.write("full \n")

def which(program):
    """Determines if and where an executable exists on the users path.
    This code was contributed by Jay at http://stackoverflow.com/a/377028

    Args:
        program (str): The name, or path for the program.

    Returns:
        The program or executable.
    """
    import os
    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(program)
    if fpath and is_exe(program):
        return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            path = path.strip('"')
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file

    return None

def create_labeling(config):
    """This routine takes a string of atomic locations as a vector and the
    atomic concentrations and returns a unique labeling.

    Args:
        config (list of int): A list of integers describing the arrangement of atoms on
          the lattice.

    Returns:
        label (str): The atomic labeling for the arrangement.
        arrow (str): The arrow directions for the arrangement.
    """
    label = ''
    arrow =  ''
    if isinstance(config[0],list):
        for i in config:
            label += str(i[1]-1)
            arrow += str(i[0]+1)
    else:
        for i in config:
            label += str(i-1)
            arrow += '0'
    
    return(label,arrow)

def read_enum_out(args):
    """Reads the enum.out file and builds a dictionary with the needed
    information to construct a POSCAR.
    
    Args:
        args (dict): The makeStr.py input arguments.

    Returns:
        system (dict): A dictionary with the keys:
          "plattice": The parrent lattice vectors.
          "dvecs": The atomic basis vectors.
          "title": The name of the system.
          "bulksurf": 'bulk' or 'surface'.
          "nD": The number of atomic basis vectors.
          "k": The number of atomic species.
          "eps": Finite precision tolerance.

        structure_data (list of dict): A list of dictionaries for each desired structure. Each
          dictionary contains:
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
     """
    from numpy import transpose
    
    # which structures are wanted
    if args["structures"] == None:
        with open(args["input"],"r") as f:
            for count, l in enumerate(f):
                pass
        structures = list(range(1,count-13))
    else:
        structures = args["structures"]
    # open the enum.out style file.
    structf = open(args["input"],"r")

    # we'll build a dictionary of the system data and a list of
    # dictionaries for the structures that are wanted.
    structure_data = []
    system = {}
    system["plattice"] = []
    system["dvecs"] = []
    line_count = 1
    system["nD"] = 0
    for line in structf:
        temp = line.rstrip()
        if not temp.startswith("#"):
            if line_count == 1:
                system["title"] = temp
            if line_count == 2:
                system["bulksurf"] = temp
            if line_count in [3,4,5]:
                vec = [float(v) for v in temp.split()]
                system["plattice"].append(vec)
            if line_count == 6:
                system["nD"] = int(temp.rstrip())
            if system["nD"] != 0 and line_count in range(7,7+system["nD"]):
                vec = [float(v) for v in temp.split()]
                system["dvecs"].append(vec)
            if line_count == 7+system["nD"]:
                system["k"] = int(temp.split('-')[0].strip())
            if line_count == 9 + system["nD"]:
                system["eps"] = float(temp.strip())
            if line_count - (14 + system["nD"]) in structures:
                data = temp.split()
                this_struct = {}
                this_struct["strN"] = int(data[0])
                this_struct["hnfN"] = int(data[1])
                this_struct["hnf_degen"] = int(data[2])
                this_struct["lab_degen"] = int(data[3])
                this_struct["tot_degen"] = int(data[4])
                this_struct["sizeN"] = int(data[5])
                this_struct["n"] = int(data[6])
                this_struct["pgOps"] = int(data[7])
                this_struct["diag"] = [int(data[8]),int(data[9]),int(data[10])]
                this_struct["HNF"] = [[int(data[11]),0,0],[int(data[12]),int(data[13]),0],
                                      [int(data[14]),int(data[15]),int(data[16])]]
                this_struct["L"] = [[int(data[17]),int(data[18]),int(data[19])],
                                    [int(data[20]),int(data[21]),int(data[22])],
                                      [int(data[23]),int(data[24]),int(data[25])]]
                this_struct["labeling"] = data[26]
                this_struct["directions"] = data[27]                
                structure_data.append(this_struct)
        line_count += 1

    system["plattice"] = transpose(system["plattice"])
        
    return (system, structure_data)

def write_POSCAR(system_data,space_data,structure_data,args):
    """Writes a vasp POSCAR style file for the input structure and system
    data.

    Args:
        system_data (dict): A dictionary of the system_data.
        space_data (dict): A dictionary containing the spacial data
        structure_data (dict):: a dictionary of the data for this structure
        args (dict): Dictionary of user supplied input.
    """

    from numpy import array
    from phenum.element_data import get_lattice_parameter
    from random import uniform
    from itertools import product

    if "{}" in args["outfile"]:
        filename = args["outfile"].format(str(structure_data["strN"]))
    else:
        filename = args["outfile"] + ".{}".format(str(structure_data["strN"]))

    labeling = structure_data["labeling"]            
    g_indx = space_data["gIndx"]
    arrows = structure_data["directions"]
    struct_n = structure_data["strN"]

    arrow_directions = [[0,0,0],[0,0,-1],[0,-1,0],[-1,0,0],[1,0,0],[0,1,0],[0,0,1]]
    directions = []

    concs = [0]*system_data["k"]
    for atom in range(structure_data["n"]*system_data["nD"]):
        concs[int(labeling[g_indx[atom]])] += 1
        
    def_title = "{} str #: {}\n".format(str(system_data["title"]),str(structure_data["strN"]))

    lattice_parameter, title = get_lattice_parameter(args["species"],concs,def_title)

    for arrow in arrows:
        directions.append(array(arrow_directions[int(arrow)]))
    slv = space_data["sLV"]
    with open(filename,"w+") as poscar:
        poscar.write(title)
        poscar.write("{0:.2f}\n".format(lattice_parameter))
        for i in range(3):
            poscar.write(" {}\n".format(" ".join(
                ["{0: .8f}".format(j) for j in slv[i]])))
        poscar.write("  ")
        if args["species"] == None:
            for ic in concs:
                poscar.write("{}   ".format(str(ic)))
        else:
            for ic in concs:
                if ic != 0:
                    poscar.write("{}   ".format(str(ic)))                    

        poscar.write("\n")
        poscar.write("D\n")
        for ilab, iAt in product(range(system_data["k"]),range(structure_data["n"]*system_data["nD"])):
            rattle = uniform(-args["rattle"],args["rattle"])
            displace = directions[iAt]*args["displace"]*lattice_parameter
            displace += displace*rattle
            if labeling[g_indx[iAt]] == str(ilab):
                out_array = array(space_data["aBas"][iAt]) + displace
                poscar.write(" {}\n".format(
                    "  ".join(["{0: .8f}".format(i) for i in out_array.tolist()])))
