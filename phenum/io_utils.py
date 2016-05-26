"""Methods for reading and writing enumeration and polya-counting results."""
def read_group(fname):
    """Reads the symmetry group in from the 'rot_perms' styled group
    output by enum.x.

    :arg fname: path to the file to read the group from.
    """
    i=0
    groupi = []
    with open(fname) as f:
        for line in f:
            if i > 5:
                if ('Perm #:') in line:
                    groupi.append(map(int, line.split()[4::]))
                else:
                    groupi[-1] += map(int, line.split())
            i += 1
    from numpy import array
    return(map(list, array(groupi)-1))

def read_lattice(filename='lattice.in', verbose=False):
    """Reads the lattice.in file; returns a dictionary with the following fields:
      "sizes": the range of cell sizes,
      "lat_vecs": lattice vectors of the parent cell,
      "nspecies": the number of atomic species in the enumeration,
      "basis_vecs": basis vectors of the parent cell,
      "is_crestricted": logical that indicates if the concentrations will be restricted,
      "arrows": logical that indicates if arrows are present,
      "concs": array of the concentrations in format [1,3].
    """
    from msg import info
    with open(filename,'r') as f:
        lines = [l for l in f.readlines() if l.strip()[0] != '#']
        
    sizes = [int(el) for el in lines[0].split(' ')]
    bulk = lines[1].strip()[0].lower() == 'b'
    lat_vecs = [map(float, l.split()) for l in lines[2:5]]
    nspecies = int(lines[5])
    nbas = int(lines[6])
    bas_vecs = [map(float, l.split()) for l in lines[7:7+nbas]]

    arrows = False
    res_conc = False
    if lines[7+nbas].strip() == 'F':
        info('No concentration restrictions given.')
        concs = []
    else:
        res_conc = True
        info('Concentrations are being restrictied.')
        arrows = lines[7+nbas+1].strip()
        if arrows == 'F':
            arrows = []
            info('No displacement directions are included')
        else:
            arrows = True
        concs = [map(float, l.split()) for l in lines[7+nbas+2:]]

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
    """    
    from numpy import loadtxt
    raw = loadtxt(filename, int, ndmin=2)
    systems = []
    for i in range(len(raw)):
        HNF = raw[i,0:6]
        sys_conc = raw[i,6:-1]
        num_wanted = raw[i,-1]
        systems.append((HNF, sys_conc, num_wanted))
    return systems

def write_enum(params, outfile="enum.out"):
    """Writes a 'struct_enum.out' style preamble to the specified output file
    using the enumeration parameters gleaned from a 'lattice.in' styled file.

    :arg params: values returned from method:read_lattice().
    :arg outfile: path to desired output file.
    """
    from datetime import datetime
    lines = []
    lines.append("Random structure enumeration: %s" % (str(datetime.now())))
    lines.append("bulk" if params["bulk"] else "surface")
    lines.extend(["{0:.7f}       {1:.7f}       {2:.7f}".format(*l) for l in params["lat_vecs"]])
    lines.append("{0:d}".format(len(params["basis_vecs"])))
    lines.extend(["{0:.7f}       {1:.7f}       {2:.7f}".format(*l) for l in params["basis_vecs"]])
    lines.append("{0:d}-nary case".format(params["nspecies"]))
    lines.append('   '.join(map(str, params["sizes"])))
    lines.append("1e-7")
    lines.append("Concentration check:")
    lines.append("T" if params["is_crestricted"] else "F")
    lines.extend(["full list of labelings (including incomplete labelings) is used ",
                  "Equivalency list: 1 #Not used in the random enumeration algorithm ",
                  "start   #tot      HNF     Hdegn   labdegn   Totdegn   #size idx    "
                  "pg    SNF             HNF                 Left "
                  "transform                          labeling                 "
                  "directions",
                  ""])
        
    with open(outfile, 'w') as f:
        f.write('\n'.join(lines))
