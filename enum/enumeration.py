"""The methods needed to take input from the user and put it in useful
forms. Also the actual executable for the code."""

#import needed modules
from enum.tree import brancher
import enum.phonons as pb
from enum.polyaburnside import polya
import random
import os
import math
import itertools as it
from enum.structures import enum_data
from copy import deepcopy
import enum.io_utils as io

# hard coded error tolerance. This will need to go away and become
# part of the input files later.
eps = 1e-7

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
            label.append(range(volTable[i][0],volTable[i][1]+1))
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
        crange = range(0,size+1)
        aranges = []
        for i in range(nspecies):
            aranges.append(crange)

        p_ranges = it.product(*aranges)
        for p in p_ranges:
            if sum(p) == size:
                # if not any([list(c) in cList for c in it.permutations(p)]) == True:
                cList.append(list(p))
        cList=cList[1:-1]
                    
    return(cList)

# arrow_concs is a method that returns the concentration string
# including the arrows
def arrow_concs(cList,aconcs):
    """Uses the concentrations of the atoms and the arrows to make a
    labeling for the system.

    :arg cListr: an integer array the concentration of the colors
    :arg aconcs: an integer array of the number of arrows for each
    color
    """

    species = 1
    conc_w_arrows = []
    for i in range(len(aconcs)):
        na = aconcs[i]
        ns = cList[i]
        if ns >= na:
            while na > 0:
                conc_w_arrows.append([1,species])
                na -= 1
                ns -= 1
            while ns > 0:
                conc_w_arrows.append([-1,species])
                ns -= 1
        else:
            while ns > 0:
                conc_w_arrows.append([-1,species])
                ns -= 1
        species += 1

    return(conc_w_arrows)

def create_labeling(config):
    """This routine takes a string of atomic locations as a vector and the
    atomic concentrations and returns a unique labeling.

    :arg config: list of integers describing the arrangement of atoms on
    the lattice.
    """
    label = ''
    if isinstance(config[0],list):
        for i in config:
            label += str(i[1]-1)
    else:
        for i in config:
            label += str(i-1)            
    
    return(label)

def _get_arrow_concs(params):
    """If the concentrations are being restricted then find the correct 
    arrow for each species included.

    :arg params: the lattice.in parameters.
    """
    a_concs = []
    if params["is_crestricted"]:
        for c in params["concs"]:
            if (len(c) == 4 and params["arrows"]):
                a_concs.append(c[3])
            else:
                a_concs.append(0)
    else:
        for i in range(params["nspecies"]):
            a_concs.append(0)
    return a_concs

def _polya_out(args):
    """Generates the 'polya.out' files for the cell sizes specified in 'lattice.in'
    (or other specified input file).
    """
    from os import path
    from numpy import zeros
    params = io.read_lattice(args["lattice"])

    for s in range(params["sizes"][0], params["sizes"][1]+1):
        celldir = args["dataformat"].format(s)
        out = open(path.join(celldir, args["outfile"]), 'w+')

        # get HNFs, SNFs, and LTs
        edata = enum_data(s)
        # find the concentrations available for the desired cell sizes.
        cList = get_concs_for_size(s, params["nspecies"], params["is_crestricted"],
                                   len(params["basis_vecs"]), params["concs"])
        
        # We need to write the concs in cList to the output file.
        out.write('{0: <28}'.format("# HNF"))
        for conc in cList:
            out.write("{0: <10}".format(':'.join(map(str, conc))))
        out.write('{0: <10}\n'.format("Total"))

        a_concs = _get_arrow_concs(params)
        conc_totals = zeros(len(cList), int)
        for idata, edict in enumerate(edata):
            HNF = edict["HNF"]
            out.write("  {0: <26}".format(' '.join(map(str, HNF))))
            group = io.read_group(edict["group"])
            # TODO: for arrow enumerations, this is much more complicated, and
            # still needs to be done...
            agroup = [[g,[0]] for g in group]

            # we need to loop over the concentrations and find the
            # number of arrangements possible for each cell size
            total = 0
            for iconc, conc in enumerate(cList):
                if len(conc) > 0:
                    decorations = arrow_concs(conc,a_concs)
                    decorations = pb.col_sort(decorations)
                
                    # we need to know the concentrations of the
                    # species with and without arrows
                    concs_w_arrows = pb.find_concentrations(decorations)
                    # we also need to know the number of arrows and
                    # their species so we can undo the previous step
                    # later
                    (n_arrows,arrow_types) = pb.how_many_arrows(decorations)
                    
                    # now find the number of unique arrangements using Polya.
                    if arrow_types != 0:
                        total_num = polya(concs_w_arrows,agroup,arrowings=arrow_types)
                    else:
                        total_num = polya(conc, agroup)

                    out.write("{0: <10d}".format(total_num))
                    total += total_num
                    conc_totals[iconc] += total_num
            out.write('{0: <10d}\n'.format(total))
        out.write("# " + ''.join(['-' for i in range(len(cList)*10 + 10 + 30)]) + '\n')
        out.write("{0: <28}".format("  0 0 0 0 0 0"))
        for ctotal in conc_totals:
            out.write("{0: <10d}".format(ctotal))
        out.write("{0: <10d}\n".format(sum(conc_totals)))
        out.close()
        
# check which files are present to see if we are finding the HNFs with
# the polya count or enumerating subsets within the HNFs.
def _enum_out(args):
    """Produce the enumerations of a subset of the total number of possible 
    arrangements for the desired HNFs. It assumes that all the information 
    used in the polya part above is still available.
    """
    params = io.read_lattice(args["lattice"])
    systems = io.read_enum(args["input"])
    io.write_enum(params, outfile="enum.out")    

    count_t = 1
    from numpy import unique
    def cellsize(sHNF):
        return sHNF[0]*sHNF[2]*sHNF[5]
    cellsizes = unique([cellsize(sys[0]) for sys in systems])
    datadicts = {}
    sfmt = ("{0: >10d}{1: >10d}{2: >8d}{3: >9d}{4: >9d}{5: >12d}{6: >4d}{7: >6d}"
            "{8: >10}  {9: >18}  {10: >44}    {11}\n")
    def fmtn(l, n):
        return (''.join(["{{{0:d}: >{1:d}d}}".format(i, n) for i in range(len(l))])).format(*l)
    
    with open(args["outfile"], 'a') as f:
        for s in cellsizes:
            dataset = enum_data(s)
            datadicts.update({tuple(d["HNF"]): d for d in dataset})

        for HNF, conc, num_wanted in systems:
            edata = datadicts[tuple(HNF)]
            SNF = edata["SNF"]
            LT = edata["L"]
            
            a_concs = _get_arrow_concs(params)
            configs = enum_sys(edata["group"], list(conc), a_concs, num_wanted)

            for config in configs:
                labeling = create_labeling(config)
                o = sfmt.format(count_t, 1, 1, 1, 1, 1, sum(conc), 1,
                                fmtn(SNF, 3), fmtn(HNF, 3),
                                fmtn(LT, 5), labeling)
                f.write(o)
                count_t += 1

def enum_sys(groupfile, concs, a_concs, num_wanted):
    """Enumerates a random subset of the unique structures that have the shape
    defined by the symmetry group and the specified concentration.

    :arg groupfile: path to the file containing the symmetry group.
    :arg concs: list of integer concentrations for each species.
    :arg a_concs: list of integer *arrow* concentrations for each species.
    :arg num_wanted: the number of structures to pick randomly from the enumerated
      list.
    """
    decorations = arrow_concs(concs, a_concs)
    decorations = pb.col_sort(decorations)
    # get the symmetry group for this HNF. Assumes the group can be
    # found in the file labeled by (this_HNF)_sym_group.out
    group = io.read_group(groupfile)
    # get symgroup from HNF and lat_vecs
    # add [0] to each element of the symmetry group
    agroup = [[g,[0]] for g in group]

    # we need to know the concentrations of the
    # species with and without arrows
    concs_w_arrows = pb.find_concentrations(decorations)
    # we also need to know the number of arrows and
    # their species so we can undo the previous step
    # later
    (n_arrows, arrow_types) = pb.how_many_arrows(decorations)

    # now find the number of unique arrangements using
    # polya
    if arrow_types != 0:
        total = polya(concs_w_arrows, agroup, arrowings=arrow_types)
    else:
        total = polya(concs, agroup)
        
    # generate the random subset to be used
    if num_wanted < total:
        from random import shuffle
        subset = range(1, total+1) 
        shuffle(subset)
        subset = subset[0:num_wanted]
    else:
        from msg import warn
        warn("number of configurations requested exceeds the number of "
            "unique configurations available.")
        subset = []

    # here we get the configs and the len of the stabilizers if we're
    # doing not doing a purely arrowed enumeration. Getting the
    # stabilizers allows us to remove the superperoidic structures.
    n_stabs = []
    if len(concs) == 1 and all(decorations) >=0:
        configs = []
        a_configs = pb.add_arrows(decorations, agroup, 6)
        count = 1
        for config in a_configs:
            if count in subset:
                configs.append(config)
            count += 1
    else:
        (configs,n_stabs) = brancher(concs, agroup, decorations, 6, subset)

    # reduced_configs is the list of configurations with the
    # superperiodic configurations removed.
    reduced_configs = []
    # need to find a way to remove the superperiodic arrangements
    reduced_configs = configs                

    return reduced_configs

def _examples():
    """Print some examples on how to use this python version of the code."""
    helptext = ("For all the examples below, it is assumed that you have already "
                "compiled and ran the modified enumlib code as described in the "
                "README or in some other manner obtained the HNFs (supercells) and "
                "their symmetry groups. You will then need to specify if you are "
                "obtaining the number of unique arrangements for each supercell and "
                "concentration allowed for your system or enumerating (finding) the "
                "desired number of configurations for each HNF and concentration. "
                "Additionally you way change the default input file names to ones of "
                "your own creation.")
    egs = [("Find the Polya distribution",
            "The code below finds the number of unique arrangements for each "
            "supercell (HNF) for a binary system on an fcc lattice which can have "
            "1 to 11 atoms in the supercell as described in the lattice.in found in "
            "input/fcc sample directory. The files labeled cell.n contain the HNFs "
            "and symmetry group for the supercells of size n and were generated by "
            "the modified enum.x code. For more information on the contents of this "
            "folder please see the README. To a different input file to use rather "
            "than lattice.in use the -lattice option or if you have the HNF and "
            "symmetry group data in a different file then cells.n then use "
            "-dataformat.","./enumeration.py -polya"),
           ("Find the desired number of unique structures",
            "This code also uses the sample system and files found in "
            "input/fcc. It finds the desired number of unique arrangements for each "
            "HNF and concentration range listed in the enum.in file. This file has "
            "the format:\n # HNF                           Conc.       Number\n"
            "   1 0 1 0 0 4                   1 3         1\n "
            "Where HNF is the lower triangle of the HNF Conc is the concentrations "
            "of the 2 atoms on the lattice and Number specifies the number of "
            "arrangements for that HNF and concentration range the user desires."
            , "./enumeration.py -enum")]

    print("POLYA ENUMERATION THEOREM SOLVER\n")
    for eg in egs:
        title, desc, code = eg
        print("--" + title + '--\n')
        print(desc + '\n')
        print('  ' + code + '\n')

def _parser_options(phelp=False):
    """Parses the options and arguments from the command line."""
    import argparse
    parser = argparse.ArgumentParser(description="Partial Superstructure Enumeration Code")
    parser.add_argument("-debug", action="store_true",
                        help="Print verbose calculation information for debugging.")
    parser.add_argument("-examples", action="store_true",
                        help="Print some examples for how to use the enumeration code.")
    parser.add_argument("-polya", action="store_true",
                        help=("Predict the total number of unique superstructures for each "
                              "cell shape and possible concentration."))
    parser.add_argument("-enum", action="store_true",
                        help=("Enumerate the number of superstructures of specific shape "
                              "as specified in the 'enum.in' file."))
    parser.add_argument("-lattice",
                        help=("Override the default input file name: 'lattice.in' for "
                              "enumeration parameters."))
    parser.add_argument("-input",
                        help=("Override the default 'enum.in' file name."))
    parser.add_argument("-dataformat", default="cells.{}",
                        help=("Specify the default folder name for any cell size that contains "
                              "the matrices and groups generated by 'enum.x'. Format is: cells.{} "
                              "where {} is a placeholder for the integer cell size."))
    parser.add_argument("-verbose", type=int,
                        help="Specify the verbosity level (1-3) for additional computation info.")
    parser.add_argument("-outfile",
                        help=("Override the default output file names: 'polya.out' for polya counting; "
                              "'enum.out' for structure enumeration."))
    vardict = vars(parser.parse_args())
    if phelp or vardict["examples"]:
        _examples()
        exit(0)

    if vardict["verbose"]:
        from enum.msg import set_verbosity
        set_verbosity(vardict["verbose"])

    if not vardict["lattice"]:
        vardict["lattice"] = "lattice.in"
    if not vardict["input"]:
        vardict["input"] = "enum.in"
    if not vardict["outfile"]:
        vardict["outfile"] = "polya.out" if vardict["polya"] else "enum.out"
    return vardict

def script_enum(args):
    """Generates the 'polya.out' or 'enum.out' files depending on the script arguments.
    """
    from os import path
    if args["polya"]:
        #Perform validation for running polya.
        if not path.isfile(args["input"]):
            from enum.msg import err
            err("The input file {} does not exist.".format(args["lattice"]))
            exit()
            
    if args["enum"]:
        #Perform validation for running enum.
        if not path.isfile(args["input"]):
            from enum.msg import err
            err("The input file {} does not exist.".format(args["input"]))
            exit()
            
    if args["enum"] or args["polya"]:
        from glob import glob
        #Perform validation for running enum.
        if len(glob(args["dataformat"].split('.')[0]+'.*')) < 1:
            from enum.msg import err
            err("The input folders {} do not exist.".format(args["dataformat"]))
            exit()
            
    if args["polya"]:
        _polya_out(args)
    if args["enum"]:
        _enum_out(args)
    
if __name__ == '__main__':
    script_enum(_parser_options())
