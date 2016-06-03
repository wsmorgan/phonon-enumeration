"""The methods needed to take input from the user and put it in useful
forms. Also the actual executable for the code."""
# hard coded error tolerance. This will need to go away and become
# part of the input files later.

def _enum_in(args):
    """Makes an enum.in file that contains the desired distribution of
    structures from the polya distribution.

    :arg args: The command line inputs.
    """

    from phenum.structures import make_enum_in

    distribution = args["distribution"][0].lower()
    if args["distribution"][1].lower() == "all":
        make_enum_in(distribution,dataformat=args["dataformat"])
    else:
        make_enum_in(distribution,number=int(args["distribution"][1]),
                     dataformat=args["dataformat"])

def _polya_out(args):
    """Generates the 'polya.out' files for the cell sizes specified in 'lattice.in'
    (or other specified input file).
    """
    from os import path, makedirs
    from numpy import zeros
    from phenum.grouptheory import get_sym_group, get_full_HNF
    from phenum.msg import warn
    from phenum.symmetry import get_concs_for_size
    from phenum.io_utils import read_lattice, read_group
    from phenum.structures import enum_data
    import phenum.phonons as pb
    from phenum.polyaburnside import polya
    params = read_lattice(args["lattice"])

    for s in range(params["sizes"][0], params["sizes"][1]+1):
        celldir = args["dataformat"].format(s)

        # get HNFs, SNFs, and LTs
        edata = enum_data(s,args,params)
        
        out = open(path.join(celldir, args["outfile"]), 'w+')
            

        # find the concentrations available for the desired cell sizes.
        cList = get_concs_for_size(s, params["nspecies"], params["is_crestricted"],
                                   len(params["basis_vecs"]), params["concs"])
        
        # We need to write the concs in cList to the output file.
        out.write('{0: <28}'.format("# HNF"))
        for conc in cList:
            out.write("{0: <10}".format(':'.join(map(str, conc))))
        out.write('{0: <10}\n'.format("Total"))

        a_concs = pb.get_arrow_concs(params)
        conc_totals = zeros(len(cList), int)
        for idata, edict in enumerate(edata):
            HNF = edict["HNF"]
            out.write("  {0: <26}".format(' '.join(map(str, HNF))))

            if params["arrows"]:
                sym_g = get_sym_group(params["lat_vecs"],params["basis_vecs"],get_full_HNF(HNF),3)
                agroup = []
                for i in range(len(sym_g.perm.site_perm)):
                    agroup.append([sym_g.perm.site_perm[i],sym_g.perm.arrow_perm[i]])
            else:
                group = read_group(edict["group"])
                agroup = [[g,[0]] for g in group]

            # we need to loop over the concentrations and find the
            # number of arrangements possible for each cell size
            total = 0
            for iconc, conc in enumerate(cList):
                if len(conc) > 0:

                    decorations = pb.arrow_concs(conc,a_concs)
                
                    # we need to know the concentrations of the
                    # species with and without arrows, we also need to
                    # know the number of arrows and their species so
                    # we can undo the previous step later
                    (n_arrows,arrow_types,concs_w_arrows) = pb.how_many_arrows(decorations)

                    # now find the number of unique arrangements using Polya.
                    if arrow_types != 0:
                        # Since enumlib doesn't write the arrow group
                        # out we have to recompute the group actions
                        # paired with their effects on the arrows
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

    import phenum.io_utils as io
    from phenum.structures import enum_data
    import phenum.phonons as pb
    from numpy import unique
    from phenum.grouptheory import get_full_HNF, SmithNormalForm
    
    params = io.read_lattice(args["lattice"])
    systems = io.read_enum(args["input"])
    io.write_enum(params, outfile="enum.out")    

    count_t = 1
    def cellsize(sHNF):
        return sHNF[0]*sHNF[2]*sHNF[5]
    cellsizes = unique([cellsize(sys[0]) for sys in systems])
    datadicts = {}
    sfmt = ("{0: >10d}{1: >10d}{2: >8d}{3: >9d}{4: >9d}{5: >12d}{6: >4d}{7: >6d}"
            "{8: >10}  {9: >18}  {10: >44}    {11}    {12: >21}\n")
    def fmtn(l, n):
        return (''.join(["{{{0:d}: >{1:d}d}}".format(i, n) for i in range(len(l))])).format(*l)
    
    with open(args["outfile"], 'a') as f:
        if not params["arrows"]:
            for s in cellsizes:
                dataset = enum_data(s,args,params)
                datadicts.update({tuple(d["HNF"]): d for d in dataset})

        for HNF, conc, num_wanted in systems:
            if not params["arrows"]:
                edata = datadicts[tuple(HNF)]
                SNF = edata["SNF"]
                LT = edata["L"]
            
                a_concs = pb.get_arrow_concs(params)
                configs = pb.enum_sys(edata["group"], list(conc), a_concs, num_wanted,HNF,params)
                
            else:
                (SNF,L,R) = SmithNormalForm(get_full_HNF(HNF))
                SNF = [SNF[0][0],SNF[1][1],SNF[2][2]]
                LT = [item for row in L for item in row]
                a_concs = pb.get_arrow_concs(params)
                configs = pb.enum_sys(None, list(conc), a_concs, num_wanted,HNF,params)

            for config in configs:
                (labeling,arrowing) = io.create_labeling(config)
                o = sfmt.format(count_t, 1, 1, 1, 1, 1, sum(conc), 1,
                                fmtn(SNF, 3), fmtn(HNF, 3),
                                fmtn(LT, 5), labeling, arrowing)
                f.write(o)
                count_t += 1

def _examples():
    """Print some examples on how to use this python version of the code."""
    helptext = ("For all the examples below, it is assumed that you have already "
                "compiled the modified enumlib code as described in the "
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
            "-dataformat.","enumeration.py -polya"),
           ("Construct an enum.in file before running the -enum mode",
            "This code assumes that the -polya mode has arleady been run. It takes two"
            " arguments; the first is the disered distribution type ('size', 'conc',"
            " 'shape', 'all'), the second is the desired number of structures, if all"
            " the structures are wanted then the second argument should be 'all'.",
            "enumeration.py -distribution all all \n  enumeration.py -distribution "
            " size 200"),
           ("Find the desired number of unique structures",
            "This code also uses the sample system and files found in "
            "input/fcc. It finds the desired number of unique arrangements for each "
            "HNF and concentration range listed in the enum.in file. This file has "
            "the format:\n # HNF                           Conc.       Number\n"
            "   1 0 1 0 0 4                   1 3         1\n "
            "Where HNF is the lower triangle of the HNF Conc is the concentrations "
            "of the 2 atoms on the lattice and Number specifies the number of "
            "arrangements for that HNF and concentration range the user desires."
            , "enumeration.py -enum")]

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
    parser.add_argument("-exec",
                        help=("Override the default 'enum.x' executable with an executable "
                              "name or path."))
    parser.add_argument("-distribution", nargs= "+",
                        help=("Makes an enum.in file when the distribution is 'all'. Otherwise "
                              "the distribution is printed to the screen for the user. "
                              " The distributions are built from the results of the polya "
                              "run and contain the specified number of structures."))
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
        from phenum.msg import set_verbosity
        set_verbosity(vardict["verbose"])

    if not vardict["lattice"]:
        vardict["lattice"] = "lattice.in"
    if not vardict["input"]:
        vardict["input"] = "enum.in"
    if not vardict["exec"]:
        vardict["exec"] = "enum.x"
    if not vardict["outfile"]:
        vardict["outfile"] = "polya.out" if vardict["polya"] else "enum.out"
    return vardict

def _script_enum(args):
    """Generates the 'polya.out' or 'enum.out' files depending on the script arguments.
    """
    from os import path, system
    if args["polya"]:
        #Perform validation for running polya.
        if not path.isfile(args["lattice"]):
            from phenum.msg import err
            from os import listdir
            err("The input file {} does not exist.".format(args["lattice"]))
            exit()
            
    if args["enum"]:
        #Perform validation for running enum.
        if not path.isfile(args["input"]):
            from phenum.msg import err
            err("The input file {} does not exist.".format(args["input"]))
            exit()

    if args["enum"] or args["polya"]:
        from glob import glob
        #Perform validation for running enum.
        if len(glob(args["dataformat"].split('.')[0]+'.*')) < 1:
            from phenum.msg import err, warn
            from phenum.io_utils import which, read_lattice, write_struct_enum
            from os import system
            warn("The input folders {} do not exist.".format(args["dataformat"]))
            warn("Now running your enumlib executable to build the folders.")
            if which(args["exec"]) != None:
                write_struct_enum(read_lattice(args["lattice"]))
                system(args["exec"])
                system('rm symops_enum_parent_lattice.out readcheck_enum.out fort.*')
                if len(glob(args["dataformat"].split('.')[0]+'.*')) < 1:
                    err("The executable you have for {} does not produce "
                        "the needed input folders {}. In order to correct this "
                        "you need to follow the compilation instructions found "
                        "in the README.".format(args["exec"],args["dataformat"]))
                    exit()                            
            else:
                err("Could not find {} on your path. Please add it to your path "
                    "if you have already\n compiled it. Otherwise please follow the "
                    "instructions found in the README to download, make,\n and place the "
                    "executable in your path.".format(args["exec"]))
                exit()

    if args["distribution"]:
        if len(args["distribution"]) != 2:
            from phenum.msg import err
            err("The distribution option takes two arguments. The parameters the distribution "
                " is over ('shape', 'conc', 'size', 'all') and the number of structures desired."
                "If all the structures are wanted then the second argument should be 'all'.")
            exit()
            
    if args["polya"]:
        _polya_out(args)
    if args["enum"]:
        _enum_out(args)
    if args["distribution"]:
        _enum_in(args)
        
if __name__ == '__main__':
    _script_enum(_parser_options())
