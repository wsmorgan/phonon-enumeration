#!/Users/trunks/envs/py3/bin/python3.4
from phenum import msg

"""The methods needed to take input from the user and put it in useful
forms. Also the actual executable for the code."""
# hard coded error tolerance. This will need to go away and become
# part of the input files later.

def _plot_HNFs(args,testmode=False):
    """Makes plots of the shapes of the HNF's provided in the enum.in input file.

    :arg args: The command line inputs.
    """
    if args["shapes"]:
        from phenum.visualize import HNF_shapes
        HNF_shapes(args["input"],args["lattice"],args["show"],testmode=testmode)
    else:
        from phenum.visualize import HNF_atoms
        HNF_atoms(args["input"],args["lattice"],args["show"],testmode=testmode)    

def _enum_in(args):
    """Makes an enum.in file that contains the desired distribution of
    structures from the polya distribution.

    :arg args: The command line inputs.
    """

    from phenum.structures import make_enum_in

    distribution = args["distribution"][0].lower()
    if args["sizes"]:
        sizes = list(range(*args["sizes"]))
    else:
        sizes = None
    if args["distribution"][1].lower() == "all":
        make_enum_in(distribution,args["cellsdir"],dataformat=args["dataformat"],sizes=sizes,
                     outfile=args["outfile"],save= True if args["savedist"] else False,seed=args["seed"],
                     restrict=args["filter"])
    else:
        make_enum_in(distribution,args["cellsdir"],number=int(args["distribution"][1]),
                     dataformat=args["dataformat"],sizes=sizes,outfile=args["outfile"],
                     save= True if args["savedist"] else False,seed=args["seed"],restrict=args["filter"])

def _polya_out(args):
    """Generates the 'polya.out' files for the cell sizes specified in 'lattice.in'
    (or other specified input file).
    """
    from os import path, makedirs
    from phenum.grouptheory import get_sym_group, get_full_HNF
    from phenum.msg import warn
    from phenum.symmetry import get_concs_for_size
    from phenum.io_utils import read_lattice, read_group
    from phenum.structures import enum_data
    import phenum.phonons as pb
    from phenum.polyaburnside import polya
    params = read_lattice(args["lattice"])

    for s in range(params["sizes"][0], params["sizes"][1]+1):
        celldir = "{}{}".format(args["cellsdir"],args["dataformat"].format(s))
        # get HNFs, SNFs, and LTs
        edata = enum_data(s,args,params)
        
        out = open(path.join(celldir, args["outfile"]), 'w+')
            
        # find the concentrations available for the desired cell sizes.
        cList = get_concs_for_size(s, params["nspecies"], params["is_crestricted"],
                                   len(params["basis_vecs"]), params["concs"])
        
        # We need to write the concs in cList to the output file.
        out.write('{0: <28}'.format("# HNF"))
        for conc in cList:
            out.write("{0: <50}".format(':'.join(map(str, conc))))
        out.write('{0: <50}\n'.format("Total"))

        a_concs = pb.get_arrow_concs(params)
        conc_totals = [0 for i in range(len(cList))]
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

                    out.write("{0: <50d}".format(total_num))
                    total += total_num
                    conc_totals[iconc] += total_num
            out.write('{0: <10d}\n'.format(total))
        out.write("# " + ''.join(['-' for i in range(len(cList)*10 + 10 + 30)]) + '\n')
        out.write("{0: <28}".format("  0 0 0 0 0 0"))
        for ctotal in conc_totals:
            out.write("{0: <50d}".format(ctotal))
        out.write("{0: <50d}\n".format(sum(conc_totals)))
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

    keep_supers = False
    if args["super"]:
        keep_supers = True
    
    params = io.read_lattice(args["lattice"])
    systems = io.read_enum(args["input"])
    io.write_enum(params, outfile=args["outfile"])    

    count_t = 1
    count_s = 0
    from operator import itemgetter
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

        enumlist = []    
        for HNF, conc, num_wanted in systems:
            if not params["arrows"]:
                edata = datadicts[tuple(HNF)]
                SNF = edata["SNF"]
                LT = edata["L"]
            
                a_concs = pb.get_arrow_concs(params)
                configs = pb.enum_sys(edata["group"], list(conc), a_concs, num_wanted,HNF,params,
                                      keep_supers, args["acceptrate"])
            else:
                (SNF,L,R) = SmithNormalForm(get_full_HNF(HNF))
                SNF = [SNF[0][0],SNF[1][1],SNF[2][2]]
                LT = [item for row in L for item in row]
                a_concs = pb.get_arrow_concs(params)
                configs = pb.enum_sys(None, list(conc), a_concs, num_wanted,HNF,params, keep_supers, args["acceptrate"])

            for config in configs:
                labeling, arrowing = io.create_labeling(config)
                enumlist.append((sum(conc), HNF, SNF, LT, labeling, arrowing))
                

            sortenum = sorted(enumlist, key=itemgetter(0, 4))
            last_sz = sortenum[0][0]
            
        for (size, HNF, SNF, LT, labeling, arrowing) in sortenum:
            if size != last_sz:
                count_s = 1
                last_sz = size
            else:
                count_s += 1
            o = sfmt.format(count_t, 1, 1, 1, 1, count_s, size, 1, fmtn(SNF, 3), fmtn(HNF, 3),
                            fmtn(LT, 5), labeling, arrowing)
            f.write(o)
            count_t += 1

def examples():

    """Prints examples of using the script to the console using colored output.
    """
    script = "PHENUM: Finds the symmetrically unique arrangements of atoms in a lattice."
    explain = ("For simple 1D potentials such as the infinite square well, " 
               "kronig-penny, ect. This code produces a numerical solution "
               "using a bisis expansion.")
    contents = [("Find the Polya distribution",
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

    required = ("REQUIRED: the system to be enumerated in a `enum.in` file "
                "and a `lattice.in` file that contains the lattice data.")
    output = ("RETURNS: a file `struct_enum.out` for the full enumeration, "
              "an `enum.in` file if the -distribution flag is used, or a "
              "`polya.out` file for each cell size if the -polya flag is used.")
    details = ("For all the examples below, it is assumed that you have already "
                "compiled the modified enumlib code as described in the "
                "README or in some other manner obtained the HNFs (supercells) and "
                "their symmetry groups. You will then need to specify if you are "
                "obtaining the number of unique arrangements for each supercell and "
                "concentration allowed for your system or enumerating (finding) the "
                "desired number of configurations for each HNF and concentration. "
                "Additionally you way change the default input file names to ones of "
                "your own creation.")
    outputfmt = ("")

    msg.example(script, explain, contents, required, output, outputfmt, details)

script_options = {
    "-polya": dict(action="store_true",
                   help=("Predict the total number of unique superstructures for each "
                         "cell shape and possible concentration.")),
    "-enum": dict(action="store_true",
                  help=("Enumerate the number of superstructures of specific shape "
                        "as specified in the 'enum.in' file.")),
    "-lattice": dict(default="lattice.in", type=str,
                     help=("Override the default input file name: 'lattice.in' for "
                           "enumeration parameters.")),
    "-input": dict(default="enum.in", type=str,
                   help=("Override the default 'enum.in' file name.")),
    "-exec": dict(default="enum.x", type=str,
                  help=("Override the default 'enum.x' executable with an executable "
                        "name or path.")),
    "-distribution": dict(default=None, nargs= "+", type=int,
                        help=("Makes an enum.in file when the distribution is 'all'. Otherwise "
                              "the distribution is printed to the screen for the user. "
                              " The distributions are built from the results of the polya "
                              "run and contain the specified number of structures.")),
    "-savedist": dict(action="store_true",
                      help=("Makes an enum.in file for the distribution over 'size', "
                            "'shape', or 'conc'.")),
    "-dataformat": dict(default="cells.{}", type=str,
                        help=("Specify the default folder name for any cell size that contains "
                              "the matrices and groups generated by 'enum.x'. Format is: cells.{} "
                              "where {} is a placeholder for the integer cell size.")),
    "-cellsdir": dict(default="", type=int,
                      help=("Specify the path to the cells.{} type folders for the run. This "
                            "can be a relative or an absolute path.")),
    "-outfile": dict(default=None, type=str,
                     help=("Override the default output file names: 'polya.out' for "
                           "polya counting; 'enum.out' for structure enumeration; "
                           "'enum.in' for the distributions.")),
    "-sizes": dict(default= None, nargs=2, type=int,
                   help=("Specify the start and stop cell sizes over which to distribute "
                         "the structure selection, weighted by the Polya distribution.")),
    "-acceptrate": dict(default=None, type=float,
                        help=("For large numbers of unique arrangements, specify how often "
                              "one should be kept in the final list.")),
    "-profile": dict(action="store_true",
                     help=("Profile the code as it runs. Requires vprof (`pip install "
                           "vprof`). Argument is a combination of 'c', 'm' and 'h'. See "
                           "vprof docs.")),
    "-seed": dict(default=None, type=int,
                  help=("The integer seed for the random number generator.")),
    
    "-super": dict(action="store_true",
                   help=("Overrides the exclusion of the superperiodic cells from the output.")),
    "-filter": dict(default=None, nargs = 2,
                    help=("Applys a filter over the 'shape' or the 'conc' option of the "
                          "distributions. The first entry should specify the filter the "
                          "second should be the name of the file containing the desired "
                          "restrictions.")),
    "-visualize": dict(action="store_true",
                       help=("Makes plots of the shapes of the HNF's in the provided "
                             "enum.in style file.")),
    "-shapes": dict(action="store_true",
                    help=("The visualization excludes the atoms showing just the shapes "
                          "of the cells.")),
    "-show": dict(action="store_true",
                  help=("If present then the visualization will loop through an interactive "
                        "image for each HNF.")),
}
"""dict: default command-line arguments and their
    :meth:`argparse.ArgumentParser.add_argument` keyword arguments.
"""

def _parser_options():
    """Parses the options and arguments from the command line."""
    #We have two options: get some of the details from the config file,
    import argparse
    from phenum import base
    pdescr = "Numerical DFT code."
    parser = argparse.ArgumentParser(parents=[base.bparser], description=pdescr)
    for arg, options in script_options.items():
        parser.add_argument(arg, **options)
        
    args = base.exhandler(examples, parser)
    if args is None:
        return

    return args # pragma: no cover

def _script_enum(args, testmode=False):
    """Generates the 'polya.out' or 'enum.out' files depending on the script arguments.
    """
    from os import path, system
    from phenum.base import set_testmode
    set_testmode(testmode)
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
        if len(glob("{}{}".format(args["cellsdir"],args["dataformat"].split('.')[0]+'.*'))) < 1:
            from phenum.msg import err, warn
            from phenum.io_utils import which, read_lattice, write_struct_enum
            from os import system, chdir, getcwd
            warn("The input folders {}{} do not exist.".format(args["cellsdir"],args["dataformat"]))
            warn("Now running your enumlib executable to build the folders.")
            # We need to know where the current directory is so we can get back to it.
            initial_dir = getcwd()
            if args["cellsdir"] != "":
                system("cp {} {}.".format(args["lattice"],args["cellsdir"]))
                if args["enum"]:
                    system("cp {} {}.".format(args["input"],args["cellsdir"]))

                chdir(args["cellsdir"])
                if path.isfile("struct_enum.out"):
                    system("rm struct_enum.out")
                
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
            # Navigate back to the initial directory if we left it.
            if initial_dir != getcwd():
                chdir(initial_dir)

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
    if args["visualize"]:
        _plot_HNFs(args,testmode=testmode)
        
if __name__ == '__main__':
    args = _parser_options()
    
    if args["profile"]:
        from vprof import profiler
        profiler.run(_script_enum, args["profile"], args=(args,), host='localhost', port=8000)
    else:
        _script_enum(args)
