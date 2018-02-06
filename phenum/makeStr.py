#!/usr/bin/env python
from phenum import msg

def RepresentsInt(s):
    """Determines if a stri could be an int.

    Args:
        s (str): The string to be tested.
    
    Returns:
        True if the string is an integer.
    """    
    try: 
        int(s)
        return True
    except ValueError:
        return False

def _make_structures(args, return_euids = False):
    """Makes a VASP POSCAR file for the desired structures.

    Args:
        args (dict): The user input.
        return_euid (bool): a unique identifier for the structures in the
            POSCARS created.
    """
    from phenum.io_utils import read_enum_out, write_POSCAR, write_config
    from phenum.vector_utils import map_enumStr_to_real_space, cartesian2direct

    (system, structure_data) = read_enum_out(args)

    euids = []
    # for each structure write the vasp POSCAR
    for structure in structure_data:
        # space_data is a dictionary containing the spacial data for
        # the structure
        euid = "".join((str(i) for i in system["plattice"]))
        euid += "".join(str(i) for i in structure["HNF"])
        euid += structure["labeling"]
        euid = (euid,structure["directions"],str(args["rattle"]),str(args["displace"]))
        from hashlib import sha1 
        euids.append(sha1(''.join(euid).encode('utf-8')))
        space_data = map_enumStr_to_real_space(system,structure,args["mink"])

        space_data["aBas"] = cartesian2direct(space_data["sLV"],
                                              space_data["aBas"],system["eps"])

        if args["config"]=="f":
            write_POSCAR(system,space_data,structure,args)
        elif args["config"]=="t":
            write_config(system,space_data,structure,args,args["mapping"])
            
    if return_euids: #pragma: no cover
        return euids

def examples():
    """Print some examples on how to use this python version of makeStr."""
    script = "makeStr: Makes a vasp style POSCAR for the desired system."
    explain = ("For all the examples bellow it is assumed you have already "
               "run the enumeration code and produced an enum.out style file.")
    contents = [("Make a single POSCAR file",
            "To make a POSCAR file for a specific structure listed in the "
            "`enum.out` style file you will need to identify the structure \n number "
            "(the first number of each row in the file) for the structure you want "
            ". For example to make a POSCAR for structure number 10 \n from an `enum.out` "
            "file.","makeStr.py 10 \n"),
           ("Make multilpe POSCARS at once",
            "To make multiple POSCARS for a range of values in the `enum.out` style "
            "file simply list the starting and ending structure numbers \n of the range. "
            "To make POSCARS for every structure in the output file use the word `all`.",
            "makeStr.py 10 20 \n  makeStr.py all \n"),
            ("Find the lattice parameter for the system",
             "To have makeStr.py predict the lattice parameter for the system using "
             "Vegard's Law use the -species option followed by a space \n separated list "
             "of the elements in the system.","makeStr.py 10 -species Al Cu \n"),
            ("Include displacements in POSCAR",
             "If `arrows` (displacement directions) were included in the enumeration "
             "then it is possible to displace them off the lattice points \n when making the "
             "POSCARS using the -displace option followed by the displacement amount "
             "expressed in terms of the lattice parameter. \n In other words if `a` is "
             "the lattice parameter and the atoms were to be displaced by `a/2` then "
             "the command would be:","makeStr.py 10 -displace 0.5 \n"),
            ("Make displacements have different lengths in POSCAR",
             "If `arrows` were included in the model and the `-displace` flag is being "
             "used it is possible to 'rattle' the displacements so that \n they are not all "
             "the same length. Using the `-rattle` option applies a random distribution "
             "to the displacements with the larges change \n in the displacements specified "
             "by the user as a fraction of the displacement given. So if a displacement of "
             "0.5 was given and the \n displacements were to be randomized by 1/4 of that total "
             "displacement the the command would be:",
             "makeStr.py 10 -displace 0.5 -rattle 0.25"),
            ("Make an MTP config file for a ternary endge from a binary enumeration",
             "MTP is an excellent interatomic potential capable of quickly relaxing "
             "a cell as such we would like to be able to output the structures we've "
             "selected to an MTP compatible config file rather than a VASP POSCAR. "
             "There are also times when a binary enumeration is desired to form the edges "
             "of a ternary phase diagram, to do so the species present have to be able to "
             "be mapped to other species in use. To make a binary edge of a ternary phase "
             "diagram for an MTP config file use:",
             "makeStr.py 10 -config t -species Al Cu Ni -outfile=to-relax.cfg "
             "-species_mapping 1 2 (or 0 2)")]
    required = ("REQUIRED: A `enum.out` file.")
    output = ("RETURNS: A vasp style POSCAR labeled vasp.* where the `*` is replaced "
              "with the structure number for the `enum.out` file.")
    details = ("")
    outputfmt = ("")

    msg.example(script, explain, contents, required, output, outputfmt, details)

script_options = {
    "structures": dict(nargs="+",
                        help=("The desired structure numbers from the enum.out file. This "
                              "can be either a single value or a desired range indicated by "
                              "the starting and stopping structure numbers.")),
    "-displace": dict(default=0.0, type=float,
                        help=("The displacement amount for the arrows in units of the lattice "
                               "parameter. Default is 0.")),
    "-input": dict(default="struct_enum.out",type=str,
                        help=("Override the default 'struct_enum.out' file name.")),
    "-mink": dict(default="t", choices=["t","f"],
                        help=("Sets flag to perform minkowski reduction of the basis (T/F)."
                              " Default is True.")),
    "-species": dict(default=None, nargs="+",type=str,
                        help=("Specify the atomic species present in the system.")),
    "-species_mapping": dict(default=None, nargs="+",type=int,
                        help=("Specify the atomic species numbers for the system. This option "
                              "is only used for making MTP config files.")),
    "-outfile": dict(default="vasp.{}",type=str,
                        help=("Override the default output file names: 'vasp.{structure#}'" 
                              "for the structures.")),
    "-rattle": dict(default=0.0, type=float,
                        help=("Randomizes the positions of the atoms in the POSCAR by no "
                              "more than the fraction of the displacement provided.")),
    "-config" : dict(default="f",choices=["t","f"],
                   help=("make an MTP config file instead of a VASP POSCAR. If the "
                         "MTP config file already exists it will be appended to, not "
                         "overwritten.")),
    "-remove_zeros" : dict(default="f",choices=["t","f"],
                   help=("Remove the zeros from the concentrations string in the 'POSCAR'."))    
}
"""dict: default command-line arguments and their
    :meth:`argparse.ArgumentParser.add_argument` keyword arguments.
"""

def _parser_options():
    """Parses the options and arguments from the command line."""
    #We have two options: get some of the details from the config file,
    import argparse
    from phenum import base
    pdescr = "POSCAR contstruction."
    parser = argparse.ArgumentParser(parents=[base.bparser], description=pdescr)
    for arg, options in script_options.items():
        parser.add_argument(arg, **options)
        
    args = base.exhandler(examples, parser)
    if args is None:
        return

    return args #pragma: no cover

def run(args):
    """Generates the vasp output file for the desired structure.

    Args:
        args (dict): The user input.
    """

    if args is None:
        raise ValueError("Please enter a single structure number, two structures that "
                         "indicate the first and last structure to be used in the input "
                         "file, or all.")
    elif args["structures"] != None :
        if not RepresentsInt(args["structures"][0]) and args["structures"][0].lower() == "all":
            args["structures"] = None
        elif len(args["structures"])  == 1 and RepresentsInt(args["structures"][0]):
            args["structures"] = [int(args["structures"][0])]
        elif len(args["structures"]) == 2:
            args["structures"] = list(range(int(args["structures"][0]),
                                            int(args["structures"][1])+1))
        else:
            raise ValueError("Please enter a single structure number, two structures that "
                             "indicate the first and last structure to be used in the input "
                             "file, or all. The values {} don't match this "
                             "format.".format(args["structures"]))
    else:
        raise ValueError("Please enter a single structure number, two structures that "
                         "indicate the first and last structure to be used in the input "
                         "file, or all. The values {} don't match this "
                         "format.".format(args["structures"]))

    if args["species_mapping"] is not None:
        args["mapping"] = {}
        for i in range(len(args["species_mapping"])):
            args["mapping"][i] = args["species_mapping"][i]
    else:
        args["mapping"] = None
        
    if args["remove_zeros"] == "t":
        args["remove_zeros"] = True
    else:
        args["remove_zeors"] = False

    _make_structures(args)
        
if __name__ == '__main__':
    run(_parser_options())
