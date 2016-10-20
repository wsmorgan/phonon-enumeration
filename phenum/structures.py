"""Methods for reading and writing enumeration and polya-counting results."""

"""Generates a random subset of possible structures weighted by the
Polya distribution for superstructures.
"""
def enum_data(cellsize, args, params): #pragma: no cover
    """Returns a list of dictionaries with the HNF, SNF, Left Transform
    and permutation group file names for the given cell size.
    """
    from os import path
    from numpy import loadtxt, array
    dirname = "{}{}".format(args["cellsdir"],args["dataformat"].format(cellsize))
    result = []

    if path.isdir(dirname):
        matrices = loadtxt(path.join(dirname, "matrices"), int)
        limit = matrices.shape[0] if len(matrices.shape) > 1 else 1
        matrices = array([matrices]) if limit == 1 else matrices
        for i in range(limit):
            result.append({
                "SNF": matrices[i,1:4],
                "HNF": matrices[i,4:10],
                "L": matrices[i,10:],
                "group": path.join(dirname, "group.{}".format(matrices[i,0]))})
    else:
        from phenum.msg import warn, err
        from os import system
        from phenum.io_utils import write_struct_enum, which
        warn("No data for cell size {0:d} found at {1}. Creating new {2} files now using "
             "{3}.".format(cellsize, dirname,args["dataformat"],args["exec"]))
        if which(args["exec"]) != None:
            write_struct_enum(params)
            system(args["exec"])
            system('rm symops_enum_parent_lattice.out readcheck_enum.out fort.*')
            if path.isdir(dirname):
                matrices = loadtxt(path.join(dirname, "matrices"), int)
                limit = matrices.shape[0] if len(matrices.shape) > 1 else 1
                matrices = array([matrices]) if limit == 1 else matrices
                for i in range(limit):
                    result.append({
                        "SNF": matrices[i,1:4],
                        "HNF": matrices[i,4:10],
                        "L": matrices[i,10:],
                        "group": path.join(dirname, "group.{}".format(matrices[i,0]))})
            else:
                err("The executable {0} does not produce the needed input folders {1}, "
                    "or else your struct_enum.in doesn't agree with your {2} file. "
                    "\nPlease ensure that your struct_enum.in is correct. If it is then "
                    "please follow the instructions found to compile {0} in the "
                    "README.".format(args["exec"],args["dataformat"],args["lattice"]))
                exit()                            
        else:
            err("Could not find {} on your path. Please add it to your path "
                "if you have already\n compiled it. Otherwise please follow the "
                "instructions found in the README to download, make,\n and place the "
                "executable in your path.".format(args["exec"]))
            exit()
        
    return result   

def _read_struct_enum(): #pragma: no cover
    """Reads in the struct_enum.out file and returns a dictionary with key
    for each cell size, and value a dictionary keyed by HNF with the number
    of unique structures.
    """
    i = 0
    structs = {}
    with open("struct_enum.out") as f:
        for line in f:
            i += 1                        
            if i > 17:
                vals = line.split()
                sN = int(vals[6])
                HNF = tuple(map(int, vals[11:17]))
                concs = tuple([vals[-1].count(c) for c in ["0", "1"]])
                if sN not in structs:
                    structs[sN] = {}
                if HNF not in structs[sN]:
                    structs[sN][HNF] = {}
                    structs[sN][HNF][concs] = 1
                else:
                    if concs in structs[sN][HNF]:
                        structs[sN][HNF][concs] += 1
                    else:
                        structs[sN][HNF][concs] = 1
    return structs

def _write_struct_summary(structs): # pragma: no cover
    """Writes the summary of unique structure counts by HNF and cell size
    to file for verification of polya.
    """
    from numpy import zeros
    for sN, HNFS in list(structs.items()):
        out = open('enum_'+str(sN)+'.out', 'w+')
        out.write('{0: <28}'.format("# HNF"))
        first = sorted(next(iter(HNFS.values())))
        for conc in first:
            out.write("{0: <10}".format(':'.join(map(str, conc))))
        out.write('{0: <10}\n'.format("Total"))

        conc_totals = {conc: 0 for conc in first}
        sHNFs = sorted(HNFS.values)
        for HNF in sHNFs:
            out.write("  {0: <26}".format(' '.join(map(str, HNF))))
            for conc in first:
                if conc in dHNF:
                    out.write("{0: <10d}".format(dHNF[conc]))
                    conc_totals[conc] += dHNF[conc]
                else:
                    out.write("{0: <10d}".format(0))
            out.write('{0: <10d}\n'.format(sum(dHNF.values())))
            
        out.write("# " + ''.join(['-' for i in range(len(first)*10 + 10 + 30)]) + '\n')
        out.write("{0: <28}".format(""))
        for conc in first:
            out.write("{0: <10d}".format(conc_totals[conc]))
        out.write("{0: <10d}\n".format(sum(conc_totals.values())))
        out.close()

def _distribute(cellsizes, ftype, n=None, dataformat="cells.{}",seed=None, res_type=None, res_values=None):
    """Returns a dictionary specifying how many of each cell shape, size and concentration
    to enumerate in order to obtain a grand total of 'n' unique cells, distributed according
    to the abundance of unique structures predicted by Polya. See the calling signature for
    method:distribution_fxn(). However, while that method does not guarantee the total 'n',
    this one does; if it is short a few at the end, it tops up the list according to the
    abundances across all cell sizes, shapes, etc.

    :arg cellsizes: a list of the cell sizes to include in the data set summary.
    :arg ftype: one of ['shape', 'conc', 'size', 'all']; determines the number of 
      arguments that the distr. function will accept and average over, thus controlling
      the granularity of the predictions.
    :arg n: the number of structure types to return, weighted by their polya distribution.
      If this value is None, then *all* structure types are returned with their unique
      number.
    :arg dataformat: the name of the directories to search for 'polya.out' files in.
    :arg seed: The seed for the random number generator.
    """
    (f, dataset, gtotal) = _distribution(ftype, None, None, cellsizes=cellsizes, dataformat=dataformat,
                                         res_values=res_values,res_type=res_type)
    #If they didn't specify an 'n', then we also return all the structures.
    if n is None:
        n = gtotal
    elif n > gtotal:
        from .msg import warn
        warn("The number of unique structures you requested ({}) ".format(n) +
             "exceeds the total number of unique structures ({}). ".format(gtotal) +
             "The code will return *all* structures for the distribution.")
        n = gtotal

    result = {}
    rtotal = [0]
    def vinsert(key, value, rtotal, n, result):
        if rtotal[0] < n: #Only update if we are deficient in structures.
            if value > 0:
                if key in result:
                    result[key] += value
                else:
                    result[key] = value
                rtotal[0] += value
            return True
        else:
            return False   
    
    def assign(relvals, f, rtotal, gtotal, n, result, ftype,seed_val=None):
        """Recurses through the relative weights specified to fill the
        result dictionary until it has the exact number of structures
        requested.
        """
        if n == gtotal:
            for key, value, limit in relvals:
                vinsert(key, limit, rtotal, n, result)
        else:
            from numpy import round
            from math import modf
            from random import random, seed, getstate
            from operator import itemgetter
            import sys
            if seed_val != None:
                seed(a=seed_val)
                
            ids = {
                "all": lambda key: key,
                "shape": lambda key: key[0:2],
                "conc": lambda key: (key[0], key[2]),
                "size": lambda key: (key[0],)
            }
            #Sort the relative values according to their weights so that
            #the most abundant ones pick up the extras.
            relvals = sorted(relvals, key=itemgetter(1), reverse=True)
            recount = 0
            while rtotal[0] < n and recount < 4:
                delta = n-rtotal[0]
                for key, value, limit in relvals:
                    #If the abundance guarantees an integer part, use that, but then
                    #use a random number to decide who to keep for fractional parts.
                    #Eventually this should terminate unless all of the fractional
                    #pieces are super small; but even then, we should get what we want
                    #by the time we have iterated through the whole list.
                    rf, qf = modf(value if recount == 0 else f(delta, *(ids[ftype](key))))
                    retval = int(qf + (1 if random() < rf else 0))
                    if not vinsert(key, retval, rtotal, n, result):
                        break
                recount += 1

            if recount == 3 and rtotal[0] < n: #pragma: no cover
                from .msg import warn
                warn("Reached maximum recursion limit in random assignment.")
                
    #First, compile a dictionary at the lowest level that has the weights
    #for the desired value of 'n'.
    relvals = []
    if ftype=="all":
        for size, data in list(dataset.items()):
            for HNF, distr in list(data["distr"].items()):
                for conc, limit in list(distr.items()):
                    key = (size, HNF, conc)
                    value = f(n, size, HNF, conc)
                    relvals.append((key, value, limit))
    elif ftype == "shape":
        for size, data in list(dataset.items()):
            for HNF, limit in list(data["stotals"].items()):
                key = (size, HNF, None)
                value = f(n, size, HNF)
                relvals.append((key, value, limit))                
    elif ftype == "conc":
        for size, data in list(dataset.items()):
            for conc, limit in list(data["ctotals"].items()):
                key = (size, None, conc)
                value = f(n, size, conc)
                relvals.append((key, value, limit))
    elif ftype == "size":
        for size, data in list(dataset.items()):
            limit = data["gtotal"]
            key = (size, None, None)
            value = f(n, size)
            relvals.append((key, value, limit))
    
    #This makes the selection according to the relative values and keeps
    #choosing, weighted by relative abundance, until we have the right
    #number of structures.
    assign(relvals, f, rtotal, gtotal, n, result, ftype,seed_val=seed)
    #Make sure we are returning exactly how many they asked for.
    if rtotal[0] > n: #pragma: no cover
        from .msg import warn
        warn("More structures were returned than asked for, should not be possible.")

    from operator import itemgetter
    return result

def _distribution(ftype, dataset, gtotal, cast=float, cellsizes=None, dataformat="cells.{}",res_type = None,
                  res_values=None):
    """Returns the Polya-weighted distribution function for the specified cell sizes.
    If dataset is None, a dataset and total are loaded using the specified cell sizes.
    This dataset and the total number of structures it contains are returned as
    arguments 2 and 3 of the tuple. Depending on the value of 'ftype', it returns 
    a tuple (f, dataset, gtotal) where 'f' is a function function with arguments:

    'all':   f(n ,size, HNF, conc)
    'shape': f(n, size, HNF)
    'conc':  f(n, size, conc)
    'size':  f(n, size)

    where 'size' is the integer cell size; 'HNF' is the array of 6-elements defining
    the HNF matrix that is the cell shape; 'conc' is a list of concentrations such as
    [1,5] for a 6-atom cell; 'n' is the *grand* total number of structures desired from
    calling the function over all possible values of the various parameters passed to
    the distribution function.

    The value returned by the function is an integer specifying how many of that
    structure type should be selected according to the total distribution across
    all structures included in cellsizes.

    :arg dataset: the dictionary returned by method:_distribution_summary().
    :arg gtotal: the total number of unique structures corresponding to the
      dictionary 'dataset'.
    :arg ftype: one of ['shape', 'conc', 'size', 'all']; determines the number of 
      arguments that the distr. function will accept and average over, thus controlling
      the granularity of the predictions.
    :arg cellsizes: a list of the cell sizes to include in the data set summary.
    :arg dataformat: the name of the directories to search for 'polya.out' files in.
    :arg res_type: 'shape' if shapes are being restricted, or 'conc' if concentrations are  being restricted.
    :arg res_values: The allowed values of the restricted parameter.
    """
    if dataset is None:
        if cellsizes is not None:
            if res_type is not None:
                if res_type == "shape":
                    (dataset, gtotal) = _distribution_summary(cellsizes, dataformat=dataformat, HNFs=res_values)
                elif res_type == "conc":
                    (dataset, gtotal) = _distribution_summary(cellsizes, dataformat=dataformat, wanted_concs=res_values)
                else:
                    raise ValueError("Cannot filter the distribution using {}. Please use 'shape' or 'conc'.".format(res_type))

            else:
                (dataset, gtotal) = _distribution_summary(cellsizes, dataformat=dataformat)
        else:
            from .msg import err
            err("No dataset or cell sizes specified.")
            return None
        
    #We just need to sum up the total number of unique structures across all the
    #given cell sizes and then create a function that weights a total number of
    #desired structures by concentration and cell shape.
    ftotal = float(gtotal)
    if ftype == "all":
        f = lambda n, size, HNF, conc: cast(dataset[size]["distr"][tuple(HNF)][conc]/ftotal*n)
    elif ftype == "shape":
        f = lambda n, size, HNF: cast(dataset[size]["stotals"][tuple(HNF)]/ftotal*n)
    elif ftype == "conc":
        f = lambda n, size, conc: cast(dataset[size]["ctotals"][tuple(conc)]/ftotal*n)
    elif ftype == "size":
        f = lambda n, size: cast(dataset[size]["gtotal"]/ftotal*n)
    else:
        raise ValueError("The parameter {} is not a valid parameter for the distribution. "
                         "Please use size, shape, conc, or all.".format(ftype))

    return (f, dataset, gtotal)
        
def _distribution_summary(cellsizes, HNFs = None, wanted_concs = None, dataformat="cells.{}"):
    """Returns a dictionary that summarizes the Polya predictions for unique structures
    by cell size, shape and concentration for use by a distribution function. Returns
    a tuple (summary dict, grand total), where grand total is the total number of structures
    predicted for *all* cell sizes.

    :arg cellsizes: a list of the cell sizes to include in the data set summary.
    :arg dataformat: the name of the directories to search for 'polya.out' files in.
    :arg HNFs: a list of the allowed HNFs.
    :arg wanted_concs: a list of the allowed concentrations.
    """
    from os import path
    from numpy import loadtxt
    gtotal = 0
    dataset = {}
    for s in cellsizes:
        dirname = dataformat.format(s)
        source = path.join(dirname, "polya.out")
        if not path.isfile(source):
            from .msg import err
            err("Cannot find polya distribution for size {} at {}".format(s, source))
            continue

        #First load the possible concentrations from the commented first line.
        with open(source) as f:
            headings = f.readline().split()
            concs = [tuple(map(int, h.split(":"))) for h in headings[2:-1]]
        polya = loadtxt(source, int)
        distr = {}

        if HNFs == None and wanted_concs == None:
            for iHNF, HNF in enumerate(polya[0:-1,0:6]):
                distr[tuple(HNF)] = {tuple(c): v for c, v in zip(concs, polya[iHNF,6:-1])}
            stotals = {tuple(m): v for m, v in zip(polya[0:-1,0:6], polya[:-1,-1])}
            ctotals = {c: v for c, v in zip(concs, polya[-1, 6:-1])}
        elif HNFs != None:
            ctotals = {c: v for c, v in zip(concs,[0]*len(concs))}
            stotals = {}
            for iHNF, HNF in enumerate(polya[0:-1,0:6]):
                if tuple(HNF) in HNFs:
                    distr[tuple(HNF)] = {tuple(c): v for c, v in zip(concs, polya[iHNF,6:-1])}
                    for c, v in zip(concs, polya[iHNF,6:-1]):
                        ctotals[c] += v
            for m, v in zip(polya[0:-1,0:6], polya[:-1,-1]):
                if tuple(m) in HNFs:
                    stotals[tuple(m)] = v
        elif wanted_concs != None:
            ctotals = {c: v for c, v in zip(concs,[0]*len(concs))}
            stotals = {}
            for iHNF, HNF in enumerate(polya[0:-1,0:6]):
                cs = []
                vs = []
                count = 0

                for c, v in zip(concs, polya[iHNF,6:-1]):
                    if c in wanted_concs:
                        ctotals[c] += v
                        cs.append(c)
                        vs.append(v)
                        count += 1

                distr[tuple(HNF)] = {tuple(c): v for c, v in zip(cs, vs)}
            for m, v in zip(polya[0:-1,0:6], polya[:-1,-1]):
                if tuple(m) in distr:
                    stotals[tuple(m)] = v

        dataset[s] = {
            "distr": distr,
            "concs": concs,
            "ctotals": ctotals,
            "stotals": stotals,
            "gtotal": sum([ctotals[key] for key in ctotals])
        }
        gtotal += dataset[s]["gtotal"]

    return (dataset, gtotal)    

def _print_distribution(distr, distribution, filename=None, header=True, append=False, show=False):
    """Prints the specified distribution to screen or file.

    :arg distr: the distribution returned by method:distribute().
    :arg distribution: The type of distribution, i.e., 'shape', 'size', 'conc', or 'all'.
    :arg filename: The output file name.
    :arg header: True if the header is to be included in the file.
    :arg append: True if the file is to be appended to.
    :arg show: True if the distribution is to be printed to the screen.
    """
    if show:
        from .msg import arb, cenum
        bysize = {}
        sfmt = " {0: <5d} | {1: <15} | {2: <5} | {3: <5d} "
        for key, value in list(distr.items()):
            size, HNF, conc = key
            skey = sfmt.format(size, "" if HNF is None else ' '.join(map(str, HNF)),
                               "" if conc is None else ':'.join(map(str, conc)), value)
            cols = (cenum["cwarn"], cenum["cinfo"], cenum["cgens"], cenum["cokay"])
            if (size, value) in bysize:
                bysize[(size, value)].append((skey, cols, tuple(HNF) if HNF is not None else None))
            else:
                bysize[(size, value)] = [(skey, cols, tuple(HNF) if HNF is not None else None)]

        from operator import itemgetter
        for size, value in sorted(bysize.keys()):
            for skey, cols, HNF in sorted([[(i or "") for i in x] for x in bysize[(size, value)]], key=itemgetter(2)):
                arb(skey, cols, "|")

    if filename is not None:
        #We don't worry about ordering it, just write them to file in whatever
        #order the keys are in.
        with open(filename, 'w' if not append else 'a') as f:
            if header:
                if distribution == "all":
                    f.write("# {0: <28}  {1: <10}  {2}\n".format("HNF", "Conc.", "Number"))
                elif distribution == "shape":
                    f.write("# {0: <18}  {1}\n".format("HNF", "Number"))
                elif distribution == "conc":
                    f.write("# {0: <6}  {1: <6}  {2}\n".format("Size", "Conc.", "Number"))
                else:
                    f.write("# {0: <6}  {1}\n".format("Size", "Number"))
            for key, value in list(distr.items()):
                size, HNF, conc = key
                if conc == None:
                    conc = []
                if distribution == "all":
                    f.write("  {0: <28}  {1: <10}  {2:d}\n".format(' '.join(map(str, HNF)), ' '.join(map(str, conc)), value))
                elif distribution == "shape":
                    f.write("  {0: <18}  {1}\n".format(' '.join(map(str, HNF)), value))
                elif distribution == "conc":
                    f.write("  {0: <6}  {1: <6}  {2:d}\n".format(size, ' '.join(map(str, conc)), value))
                else:
                    f.write("  {0: <6}  {1:d}\n".format(size, value))
                    
def make_enum_in(distribution,directory,outfile,number=None,dataformat="cells.{}",sizes=None,save=True,seed=None,
                 restrict=None):
    """Makes an enum.in file if the distrubiton type is all with the
    desired number of structures. Otherwise prints the distribution
    information to the screen for the user.

    :arg distribution: The parameters that the distribution is over
    ('shape', 'conc', 'size', 'all').
    :arg n: The number of structures or 'all'.
    :arg dataformat: The folder name for the cell sizes.
    :arg sizes: when specified, limit the distribution to these integer cell sizes;
      otherwise, look for all cell sizes we have data for.
    :arg directory: The directory that contains the folders with the cell sizes.
    :arg outfile: The name of the output file for the distribution
    :arg save: True if the data is to be saved to file.
    :arg seed: The seed for the random number generator.
    :arg restrict: A list containing the restriction type ('shape','conc') and the file of allowed values.
    """

    from os import listdir, chdir, getcwd

    initial_directory = getcwd()
    if directory != "":
        chdir(directory)
    
    files = listdir(".")
    if sizes is None:
        sizes = []
        form = dataformat.split(".")[0]

        ready = False
        for f in files:
            if form in f:
                sizes.append(int(f.split(".")[1]))
                ready = True
    else:
        for i in sizes:
            celldir = dataformat.format(i)
            if celldir not in files:
                raise ValueError("Cannot find '{}' in current directory.".format(celldir))
        ready = True
                
    if ready == False:
        raise ValueError("The files {} don't exist in this directory. You either need to run"
            " the -polya option or else navigate into the folder that contains "
            "the output {} folders.".format(dataformat))

    if restrict is not None:
        # from numpy import loadtxt
        res_values = []
        with open(restrict[1],"r") as resf:
            for values in resf:
                if "#" in values:
                    pass
                elif restrict[0] == "shape":
                    res_values.append(tuple([int(i) for i in values.split()[0:6]]))
                elif restrict[0] == "conc":
                    res_values.append(tuple([int(i) for i in values.split()[1:3]]))

        distr = _distribute(sizes,distribution,n=number,dataformat=dataformat,seed=seed,res_type=restrict[0],res_values=res_values)        
    else:
        distr = _distribute(sizes,distribution,n=number,dataformat=dataformat,seed=seed)

    if initial_directory != getcwd():
        chdir(initial_directory)
        
    if distribution.lower() == "all":
        _print_distribution(distr,distribution,filename=outfile)
    else:
        _print_distribution(distr,distribution,filename= outfile if save else None,show=True)
