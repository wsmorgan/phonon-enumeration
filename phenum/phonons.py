"""Methods needed for the enumeration of arrows."""

#import neede modules
from numpy import array, dot
#Find concentration takes a 1D array, col, and returns a 1D array, 
#Concs, containing the number of times each element in col appears 
#in the array.
def _find_concentrations(col):
    """Finds the concentration of a given array of colors.

    Args:
        col (list): an integer array of the labeling.

    Returns:
        concs (list): The concentrations of the colors.
    """
    concs=[]

    tcol= list(col)
    for i, tcol_i in enumerate(tcol):
        if tcol_i != 0:
            concs.append(tcol.count(tcol_i))
            for j, tcol_j in enumerate(tcol):
                if tcol_j == tcol_i and i != j:
                    tcol[j] =0
            tcol[i] =0
    return(concs)

#Takes an array of colors and arrows and tells us how many different
#colors have arrows on them.
def how_many_arrows(tcol):
    """Determines the number of colors that have arrows on them.

    Args:
        tcol (list): A 2D array of the labeling that contains the colors
          and arrows for each site.

    Returns:
        arrows (int): The number of arrows in the system.
        n_species (int): The number of different atomic species with arrows.
        aconcs (list): The concentration of each arrow species.
    """
    arrows = 0
    
    species = []
    for i in tcol:
        if i[0] >= 0:
            arrows += 1
            if i[1] not in species:
                species.append(i[1])

    aconcs = _find_concentrations(tcol)
    return(arrows,len(species),aconcs)

#An algorithm for sorting the col array so that colors without arraws
#appear first in the list by concentration followed by the colors with
#arraws that are also sorted by concentration.
def _col_sort(col_list):
    """Sorts the labeling so that colors with arrows appear last in the
    list and so that the arrowed and non-arrowed colors are sorted
    from lowest to highest concentration.

    Args:
        col_list (list): A 2D integer array of the full labeling of the
          system.

    Returns:
        colt (list): A 2D integer array of the sorted colors.
    """

    col1 = []
    col2 = []
    colt=[]
    
    # seperate the arrays into colors with arrrows and colors without
    for i in col_list:
        if i[0] < 0:
            col1.append(i)
        elif i[0] >= 0:
            col2.append(i)

    # sort each array by concentration
    col1 = sorted(col1, key = col1.count)
    col2 = sorted(col2, key = col2.count)
    # put them back together again
    for i in col1:
        colt.append(i)
    for i in col2:
        colt.append(i)
    return(colt)

#this method takes a configuration and finds the unique ways of
#placing the arrows on that configuration.
#col is the initial configuration.
#agroup is the group operations with their effects on the arrows.
#dim is the number of arrow directions that are possible for this system.
def add_arrows(col,agroup,dim, translations, accept=None,nested=False,num_wanted=None, small=False, supers = False):
    """Finds the unique arrangements of arrows for a given configuration.

    Args:
        col (list): A 2D integer array of the initial labeling.
        dim (int): The number of directions the arrows can point.
        translations (list): The translations of the lattice.
        accept (float, optional): acceptance rate of configurations for large enumerations.

        agroup (list): The stabilizers for the colors only with the arrow
          permutations.

        nested (bool, optionl): Set to True if this is called from within the context
          of another progress bar.

        supers (bool, optional): True if we want to include super periodic 
          arrangements in the enumeration.

    Returns:
        arsurvivors (list): The unique arrow arrangements.
    """
    from random import random
    #Find out how many arrows there are in the col array
    (narrows,arrow_types,conc_w_arrows) = how_many_arrows(col)
    #largest_arrow is largest value any arrow could have.
    largest_arrow = [dim-1]*narrows
    #this is an array for storing the unique arrangements.
    arsurvivors = []
    len_col = len(col)

    maxpossible = _ahash(largest_arrow,dim)
    from phenum.msg import verbosity
    if verbosity is not None and verbosity >= 1 and not nested: #pragma: no cover
        from tqdm import tqdm 
        #If accept is specified, then we should get on average that many out.
        if num_wanted is not None:
            ntotal = num_wanted
        else:
            ntotal = maxpossible
        pbar = tqdm(total=ntotal)

    #this loop runs over every possible configuration of arrows
    #staring from the one with the smallest hash number and ending
    #with the one with the largest hash number to see if they are
    #unique. Each hash will be stored in orighash (original hash) then
    #compared to the permuted arrows hash (permhash) at the end.
    orighash = 0
    if not small:
        while orighash <= maxpossible:
            # for orighash in range(maxpossible+1):
            arrow_config = _ainvhash(orighash,narrows,dim)
            unique = True
            temp_coloring_with_arroms = [0]*len(col)
            #we need to cycle through each of the operations in the
            #permutation group that got passed in to see if they will turn
            #the arrow configuration into one we've already seen.
            for perm in agroup:
                l = 0
                #this loop constucts the coloring with the correct arrows
                #for the permutation group to be applied to.
                for z in range(len_col):
                    if col[z][0] < 0:
                        temp_coloring_with_arroms[z] = list(col[z])
                    elif col[z][0] >= 0:
                        temp_coloring_with_arroms[z] = list(col[z])
                        temp_coloring_with_arroms[z][0] = arrow_config[l]
                        l += 1

                new_arrow_config=[]
                arrow_perm = perm[1]
                color_perm = perm[0]
                new_coloring_with_arrows = []
                #here we use the permutations of the colors to permute the
                #colors for our configuration.
                new_coloring_with_arrows.extend([temp_coloring_with_arroms[color_perm[j]]
                                                 for j in range(len(temp_coloring_with_arroms))])
                #if there is an arrow on any lattice site it needs to
                #be updated. Otherwise just permute the colors.
                new_arrow_config.extend([arrow_perm[site[0]] for site in new_coloring_with_arrows
                                         if site[0] >= 0])
                #if the new configuration has a smaller hash number than
                #the current configuration then we have seen it before and
                #it is not unique.
                permhash = _ahash(new_arrow_config,dim)
                if permhash < orighash:
                    unique = False
                    break
                #if none of the permutations have had a smaller hash number
                #than the current configuration then the configuration is
                #unique and we need to add it to the list.
            if unique == True:
                #this loop make the coloring with the arrows that will be
                #added to the list of unique configurations.
                coloring_with_arrows = [0]*len(col)
                i = 0
                for z in range(len_col):
                    coloring_with_arrows[z] = list(col[z])
                    if col[z][0] >= 0:
                        coloring_with_arrows[z][0] = arrow_config[i]
                        i += 1

                # We need to check the superperiodic structures.
                if not supers:
                    for trans in translations:
                        action = trans[0]
                        permutation = trans[1]
                        if action != list(range(len(action))) or (permutation != list(range(len(permutation)))):
                            orig_sites = [i[1] for i in coloring_with_arrows]
                            arrow_sites = [i[0] for i in coloring_with_arrows]
                            perm_sites = [orig_sites[i] for i in action]
                            rot_arrows = [arrow_sites[i] for i in action]
                            perm_arrows = [permutation[site] for site in rot_arrows
                                           if site >= 0]
                            if orig_sites == perm_sites and perm_arrows == arrow_config:
                                unique = False
                                break
                if (accept is None or random() < accept) and unique == True:
                    arsurvivors.append(coloring_with_arrows)
                    if verbosity is not None and verbosity >= 1 and not nested: #pragma: no cover
                        pbar.update(1)
                    if num_wanted is not None and len(arsurvivors) == num_wanted:
                        break
            orighash += 1
    else:
        from random import randrange
        visited = []
        orighash = randrange(maxpossible)
        while len(arsurvivors) < num_wanted:
            while orighash in visited:
                orighash = randrange(maxpossible)

            visited.append(orighash)
            
            arrow_config = _ainvhash(orighash,narrows,dim)
            unique = True
            temp_coloring_with_arroms = [0]*len(col)
            #we need to cycle through each of the operations in the
            #permutation group that got passed in to see if they will turn
            #the arrow configuration into one we've already seen.
            for perm in agroup:
                l = 0
                #this loop constucts the coloring with the correct arrows
                #for the permutation group to be applied to.
                for z in range(len_col):
                    if col[z][0] < 0:
                        temp_coloring_with_arroms[z] = list(col[z])
                    elif col[z][0] >= 0:
                        temp_coloring_with_arroms[z] = list(col[z])
                        temp_coloring_with_arroms[z][0] = arrow_config[l]
                        l += 1

                new_arrow_config=[]
                arrow_perm = perm[1]
                color_perm = perm[0]
                new_coloring_with_arrows = []
                #here we use the permutations of the colors to permute the
                #colors for our configuration.
                new_coloring_with_arrows.extend([temp_coloring_with_arroms[color_perm[j]]
                                                 for j in range(len(temp_coloring_with_arroms))])
                #if there is an arrow on any lattice site it needs to
                #be updated. Otherwise just permute the colors.
                new_arrow_config.extend([arrow_perm[site[0]] for site in new_coloring_with_arrows
                                         if site[0] >= 0])
                #if the new configuration has a smaller hash number than
                #the current configuration then we have seen it before and
                #it is not unique.
                permhash = _ahash(new_arrow_config,dim)
                if permhash in visited and permhash != orighash: #pragma: no cover
                    #This happens rarely and so it's really hard to access in tests.
                    unique = False
                    break
                #if none of the permutations have had a smaller hash number
                #than the current configuration then the configuration is
                #unique and we need to add it to the list.
            if unique == True:
                #this loop make the coloring with the arrows that will be
                #added to the list of unique configurations.
                coloring_with_arrows = [0]*len(col)
                i = 0
                for z in range(len_col):
                    coloring_with_arrows[z] = list(col[z])
                    if col[z][0] >= 0:
                        coloring_with_arrows[z][0] = arrow_config[i]
                        i += 1

                # We need to check the superperiodic structures.
                if not supers:
                    for trans in translations:
                        action = trans[0]
                        permutation = trans[1]
                        orig_sites = [i[1] for i in coloring_with_arrows]
                        perm_sites = [orig_sites[i] for i in action]
                        perm_arrows = [permutation[site] for site in arrow_config
                                         if site >= 0]
                        if orig_sites == perm_sites and perm_arrows == arrow_config and action != list(range(len(action))) and permutation != list(range(len(permutation))): #pragma: no cover
                            #This happens rarely and so it's really hard to access in tests.
                            unique == False

                if unique == True:
                    arsurvivors.append(coloring_with_arrows)
                
                if verbosity is not None and verbosity >= 1 and not nested: #pragma: no cover
                    pbar.update(1)
                # if num_wanted is not None and len(arsurvivors) == num_wanted:
                #     break
                        
    if verbosity is not None and verbosity >= 1 and not nested: #pragma: no cover
        pbar.close()

    return(arsurvivors)

# arrow_concs is a method that returns the concentration string
# including the arrows
def arrow_concs(cList,aconcs):
    """Uses the concentrations of the atoms and the arrows to make a
    labeling for the system. It then sorts them to be in the correct
    order for the code.

    Args:
        cListr (list): An integer array the concentration of the colors.

        aconcs (list): An array of floats of the fractional number of arrows
          for each color

    Returns:
        decoration_w_arrows (list): The labeling of both colors and arrows.
    """
    aconcs = [int(cList[i]*aconcs[i]) for i in range(len(cList))]
    
    decoration_w_arrows = []
    labels = range(1,len(cList)+1)
    for i, ac_i in enumerate(aconcs):
        na = ac_i
        ns = cList[i]
        if ns >= na:
            while na > 0:
                decoration_w_arrows.append([1,labels[i]])
                na -= 1
                ns -= 1
            while ns > 0:
                decoration_w_arrows.append([-1,labels[i]])
                ns -= 1
        else:
            while ns > 0:#pragma: no cover
                #I've never been able to access this code without
                #something going wrong before this
                decoration_w_arrows.append([-1,labels[i]])
                ns -= 1

    decoration_w_arrows = _col_sort(decoration_w_arrows)
    
    return(decoration_w_arrows)

def get_arrow_concs(params):
    """If the concentrations are being restricted then find the correct 
    arrow for each species included.

    Args:
        params (dict): The parameters read in from lattice.in.

    Returns:
        a_concs (list): The arrow concentration for each species.
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

def enum_sys(groupfile, concs, a_concs, num_wanted, HNF, params, supers, accept=None):
    """Enumerates a random subset of the unique structures that have the shape
    defined by the symmetry group and the specified concentration.

    Args:
        groupfile (str): Path to the file containing the symmetry group.
        concs (list): Integer array of the concentrations for each species.
        a_concs (list): Integer array of the arrow concentrations for each species.
        HNF (list): The HNF matrix for the system we're currently enumerating
        params (dict): The dictionary of parameters read in from lattice.in
        supers (bool): True if superperiodic structures are to be kept.
        accept (float, optional): The acceptance rate for this enumeration.

        num_wanted (int): The number of structures to pick randomly from the enumerated
          list.

    Returns:
        configs (list): A list of the unique labelings in this system.

    Raises:
        ValueError: if the number of configurations found doesn't match the number requested.
    """

    from phenum.grouptheory import get_full_HNF, get_sym_group
    from phenum.polyaburnside import polya
    from phenum.tree import brancher, guess_and_check_brancher
    import phenum.io_utils as io

    cellsize = sum(concs)/len(params["basis_vecs"])
    decorations = arrow_concs(concs, a_concs)
    # get the symmetry group for this HNF. Assumes the group can be
    # found in the file labeled by (this_HNF)_sym_group.out
    if groupfile is  None:
        if sum(a_concs) == 0:
            arrows = False
        else:
            arrows = True

        sym_g = get_sym_group(params["lat_vecs"],params["basis_vecs"],
                              get_full_HNF(HNF),3,arrows=arrows)

        agroup = [list(sym_p) for sym_p in zip(sym_g.perm.site_perm,sym_g.perm.arrow_perm)]
        
    else:
        group = io.read_group(groupfile)
        # get symgroup from HNF and lat_vecs
        # add [0] to each element of the symmetry group
        agroup = [[g,[0]] for g in group]

    # we need to know the concentrations of the species with and
    # without arrows, we also need to know the number of arrows and
    # their species so we can undo the previous step later
    (n_arrows, arrow_types, sorted_concs) = how_many_arrows(decorations)

    # if we're enumerating a relatively large system (n>=10) but only
    # want a relatively small number of unique configurations (n<=100)
    # then we don't need to run the polya algorithm.
    if sum(concs) >=10 and num_wanted <= 100:
        total = 1e10
    else:
        # now find the number of unique arrangements using
        # polya
        if arrow_types != 0:
            total = polya(sorted_concs, agroup, arrowings=arrow_types)
        else:
            total = polya(concs, agroup)
    
        # generate the random subset to be used. WARNING! for phonon enumerations
        #we have seen values that are *50* digits long! If the number exceeds
        #1e9, we change the approach to randomization.
        from phenum.msg import warn, err
        if num_wanted < total:
            if total < 1e6:
                from random import shuffle
                subset = list(range(1, total+1)) 
                shuffle(subset)
                subset = subset[0:num_wanted]
            else: # pragma: no cover
                # only used for data sets to large for test cases.
                subset = num_wanted
        elif num_wanted == total:
            subset = []
        else:
            warn("number of configurations requested exceeds the number of "
                 "unique configurations available.")
            subset = []

    n_stabs = []
    # if we're doing a purely arrow enumeration then we don't need to
    # do the tree search but instead perform the final step of the
    # algorithm to find the possible unique displacements of the atoms
    from random import random
    if len(concs) == 1 and all(decorations) >=0:
        configs = []
        if float(num_wanted)/total < 0.001 and accept is None:
            small = True
            a_configs = add_arrows(decorations, agroup, 6, agroup[0:int(cellsize)], accept=accept,
                                   num_wanted=num_wanted, small=small,supers=supers)
        else:
            small = False
            a_configs = add_arrows(decorations, agroup, 6, agroup[0:int(cellsize)], accept=accept,
                                   num_wanted=num_wanted, supers=supers)
            
        count = 1
        for config in a_configs:
            if isinstance(subset, list) and count in subset and not small:
                configs.append(config)
            else:# accept is not None or small:
                configs.append(config)
            count += 1
    else:
        if sum(concs) >= 10 and num_wanted <= 100:
            configs = guess_and_check_brancher(sorted_concs, agroup, decorations, 6, supers, cellsize, num_wanted)
        else:
            configs = brancher(sorted_concs, agroup, decorations, 6, supers, cellsize, total, subset, accept)

    if len(configs) != num_wanted and not super: #pragma: no cover
        #I've never been able to trigger this error.
        raise ValueError("Warning the enumeration code returned {0} structures when {1} "
                         "were asked for. This should not happen. Please submit a bug "
                         "report on https://github.com/wsmorgan/phonon-enumeration "
                         "including your input files so that this error may be "
                         "corrected.".format(str(len(configs)),str(num_wanted)))
        if len(configs) == 0: #pragma: no cover
            #I've neven been able to trigger this error.
            exit()

    return configs

def _ahash(coloring,dim):
    """Produces a unique number for each possible configuration of
    arrows. This number is used to compare the order that the arrows
    occure in.

    Args:
        coloring (list): A 2D integer array of the full coloring with colors and arrows.
        dim (int): The number of directions the arrows can point.

    Returns:
        hash (int): The hash for the arrow configuration.
    """
    return sum([coloring[i]*dim**i for i in range(len(coloring))])

#anum is a unique number that is associated with an array of arrows.
#num_of_arrows in the number of arrows that are in the array.
#dim is the number of directions the arrows can point.
def _ainvhash(anum,num_of_arrows,dim):
    """Turns an arrow hash back into the array of arrow directions.

    Args:
        anum (int): The arrow hash number.
        num_of_arrows (int): The number of arrows in the system.
        dim (int): The number of directions the arrows can point.

    Returns:
        arrows (list): The arrow labeling.
    """
    arrows = [0]*num_of_arrows
    for i in range(num_of_arrows):
        base = dim**(num_of_arrows-1-i)
        arrows[num_of_arrows-1-i] = anum//base
        anum -= base*arrows[num_of_arrows-1-i]
    return(arrows)
