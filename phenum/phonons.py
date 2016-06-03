"""Methods needed for the enumeration of arrows."""

#import neede modules
from copy import deepcopy

#Find concentration takes a 1D array, col, and returns a 1D array, 
#Concs, containing the number of times each element in col appears 
#in the array.
def _find_concentrations(col):
    """Finds the concentration of a given array of colors.

    :arg col: an integer array of the labeling
    """
    Concs=[]

    tcol=deepcopy(col)
    for i in range(len(tcol)):
        if tcol[i] != 0:
            Concs.append(tcol.count(tcol[i]))
            for j in range(len(tcol)):
                if tcol[j] == tcol[i] and i != j:
                    tcol[j] =0
            tcol[i] =0
    return(Concs)

#Takes an array of colors and arrows and tells us how many different
#colors have arrows on them.
def how_many_arrows(tcol):
    """Determines the number of colors that have arrows on them.

    :arg tcol: a 2D array of the labeling that contains the colors
    and arrows for each site.
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

    :arg col_list: a 2D integer array of the full labeling of the
    system
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
    col1 = sorted(col1, key = col1.count, reverse=True)
    col2 = sorted(col2, key = col2.count, reverse=True)

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
def add_arrows(col,agroup,dim):
    """Finds the unique arrangements of arrows for a given configuration.

    :arg col: a 2D integer array of the initial labeling
    :arg agroup: the stabilizers for the colors only with the arrow
    permutations
    :arg dim: the number of directions the arrows can point
    """

    with open("add_arrows_group.out","w+") as af:
        for i in agroup:
            af.write(str(i) + "\n")
    
    #survivors is the array that contains the end result of the permutations
    survivor = []

    #Find out how many arrows there are in the col array
    (narrows,arrow_types,conc_w_arrows) = how_many_arrows(col)
    #largest_arrow is largest value any arrow could have.
    largest_arrow = [dim-1]*narrows
    #this is an array for storing the unique arrangements.
    arsurvivors = []

    #this loop runs over every possible configuration of arrows
    #staring from the one with the smallest hash number and ending
    #with the one with the largest hash number to see if they are
    #unique. Each hash will be stored in orighash (original hash) then
    #compared to the permuted arrows hash (permhash) at the end.
    for orighash in range(_ahash(largest_arrow,dim)+1):
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
            for z in range(len(col)):
                if col[z][0] < 0:
                    temp_coloring_with_arroms[z] = deepcopy(col[z])
                elif col[z][0] >= 0:
                    temp_coloring_with_arroms[z] = deepcopy(col[z])
                    temp_coloring_with_arroms[z][0] = arrow_config[l]
                    l += 1

            new_arrow_config=[]
            arrow_perm = perm[1]
            color_perm = perm[0]
            new_coloring_with_arrows = []
            #here we use the permutations of the colors to permute the
            #colors for our configuration.
            for j in range(len(temp_coloring_with_arroms)):
                new_coloring_with_arrows.append(temp_coloring_with_arroms[color_perm[j]])
            #if there is an arrow on any lattice site it needs to
            #be updated. Otherwise just permute the colors.
            for site in new_coloring_with_arrows:
                if site[0] >= 0:
                    site[0] = arrow_perm[site[0]]
                    new_arrow_config.append(site[0])
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
            for z in range(len(col)):
                if col[z][0] < 0:
                    coloring_with_arrows[z] = deepcopy(col[z])
                elif col[z][0] >= 0:
                    coloring_with_arrows[z] = deepcopy(col[z])
                    coloring_with_arrows[z][0] = arrow_config[i]
                    i += 1
                    
            arsurvivors.append(coloring_with_arrows)

    return(arsurvivors)

# arrow_concs is a method that returns the concentration string
# including the arrows
def arrow_concs(cList,aconcs):
    """Uses the concentrations of the atoms and the arrows to make a
    labeling for the system. It then sorts them to be in the correct
    order for the code.

    :arg cListr: an integer array the concentration of the colors 
    :arg aconcs: an array of floats of the fractional number of arrows
    for each color
    """

    aconcs = [int(cList[i]*aconcs[i]) for i in range(len(cList))]
    
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

    conc_w_arrows = _col_sort(conc_w_arrows)

    return(conc_w_arrows)

def get_arrow_concs(params):
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

def enum_sys(groupfile, concs, a_concs, num_wanted, HNF, params):
    """Enumerates a random subset of the unique structures that have the shape
    defined by the symmetry group and the specified concentration.

    :arg groupfile: path to the file containing the symmetry group.
    :arg concs: list of integer concentrations for each species.
    :arg a_concs: list of integer *arrow* concentrations for each species.
    :arg num_wanted: the number of structures to pick randomly from the enumerated
      list.
    :args HNF: the HNF for the system we're currently enumerating
    :args params: the dictionary of parameters read in from lattice.in
    """

    from phenum.grouptheory import get_full_HNF, get_sym_group
    from phenum.polyaburnside import polya
    from phenum.tree import brancher
    import phenum.io_utils as io

    decorations = arrow_concs(concs, a_concs)
    # get the symmetry group for this HNF. Assumes the group can be
    # found in the file labeled by (this_HNF)_sym_group.out
    if groupfile == None:
        sym_g = get_sym_group(params["lat_vecs"],params["basis_vecs"],get_full_HNF(HNF),3)
        agroup = []
        for i in range(len(sym_g.perm.site_perm)):
            agroup.append([sym_g.perm.site_perm[i],sym_g.perm.arrow_perm[i]])
    else:
        group = io.read_group(groupfile)
        # get symgroup from HNF and lat_vecs
        # add [0] to each element of the symmetry group
        agroup = [[g,[0]] for g in group]

    # we need to know the concentrations of the species with and
    # without arrows, we also need to know the number of arrows and
    # their species so we can undo the previous step later
    (n_arrows, arrow_types,concs_w_arrows) = how_many_arrows(decorations)

    # now find the number of unique arrangements using
    # polya
    if arrow_types != 0:
        # Since enumlib doesn't write the arrow group out we have to
        # recompute the group actions paired with their effects on the
        # arrows
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
        from phenum.msg import warn
        warn("number of configurations requested exceeds the number of "
            "unique configurations available.")
        subset = []

    n_stabs = []
    # if we're doing a purely arrow enumeration then we don't need to
    # do the tree search but instead perform the final step of the
    # algorithm to find the possible unique displacements of the atoms
    if len(concs) == 1 and all(decorations) >=0:
        configs = []
        a_configs = add_arrows(decorations, agroup, 6)
        count = 1
        for config in a_configs:
            if count in subset:
                configs.append(config)
            count += 1
    else:
        if arrow_types != 0:
            configs = brancher(concs_w_arrows, agroup, decorations, 6, subset)
        else:
            configs = brancher(concs, agroup, decorations, 6, subset)

    # reduced_configs is the list of configurations with the
    # superperiodic configurations removed.
    reduced_configs = []
    # need to find a way to remove the superperiodic arrangements
    reduced_configs = configs                
    return reduced_configs

def _ahash(coloring,dim):
    """Produces a unique number for each possible configuration of
    arrows. This number is used to compare the order that the arrows
    occure in.

    :arg coloring: is a 2D integer array of the full coloring with colors and arrows	    
    :arg dim: the number of directions the arrows can point
    """
	
    narrow = 0
    for i in range(len(coloring)):
	narrow = narrow + coloring[i]*dim**i
    return(narrow)

#anum is a unique number that is associated with an array of arrows.
#num_of_arrows in the number of arrows that are in the array.
#dim is the number of directions the arrows can point.
def _ainvhash(anum,num_of_arrrows,dim):
    """Turns an arrow hash back into the array of arrow directions.

    :arg anum: the arrow hash number
    :arg num_of_arrows: the number of arrows in the system
    :arg dim: the number of directions the arrows can point
    """
    arrows = [0]*num_of_arrrows
    for i in range(num_of_arrrows):
	base = dim**(num_of_arrrows-1-i)
	arrows[num_of_arrrows-1-i] = anum//base
	anum = anum -base*(anum//base)
    return(arrows)
