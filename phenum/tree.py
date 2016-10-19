"""The heart of the code is here. These methods build the tree and
check the canfigurations for uniqueness."""

def _coefficients(number_of_each_color,number_of_slots):
    """This method generates the radix numbers, or total number of
    combinations of each color, for a given aray.

    :arg number_of_each_color: is an integer array integer indicating the
    number of times each color is found in the array.
    :arg number_of_slots: an integer equal to the sum of the
    number_of_each_color array.
    """
    from .numerics import binomial_coefficient

    #renames variables for convience and creates an empty array
    coe = []
    for i in range(0,len(number_of_each_color)):
    #Uses the binomial_coefficient method to build an array of radix numebrs
        coe.append(binomial_coefficient(number_of_slots,number_of_each_color[i]))
        number_of_slots -= number_of_each_color[i]
    return coe

def _hash(listi,color):
    """Generates a hash for each color in the ladeling returning a array
    of the hash for each color.
    
    :arg listi: An integer array that contains the labeling.
    :arg color: An integer array of the concentrations of the
    colors.

    For details on this method see: http://msg.byu.edu/papers/enum3.pdf
    """
    from .numerics import binomial_coefficient

    m = len(listi)
    li = listi[::-1]
    rm1 = m-li.index(1)
    y = 0
    z = listi[:rm1].count(0)#m - listi[:rm1].count(1)
    li = listi
    for i in range(z):
        p0 = li.index(0) #next(x[0] for x in enumerate(li) if x[1] > j)
        k = li[p0:].count(1)
        li = li[p0+1:]
        y += binomial_coefficient(m-(p0+1),k-1)
        m -= (p0 + 1)

    return y

def _invhash(branch,colors,n):
    """Turns a hash array into a labeling, i.e., undoes the hash()
    subroutine.
    
    :arg branch: an integer list of the radix number/hash array
    :arg colors: an integer array of the concentration of the colors
    :arg n: the length of the labeling array, i.e., the
    number of sites in the system

    For the details on this method see: http://msg.byu.edu/papers/enum3.pdf
    """

    from .numerics import binomial_coefficient

    new = [0 for y in range(n)]
    coluse = 1
    for i in range(0,len(branch)):
        newi = [0 for y in range(n)]
        I = branch[i]
        t = colors[i]
        for l in range(n,0,-1):
            binom = binomial_coefficient(l-1, t-1)
            if binom <= I:
                I -= binom
            else:
                newi[n-l] = coluse
                t -= 1
        coluse += 1
        n -= colors[i]
        si = 0
        for s in range(0,len(new)):
            if new[s] == 0:
                new[s] = newi[si]
                si += 1
    return(new)

#color_list takes the configuration and returns a unique list of the
#colors used without duplicates.
def _color_list(col):
    """Finds the unique colors in a labeling.

    :arg col: an integer array of the labeling.
    """
    colors=[]

    tcol=col[:]
    for i in range(len(tcol)):
        if tcol[i] != 0:
            colors.append(tcol[i])
            for j in range(len(tcol)):
                if tcol[j] == tcol[i] and i != j:
                    tcol[j] =0
            tcol[i] =0
    return(colors)



def _perm(casei,colors,length,index,gen,stab,order,ast):
    """This method applies a cyclical permutation to an array to determine
    if it is unique. 
		
    :arg casei: the integer array hash for the current branch of
    the tree that will be permuted
    :arg colors: an integer array indicates the concentration of
    each color in the system
    :arg ast: the stabalizer for the final configuration that will
    be passed to the arrow permutation routine.
    :arg length: the total number of sites (colors) in the system
    :arg index: an integer that indicates which level of the tree
    we are in, i.e., which color is being added to the tree
    :arg gen: a 3D integer array that contains the permutation group
    :arg stab: a 3D integer array of the stabilizers for
    the previous level of this branch in the tree
    :arg order: is an integer array that keeps track of which
    configurations of the first level of the tree have already
    been seen

    The routin outputs the following:

    unique: an integer that indicates if the configuration is
    unique, a returned value of 0 is a unique configuration
    anything else is not
    st: a 2D integer array that stores the stabilizers for this
    configuration		   
    order and ast
    """
		
    #initialize unique to be 0 and arrays to be empty
    unique = 0
    if stab == 0:
        stab = []
    st = []
    ast = []
		
    # uses the invhash method from Inverse_radix_num to turn the id
    # number, casei, into an array
    li = list(_invhash(casei,colors,length))

    # if we're on the first level of the tree then we need to
    # follow a special procedure
    if index == 0:
	# first apply all the symmetry operations one at a time, a
	# is a single operation
        for a in gen:
	    # get the effect of the group on the colors, ingore the arrows
	    # i is the permutation of the colors for this operation
            i = a[0]
            lnew = []
	    # apply the permutation to get an equivalent configuration
            for j in i:
                lnew.append(li[j])
            lnewb = list(lnew)
                
	    # any color that was placed that isn't the color for the
	    # current level is set to zero
            for y in range(0,len(lnew)):
                if lnew[y] != index + 1:
                    lnewb[y] = 0
                elif lnew[y] == index + 1:
                    lnewb[y] = 1

	    # hash from radix_num_generator turns the array back into
	    # an integer hash array
            rtest = _hash(list(lnewb),colors)

	    # if the original hash is smaller than the new one we
	    # don't want to visit the new one again so we set its
	    # location in order to be -1
            if casei[0] < rtest:
                order[rtest] = -1
            
	    # if the original hash is larger than the new one then
	    # this configuration isn't unique so we set unique to be 1
	    # and reset the other variables
            elif casei[0] > rtest:
                unique += 1
                st = []
                ast = []
                break

	    # if the symmetry group didn't change this level
            # of the hash then it is a stabilizer and needs to
            # be saved
            elif casei[0] == rtest:
                st.append(a)

	    # if the symmetry op changed nothing in the
	    # configuration then it is a stabilizer for the
	    # arrow level and needs to be saved
            if li == lnew:
                ast.append(a)

    # if there is only a single stabilizer then the configuration is
    # unique and the stabilizer needs to be used for any later braches
    elif len(stab) == 1:
        st = stab
        ast = stab

    # if there are multiple stabilizers then we need to apply each of
    # them just like we needed to apply the entire group on the first
    # level of the tree
    elif len(stab) > 1:
        ast = []
        # first apply all the stabilizers one at a time, a is a single
        # operation
        for a in stab:
	    # get the effect of the group on the colors, ingore the
	    # arrows i is the permutation of the colors for this
	    # operation
            i = a[0]
            lnew = []
	    # apply the permutation to get an equivalent configuration
            for j in i:
                lnew.append(li[j])
                lib= list(li)
            lnewb = list(lnew)
            lnewb2 = []
	    # any color that was placed that isn't the color for
	    # the current level is set to zero
            for y in range(0,len(li)):
                if li[y] > index + 1:
                    lib[y] = 0
            for y in range(0,len(lnew)):
                if lnew[y] > index + 1:
                    lnewb[y] = 0
            for y in range(0,len(lnew)):
                if lnewb[y] == 0:
                    lnewb2.append(0)
                elif lnewb[y] == index + 1:
                    lnewb2.append(1)
                    
	    # if the the new and original configurations are the
	    # same then save the stabilizers
            if lnewb == lib:
                st.append(a)
            if li == lnew:
                ast.append(a)

	    # hash from radix_num_generator turns the array back
	    # into an integer hash array
            rtest = _hash(list(lnewb2),colors)
	    # if the original hash is larger than the new one then
	    # this configuration isn't unique so we set unique to
	    # be 1 and reset the other variables
            if casei[index] > rtest:
                unique = 1
                st = []
                ast = []
                break
            
	    # if the configuration isn't unique break from the
            # loop
            if unique == 1:
                ast = []
                st = []
                break
    else:
        unique = 1
		
    return(unique,st,order,ast)


def brancher(concs,group,colors_w_arrows, dim, supers, cellsize, total=0, subset=None, accept=None):
    """This routine navigates the tree and saves the unique configurations
    to an array survivors.
		
    :arg concs: the concentrations of the colors or atoms
    :arg group: the symmetry group for the system, including the
    effect on the arrows (displacement directions)
    :arg colors_w_arrows: an integer 2D array that indiciates
    which atoms are being displaced, i.e., where the arrows are.
    :arg dim: the number of directions the arrows can point
    :arg total: the total number predicted by polya.
    :arg subset: an integer array of the subset of unique
    arrangements wanted
    :arg accept: for large enumerations, how often to accept configurations.
    :arg supers: Logical that indicates if super periodic structures are to be kept.
    :arg cellsize: The number of cells in the system 
	
    The method returns the list of unique configurations and the
    number of stabilizers for the last level of the tree.

    survivors: is the 3D array of the unique arrangements of
    colors and arrows
    stob_len: the number of stabilizers for a given unique
    configuration
    """

    from .phonons import how_many_arrows, add_arrows
    from copy import deepcopy

    # initial setup	

    # redifine the colors so that the arrows are treated like
    # their own color
    colors = _color_list(colors_w_arrows)
    # find the total number of atoms
    n = sum(concs)

    # Find the mixed radix number counter for the system
    C = _coefficients(concs,n)

    # count how many arrows there are
    (narrows,arrow_types,concs_w_arrows) = how_many_arrows(colors_w_arrows)
    # prepare the stabalizer array to be the appropriate size
    stabalizer = [0]*(len(concs)+1)
    
    # ast is the stabalizer for the arrow permutaiton order is used to
    # track which branchs of the first layer we have already seen so
    # we can skip over the ones we already know aren't unique
    ast = []
    order = {}
    for i in range(0,C[0]+1):
        order[i] = i

    # determine if we are finding a subset of doing a full enumeration
    if subset is not None and isinstance(subset, list) and len(subset) > 0:
        use_subset = True
    elif (subset is None or isinstance(subset, int)) and accept is not None:
        use_subset = True
    else:
        use_subset = False
				
    # count is an integer counter used to keep track of where in
    # the total number of configurations we are.
    count = 1

    #create the all zero branch of lables
    #survivors is an array that stores the unique configurations
    from numpy import zeros
    branch = zeros(len(C))
    survivors = []

    #variables used for iteration and as test criteria for
    #navigating the tree
    # i keeps track of which level of the tree we are on
    # b0 lets us know when we've used all the needed colors/been
    # to all the relevant levels
    i = 0
    b0 = 0

    from random import random
    from phenum.msg import verbosity
    if verbosity is not None and verbosity >= 1: #pragma: no cover
        from tqdm import tqdm
        if isinstance(subset, list) and len(subset) > 0:
            ntotal = len(subset)
        elif isinstance(subset, int):
            ntotal = subset
        else:
            ntotal = total
        pbar = tqdm(total=ntotal)
    # Now we loop through the different possible hash arrays until
    # they have all been considered
    ncurrent = 0
    while branch[0] < C[0] and ((subset is not None and isinstance(subset, list) and len(survivors) < len(subset))
                                or use_subset == False
                                or (accept is not None and isinstance(subset, int) and len(survivors) == subset)):
	# perm determines if the new array is unique, if yes unique =
	# 0, if no then unique = 1, perm also outputs the stabilizers
	# for each level, the order for the first level, and the
	# stabilizers for the arrow configurations
        (unique,stabalizer[i+1],order,ast) = _perm(branch,concs,n,i,group,stabalizer[i],order,ast)

        if not supers:
            from operator import mul
            # if we have a unique structure then we need to check if
            # it's super periodic before saving it.
            if unique == 0 and (b0 == 1 or i == 0):
                brancht = list(_invhash(branch, concs, len(colors_w_arrows)))
                for trans in range(1,int(cellsize)):
                    action = group[trans][0]
                    trans_branch = []
                    for act in action:
                        trans_branch.append(brancht[act])
                    if trans_branch == brancht:
                        unique = 1
                        if use_subset:
                            if count in subset:
                                loc = subset.index(count)
                                nv = count
                                while nv in subset and nv < reduce(mul,C,1):
                                    nv += 1
                                subset[loc] = nv

                            count += 1
                        break
        
	# if this array is unique then we may need to append it to the
	# list of survivors
        if unique == 0:
	    # if we've visited every level of the tree
	    # then we need to append the unique
	    # configurations
            if b0 == 1 or i == 0:
		# if there are arrows we need to see
		# what their unique arrangements are
		# before saving the configurations
                if narrows > 0:
		    # first use invhash to turn the hash back to an array
                    brancht = list(_invhash(branch, concs, len(colors_w_arrows)))
		    # create a temporary copy of the branch.
                    tbrancht = []
                    for tt in range(len(brancht)):
                        tbrancht.append(colors[brancht[tt]-1][1])
		    # make a coloring with arrows to be passed to the
		    # arrow permutiation code by adding arrows back
		    # into the array where needed.
                    for z in range(len(brancht)):
                        brancht[z] = deepcopy(colors[brancht[z] -1])
		    # add_arrows from the phonon_brancher code returns
		    # the unique configurations with the unique arrow
		    # arrangements
                    arsurvivors = add_arrows(brancht,ast, dim, accept, True)
		    #write the unique confgurations to file.
                    for z in arsurvivors:
			# if we aren't using a subset write everything
			# to file.
                        if not use_subset:
                            survivors.append(z)
			# if we are then we only want to write the
			# random subset to file.
                        else:
                            if subset is not None and isinstance(subset, list) and count in subset:
                                survivors.append(z)
                            elif accept is not None:
                                survivors.append(z)
                        count += 1
                else:
		    # if there are no arrows just
		    # write the unique arrangement
		    # to file
                    if not use_subset:
                        tbranch = list(_invhash(branch, concs, len(colors_w_arrows)))
                        tbranch = [colors[leaf -1][1] for leaf in tbranch]
                        survivors.append([[-1,leaf] for leaf in tbranch])
                    else:
                        if subset is not None and isinstance(subset, list) and count in subset:
                            tbranch = list(_invhash(branch, concs, len(colors_w_arrows)))
                            tbranch = [colors[leaf -1][1] for leaf in tbranch]
                            survivors.append([[-1,leaf] for leaf in tbranch])
                        elif accept is not None and random() < accept:
                            tbranch = list(_invhash(branch, concs, len(colors_w_arrows)))
                            tbranch = [colors[leaf -1][1] for leaf in tbranch]
                            survivors.append([[-1,leaf] for leaf in tbranch])
                    count += 1

		    # if we aren't on the last contributing level
		    # of the tree then make sure b0 is 0
            if i < len(branch) - 2:
                b0 = 0

	# all that follows dictates how we navigate the tree. This can
	# be a little tricky but does work at this time. Basically try
	# not to mess with the last few lines of this code if it can
	# be helped
								
	# if b0 is 0 then we need to move up to the next level of the
	# tree at this point
        if b0 == 0:
            if i < len(branch)-1:
                i += 1
	    # if we are now looking at the last level then set b0 to 1
	    # so that any unique configurations found will be saved
            if i >= len(branch) - 2:
                b0 = 1

	# root is a variable that tells us if we've filled this
	# level of the tree yet
        root = 0
	# if root = 0 and i>0 then we need to check to see if
	# we need to back track up the tree yet
        while root == 0 and i > 0:
	    # if we've used all the configurations for this level then
	    # we need to go back a level, if not then set root to 1 so
	    # that we keep moving forward through the tree
            if branch[i] == C[i]-1:
		# first reset this level's number so that it will
		# increment properly the next time it's visited
                branch[i] = 0  
		# next reduce i by 1 to go back a level
                i -= 1
		# if i = 0 we're back to the bottom of the tree. We
		# need to reset the stabilizers, b0 and increment the
		# bottom layers index, branch[i], by 1
                if i == 0:
                    stabalizer = [0]*len(concs)
                    b0 = 0
                    branch[i] += 1
            else:
                root = 1
                
	# if b0=1 then we need to increase this level index,
	# branch[i], by 1 to consider the next configuration
        if b0 == 1:
            branch[i] += 1

	# test helps us check the first level for areas we can skip
        test = 0
        while test == 0:
	    # if we're on the first branch and the order
	    # for this configuration is negative then we
	    # skip it. Otherwise we need to consider it.
            if i == 0 and order[branch[i]] == -1:
                branch[i] += 1
            else:
                test = 1
                
        if len(survivors) > ncurrent:
            ncurrent = len(survivors)
            if verbosity is not None and verbosity >= 1: #pragma: no cover
                pbar.update(1)
                
    if verbosity is not None and verbosity >= 1: #pragma: no cover
        pbar.close()
    # done
    return survivors

