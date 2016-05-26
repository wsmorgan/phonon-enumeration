"""The heart of the code is here. These methods build the tree and
check the canfigurations for uniqueness."""

#import needed modules
import numpy as np
import radix as rng
import phonons as pb
from copy import deepcopy
import random

def perm(casei,colors,length,index,gen,stab,order,ast):
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
    li = list(rng.invhash(casei,colors,length))

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
	    rtest = rng.hash(list(lnewb),colors)

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
	    rtest = rng.hash(list(lnewb2),colors)
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


def brancher(concs,group,colors_w_arrows, dim, subset = []):
    """This routine navigates the tree and saves the unique configurations
    to an array survivors.
		
    :arg concs: the concentrations of the colors or atoms
    :arg group: the symmetry group for the system, including the
    effect on the arrows (displacement directions)
    :arg colors_w_arrows: an integer 2D array that indiciates
    which atoms are being displaced, i.e., where the arrows are.
    :arg dim: the number of directions the arrows can point
    :arg subset: an integer array of the subset of unique
    arrangements wanted
	
    The method returns the list of unique configurations and the
    number of stabilizers for the last level of the tree.

    survivors: is the 3D array of the unique arrangements of
    colors and arrows
    stob_len: the number of stabilizers for a given unique
    configuration
    """

    # initial setup	

    # redifine the colors so that the arrows are treated like
    # their own color
    colors = pb.color_list(colors_w_arrows)
    # find the total number of atoms
    n = sum(concs)

    # Find the mixed radix number counter for the system
    C = rng.coefficients(concs,n)

    # count how many arrows there are
    narrows = pb.how_many_arrows(colors_w_arrows)
    # prepare the stabalizer array to be the appropriate size
    stabalizer = [0]*len(concs)
    
    # ast is the stabalizer for the arrow permutaiton order is used to
    # track which branchs of the first layer we have already seen so
    # we can skip over the ones we already know aren't unique
    ast = []
    order = {}
    for i in range(0,rng.binomial_coefficient(sum(concs),concs[0])+1):
	order[i] = i

    # determine if we are finding a subset of doing a full enumeration
    if len(subset) > 0:
	use_subset = True
    else:
	use_subset = False
				
    # count is an integer counter used to keep track of where in
    # the total number of configurations we are.
    count = 1

    #create the all zero branch of lables
    #survivors is an array that stores the unique configurations
    branch = np.zeros(len(C))
    survivors = []

    # So that we can eliminate the superperiodic structures we
    # need to know the number of stabilizers for each survivor so
    # we'll create a list of the number of stabilizers as well.
    stab_len = []

    #variables used for iteration and as test criteria for
    #navigating the tree
    # i keeps track of which level of the tree we are on
    # b0 lets us know when we've used all the needed colors/been
    # to all the relevant levels
    i = 0
    b0 = 0
		
    # Now we loop through the different possible hash arrays until
    # they have all been considered
    while branch[0] < C[0] and (len(survivors) < len(subset) or use_subset == False):
	# perm determines if the new array is unique, if yes unique =
	# 0, if no then unique = 1, perm also outputs the stabilizers
	# for each level, the order for the first level, and the
	# stabilizers for the arrow configurations
	(unique,stabalizer[i+1],order,ast) = perm(branch,concs,n,i,group,stabalizer[i],order,ast)

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
		    brancht = list(rng.invhash(branch, concs, len(colors_w_arrows)))
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
		    arsurvivors = pb.add_arrows(brancht,ast, dim)
		    #write the unique confgurations to file.
		    for z in arsurvivors:
			# if we aren't using a subset write everything
			# to file.
			if not use_subset:
			    survivors.append(z)
			    stab_len.append(len(ast))
			# if we are then we only want to write the
			# random subset to file.
			else:
			    if count in subset:
				survivors.append(z)
				stab_len.append(len(ast))
			count += 1
		else:
		    # if there are no arrows just
		    # write the unique arrangement
		    # to file
		    if not use_subset:
			survivors.append(list(rng.invhash(branch, concs, len(colors_w_arrows))))
			stab_len.append(len(ast))
		    else:
			if count in subset:
			    survivors.append(list(rng.invhash(branch, concs, len(colors_w_arrows))))
			    stab_len.append(len(ast))
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
	    i += 1
	    # if we are now looking at the last level then set b0 to 1
	    # so that any unique configurations found will be saved
	    if i == len(branch) - 2:
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

    # done
    return(survivors,stab_len)
