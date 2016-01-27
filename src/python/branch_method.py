#binomial_unique_permutations.py
#import needed modules
import numpy as np
import radix_num_generator as rng
import Inverse_radix_num as irn 
import binomial_calculator as bc
import phonon_brancher as pb
from copy import deepcopy
import random
#this method applies a cyclical permutation to an array and removes
#duplicates from a list if found It calls for casi, the idnumber for
#the case to be permuted, survivors, a list of id numbers to be
#compared to and added to, Coefficient as an array of integers
#represendting mixed radix numbers for the systm and colors is an
#array of integers of equal length that represents the number of each
#color found in our array
#ast is the stabalizer for the final configuration that will be passed
#to the arrow permutation routine.
def perm(casei,colors,length,index,gen,stab,order,ast):
        #sets a test value, finds the length of the array using colors
	unique = 0
        if stab == 0:
                stab = []
        st = []
        ast = []
        #uses the invid method from Inverse_radix_num to turn the id
        #number, casei, into an array
	li = list(irn.invhash(casei,colors,length))
	if index == 0:
                for a in gen:
                        i = a[0]
                        lnew = []
                        for j in i:
                                lnew.append(li[j])
                                #uses idnum and hash from
                                #radix_num_generator to turn out array
                                #back into a number
                        lnewb = list(lnew)
                        for y in range(0,len(lnew)):
                                if lnew[y] != index + 1:
                                        lnewb[y] = 0
                                elif lnew[y] == index + 1:
                                        lnewb[y] = 1
			rtest = rng.hash(list(lnewb),colors)
                        if casei[0] < rtest:
                                order[rtest] = -1
			elif casei[0] > rtest:
                                #if the shifted array's number is
                                #found in the list of survivors then
                                #set unique to 1 of false
				unique += 1
                                st = []
                                ast = []
				break
                        elif casei[0] == rtest:
                                st.append(a)
                        if li == lnew:
                                ast.append(a)
        elif len(stab) == 1:
                st = stab
                ast = stab
        elif len(stab) > 1:
                ast = []
                for a in stab:
                        i = a[0]
                        lnew = []
                        for j in i:
                                lnew.append(li[j])
                        lib= list(li)
                        lnewb = list(lnew)
                        lnewb2 = []
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
                        if lnewb == lib:
                                st.append(a)
                        if li == lnew:
                                ast.append(a)

                        rtest = rng.hash(list(lnewb2),colors)
                        if casei[index] > rtest:
                                unique = 1
                                st = []
                                ast = []
                                break
                        if unique == 1:
                                ast = []
                                st = []
                                break
        else:
                unique = 1
	return(unique,st,order,ast)

# uses the Coefficients method from rdix_num_generator to find the
# radix numbers for the system
# dim is the number of directions the arrows can point.
# concs is the concentrations of the atoms
# group is the set of group operations
# colors_w_arrows is the colors popelating the lattice with theri
# displacement directions indicated
# wanted is an optional integer that indicates of a subset of the
# total number is wanted.
# total is an optional integer used to generate the desired subset
def brancher(concs,group,colors_w_arrows, dim, wanted = 0, total = 0):
        # redifine the colors so that the arrows are treated like
        # their own color
        colors = pb.color_list(colors_w_arrows)
        # find the total number of atoms
        n = sum(concs)
        # Find the mixed radix number counter for the system
        C = rng.Coefficients(concs,n)
        # count how many arrows there are
        narrows = pb.how_many_arrows(colors_w_arrows)
        # prepare the stabalizer array
        stabalizer = [0]*len(concs)
        # variables 
        # cyclicgen = []
        # coluse = 1
        # groups is the space group
        groups = []
        # ast is the stabalizer for the arrow permutaiton
        ast = []
        # strip the arrows permutation from the total group action to
        # get the space group action.
        for i in group:
                groups.append(i[0])
        order = {}
        # the order is used to determine which branch layer we are on
        # so that we can move through the tree as we nee to.
        for i in range(0,bc.binomial_coefficient(sum(concs),concs[0])+1):
                order[i] = i

        # use_subset tells the code if it's only writting out a subset
        use_subset = False
        # generate the random subset to be used
        if wanted != 0:
                use_subset = True
                subset = []
                if total == 0:
                        print("ERROR! We cannot find a subset of nothing.")
                        print("Please enter the total number of unique configurations.")
                        exit()
                while len(subset) < wanted:
                        random_num = random.randint(1,total)
                        if random_num  not in subset:
                                subset.append(random_num)

        # count is an integer counter used to keep track of where in
        # the total number of configurations we are.
        count = 1
        print('s',subset)
        #create the all zero branch of lables and use idnum from
        #radix_num_generator to turnit into a number and store it in
        #survivors
        branch = np.zeros(len(C))
        survivors = []
        i = 0
        b0 = 0
        #for the first color look at each value of it while the rest
        #of the branches are zero
        while branch[0] < C[0]:
        # sum(branch) != sum(C) - len(C) and branch[0] < C[0]-1: uses
        #perm to determine if the new array is unique, if yes unique =
        #0, if no then unique = 1
                (unique,stabalizer[i+1],order,ast) = perm(branch,concs,n,i,group,stabalizer[i],order,ast)
                #if this array is unique append it to the list of
                #survivors
                if unique == 0:
                        if b0 == 1 or i == 0:                                
                                if narrows > 0:
                                        brancht = list(irn.invhash(branch, concs, len(colors_w_arrows)))
                                        #stuff to be saved to file.
                                        f = open('phonon_out.txt', 'a')
                                        ttbranch = []
                                        for tt in range(len(branch)-1):
                                                ttbranch.append(int(branch[tt]))
                                        if narrows[0] == len(brancht):
                                                f.write('[0, '+str(ttbranch)[1:]+', ')
                                        else:
                                                f.write(str(ttbranch)+', ')
                                        tbrancht = []
                                        for tt in range(len(brancht)):
                                                tbrancht.append(colors[brancht[tt]-1][1])
                                        f.write(str(tbrancht)+'\n')
                                        #make a coloring with arrows
                                        #to be passed to the arrow
                                        #permutiation code.
                                        for z in range(len(brancht)):
                                                brancht[z] = deepcopy(colors[brancht[z] -1])
                                        arsurvivors = pb.add_arrows(brancht,ast, dim)
                                        #write the unique confgurations to file.
                                        for z in arsurvivors:
                                                # if we aren't using a
                                                # subset write
                                                # everything to file.
                                                if not use_subset:
                                                        survivors.append(z)
                                                        f.write(str(z) + '\n')
                                                # if we are then we
                                                # only want to write
                                                # the random subset to
                                                # file.
                                                else:
                                                        if count in subset:
                                                                survivors.append(z)
                                                                f.write(str(z) + '\n')
                                                count += 1
                                else:
                                        # if there are no arrows just
                                        # write the unique arrangement
                                        # to file
                                        if not use_subset:
                                                survivors.append(list(irn.invhash(branch, concs, len(colors_w_arrows))))
                                        else:
                                                if count in subset:
                                                        survivors.append(list(irn.invhash(branch, concs, len(colors_w_arrows))))
                                        count += 1
                                        
                        if i < len(branch) - 2:
                                b0 = 0
                if b0 == 0:
                        i += 1
                        if i == len(branch) - 2:
                                b0 = 1
                root = 0
                while root == 0 and i > 0:
                        if branch[i] == C[i]-1:
                                branch[i] = 0
                                i -= 1
                                if i == 0:
                                        stabalizer = [0]*len(concs)
                                        b0 = 0
                                        branch[i] += 1
                        else:
                                root = 1
                if b0 == 1:
                        branch[i] += 1
                test = 0
                while test == 0:
                        if i == 0 and order[branch[i]] == -1:
                                branch[i] += 1
                        else:
                                test = 1
        return(survivors)   
