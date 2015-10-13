#binomial_unique_permutations.py
#import needed modules
import numpy as np
import radix_num_generator as rng
import Inverse_radix_num as irn 
import binomial_calculator as bc
import phonon_brancher as pb
from copy import deepcopy
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

#uses the Coefficients method from rdix_num_generator to find the
#radix numbers for the system
#dim is the number of directions the arrows can point.
def brancher(concs,group,colors_w_arrows, dim):
        colors = pb.color_list(colors_w_arrows)
        n = sum(concs)
        C = rng.Coefficients(concs,n)
        narrows = pb.how_many_arrows(colors_w_arrows)
        stabalizer = [0]*len(concs)
        cyclicgen = []
        coluse = 1
        groups = []
        ast = []
        for i in group:
                groups.append(i[0])
        order = {}
        for i in range(0,bc.binomial_coefficient(sum(concs),concs[0])+1):
                order[i] = i

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
                                        # f.write(str(ttbranch)+', ')
                                        # f.write(str(brancht)+'\n')
                                        #make a coloring with arrows
                                        #to be passed to the arrow
                                        #permutiation code.
                                        for z in range(len(brancht)):
                                                brancht[z] = deepcopy(colors[brancht[z] -1])
                                        arsurvivors = pb.add_arrows(brancht,ast, dim)
                                        #write the unique confgurations to file.
                                        for z in arsurvivors:
                                                survivors.append(z)
                                                f.write(str(z) + '\n')
                                else:
                                        survivors.append(list(irn.invhash(branch, concs, len(colors_w_arrows))))
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
