#binomial_unique_permutations.py
#import needed modules
import numpy as np
import radix_num_generator2 as rng
import Inverse_radix_num as irn 
import time
import generator as gn
import binomial_calculator as bc

#this method applies a cyclical permutation to an array and removes duplicates from a list if found
#It calls for casi, the idnumber for the case to be permuted, survivors, a list of id numbers to
#be compared to and added to, Coefficient as an array of integers represendting mixed radix numbers
#for the systm and colors is an array of integers of equal length that represents the number of each
#color found in our array
def perm(casei,colors,length,index,gen,stab,order):
        #sets a test value, finds the length of the array using colors
	unique = 0
        if stab == 0:
                stab = []
        st = []
        #uses the invid method from Inverse_radix_num to turn the id number, casei, into an array
	li = list(irn.invhash(casei,colors,length))
	if index == 0:
                for i in gen:
                        lnew = []
                        for j in i:
                                lnew.append(li[j])
                                #uses idnum and hash from radix_num_generator to turn out array back into a number
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
                                #if the shifted array's number is found in the list of survivors then set unique to 1 of false
				unique += 1
                                st = []
				break
                        elif casei[0] == rtest:
                                st.append(i)
        elif len(stab) == 1:
                st = stab
        elif len(stab) > 1:
                for i in stab:
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
                                st.append(i)
                        rtest = rng.hash(list(lnewb2),colors)
                        if casei[index] > rtest:
                                unique = 1
                                st = []
                                break
                        if unique == 1:
                                st = []
                                break
        else:
                unique = 1
	return(unique,st,order)

#uses the Coefficients method from rdix_num_generator to find the radix numbers for the system
def brancher(col,generators):
        n = sum(col)
        C = rng.Coefficients(col,n)
        stabalizer = [0]*len(col)
        cyclicgen = []
        coluse = 1
        groups = gn.group(generators)
        order = {}
        for i in range(0,bc.binomial_coefficient(sum(col),col[0])+1):
                order[i] = i

#create the all zero branch of lables and use idnum from radix_num_generator to 
#turnit into a number and store it in survivors
        start = time.clock()
        branch = np.zeros(len(C))
        survivors = []
        i = 0
        b0 = 0
#for the first color look at each value of it while the rest of the branches are zero
        while branch[0] < C[0]:
                print(branch)
# sum(branch) != sum(C) - len(C) and branch[0] < C[0]-1:
#uses perm to determine if the new array is unique, if yes unique = 0, if no then unique = 1
                (unique,stabalizer[i+1],order) = perm(branch,col,n,i,groups,stabalizer[i],order)
#if this array is unique append it to the list of survivors
                if unique == 0:
                        if b0 == 1 or i == 0:
                                survivors.append(list(branch))
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
                                        stabalizer = [0]*len(col)
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
                print(branch)
        t = time.clock()-start
        return(groups,survivors,t)



# col = [2,1,1]
# #generators = [[1,2,3,4,0,6,7,8,9,5],[1,0,2,3,4,6,5,7,8,9]]
# # [[2,3,4,5,6,7,8,9,0,1],[2,3,0,1,4,5,6,7,8,9]]
# #generators =[[2,3,4,5,6,7,0,1],[2,3,0,1,4,5,6,7]]
# generators=[[0,1,2,3],[1,2,3,0]]
# #[[2,3,4,5,6,7,0,1],[2,3,0,1,4,5,6,7]]
# # generators= [[j -1 for j in i] for i in
# #              [[2,3,4,1,5,6,7,8,9,10,11,12,13,14,15,16],[1,2,3,4,6,7,8,5,9,10,11,12,13,14,15,16],[1,2,3,4,5,6,7,8,10,11,12,9,13,14,15,16],[1,2,3,4,5,6,7,8,9,10,11,12,14,15,16,13]]]
# # [[2,3,4,1,5,6,7,8,9,10,11,12,13,14,15,16],[1,2,3,4,6,7,8,5,9,10, 11,12,13,14,15,16],[1,2,3,4,5,6,7,8,10,11,12,9,13,14,15,16],[1,2,3, 4,5,6,7,8,9,10,11,12,14,15,16,13]]]
# # [[2,3,4,1,5,6,7,8],[2,1,3,4,5,6,7,8],[1,2,3,4,6,5,8,7],[1,2,3,4,6,7,5,8]]]
# # [[2,3,4,5,6,7,0,1],[2,3,0,1,4,5,6,7]]#,[6,7,0,1,2,3,4,5],[2,3,0,1,4,5,6,7]]

# (groups,survivors,t) = brancher(col,generators)

# f = open('string_odometer2.txt', 'w+')

# for h in survivors:
#         del h[-1]
#         for j in range(len(h)):
#                 h[j] = '{0:g}'.format(h[j])
#                 f.write(str(h[j]) + ' ')
#         f.write('\n')
# #        f.write(str(h[0]) + '\n')
# f.close()

# #print(survivors)
# print(len(groups))
# print(len(survivors))
# #print(t)
   
