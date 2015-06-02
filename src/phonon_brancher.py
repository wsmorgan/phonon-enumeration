#This module contains side programs needed for the phonon_enumeration
#code.

#import neede modules
import numpy as np
from copy import deepcopy
import radix_num_generator as rn
#Find concentration takes a 1D array, col, and returns a 1D array, 
#Concs, containing the number of times each element in col appears 
#in the array.
def find_concentrations(col):
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
    arrows=0

    for i in tcol:
        if i[1] >= 0:
            arrows += 1
    return(arrows)

#color_list takes the configuration and returns a unique list of the
#colors used without duplicating any duplicates.
def color_list(col):
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

#An algorithm for sorting the col array so that colors without arraws
#appear first in the list by concentration followed by the colors with
#arraws that are also sorted by concentration.
def col_sort(col_list):

    col1 = []
    col2 = []
    colt=[]
    for i in col_list:
        if i[1] < 0:
            col1.append(i)
        elif i[1] >= 0:
            col2.append(i)
    col1 = sorted(col1, key = col1.count, reverse=True)
    col2 = sorted(col2, key = col2.count, reverse=True)
    for i in col1:
        colt.append(i)
    for i in col2:
        colt.append(i)
    
    return(colt)

#advance_arrows takes one of the arrows as input and advances it to
#the next possible arrow.
def advance_arrows(arrow):
    x = [1,2,3,0]
    if arrow in(range(0,4)):
        arrow = x[arrow]
    else:
        print("Error, arrow not in expected range. Please ensure that the arrows range between 0 and 3.")
    return(arrow)

#this method takes a configuration and finds the unique ways of
#placing the arrows on that configuration.
#col is the initial configuration.
#agroup is the group operations with their effects on the arrows.
def add_arrows(col,agroup,dim):
    #survivors is the array that contains the end result of the permutations
    survivor = []

    #Find out how many arrows there are in the col array
    narrows = how_many_arrows(col)

    # largest_arrow = [dim-1]*narrows
    
    #This array tells us when we've applied all the possible arrow
    #directions to the colors that will take arrows
    function = np.zeros(len(col))
    for i in range(len(col)):
        if col[i][1] < 0:
            function[i] = -1
        else:
            function[i] = dim

    #Set the branch to match the configuration of the colors that we are
    #handed.
    branch=[]
    for i in col:
        branch.append(i)

    #If there are no arrows on the lattice then append it to the survivors
    #list and be done.
    print('narrows', narrows)
    i=0
    if narrows == 0:
        survivor.append(branch[:])
        visited = 1
    else:
        visited = 0
    #while there are arrows that still need to be rotated keep doing
    #this. When all the elements of funtcion are zero we are done.
    while (any(function > 0) or visited != 1):
        visited = 0
        if col[i][1] < 0:
            branch[i]=branch[i]
        elif col[i][1] >= 0:
            branch[i][1] = advance_arrows(branch[i][1])
            function[i] -= 1
        else:
            print("ERROR: cannot have an entry with more than two components.")
            
        #If we reached the last element in our configuration save it
        #to the survivors.
        if i == len(col)-1:
            survivor.append(deepcopy(branch))
            visited = 1
        #if we aren't on the last element of the configuration then
        #move to the next element
        if i != len(col) -1:
            i += 1
        #if we are now on an element that takes an arrow check that it
        #still has more rotations of the arrow that need to be
        #applied. If it does then don't switch elements. If it doesn't
        #then move to a previouse element that still has more arrows
        #that need to be applied and reset the funtion to it's maximum
        #value.
        if i >= len(col)-1 and visited==1:
            if any(function !=0):
                mv = 0
                while mv ==0 and i >= 0:
                    if function[i] <= 0 and any(function > 0):
                        if function[i] == 0:
                            function[i] =4
                        i -= 1
                    else:
                        mv = 1

    #Now that we have all the possible ways of putting the arrows on
    #the lattice we need to check to see which of them are unique.
    #arsurvivors is the final list of unique configurations
    
    #tempsurv is a list of of the effects of the group operations on
    #the configurations
    arsurvivors = []
    tempsurv = []
    print(survivor)
    for x in survivor:
        salist = rn.ahash(x,dim)
        print(salist)
        unique = True
        #apply each of the group operations to the lattice.
        for i in agroup:
            lnew=[]
            arrow_perm = i[1]
            for j in i[0]:
                nar = []
                nar.append(x[j][0])
                #if there is an arrow on any lattice site it needs to
                #be updated. Otherwise just permute the colors.
                if x[j][1] >= 0:
                    nar.append(arrow_perm[x[j][1]])
                else:
                    nar.append(x[j][1])
                lnew.append(nar)
            #if the new configuration is not already in the list of
            #permutations and is not the initial configuration add it
            #to tmpsurv
            calist = rn.ahash(lnew,dim)
            print(calist)
            if calist < salist:
                unique = False
                break
        #if the original configuration isn't in the tempsurv list or
        #the list of arrow survivors add it to the list.
        if unique == True:
            arsurvivors.append(x)
    # arsurvivors = []
    # tempsurv = []
    # for s in range(rn.ahash(largest_arrow,dim)):
    #     salist = s # rn.ahash(x)
    #     x = rn.ainvhash(s,narrows,dim)
    #     unique = True
    #     #apply each of the group operations to the lattice.
    #     for i in agroup:
    #         print(i)
    #         print(i[1])
    #         lnew=[]
    #         arrow_perm = i[1]
    #         print(arrow_perm)
    #         print(x)
    #         for j in range(len(x)):
    #             print(j)
    #             lnew.append(arrow_perm[x[j]])
    #             # nar = []
    #             # nar.append(x[j][0])
    #             #if there is an arrow on any lattice site it needs to
    #             #be updated. Otherwise just permute the colors.
    #             # if x[j][1] >= 0:
    #             #     nar.append(arrow_perm[x[j][1]])
    #             # else:
    #             #     nar.append(x[j][1])
    #             # lnew.append(nar)
    #         #if the new configuration is not already in the list of
    #         #permutations and is not the initial configuration add it
    #         #to tmpsurv
    #         calsit = rn.ahash(lnew,dim)
    #         if calsit < salist:
    #             unique = False
    #             break
    #     #if the original configuration isn't in the tempsurv list or
    #     #the list of arrow survivors add it to the list.
    #     if unique == True:
    #         nbranch = [0]*len(col)
    #         i = 0
    #         for z in range(len(col)):
    #             if col[z][1] < 0:
    #                 nbranch[z] = deepcopy(col[z])
    #             elif col[z][1] >= 0:
    #                 nbranch[z] = deepcopy(col[z])
    #                 nbranch[z][1] = x[i]
    #                 i += 1
                    
    #         arsurvivors.append(nbranch)
            
    return(arsurvivors)


