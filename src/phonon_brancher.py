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
        if i[0] >= 0:
            arrows += 1
    return(arrows)

# #color_list takes the configuration and returns a unique list of the
# #colors used without duplicating any duplicates.
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
        if i[0] < 0:
            col1.append(i)
        elif i[0] >= 0:
            col2.append(i)
    col1 = sorted(col1, key = col1.count, reverse=True)
    col2 = sorted(col2, key = col2.count, reverse=True)
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
    #survivors is the array that contains the end result of the permutations
    survivor = []

    #Find out how many arrows there are in the col array
    narrows = how_many_arrows(col)
    #largest_arrow is largest value any arrow could have.
    largest_arrow = [dim-1]*narrows
    #this is an array for storing the unique arrangements.
    arsurvivors = []

    #this loop runs over every possible configuration of arrows
    #staring from the one with the smallest hash number and ending
    #with the one with the largest hash number to see if they are
    #unique. Each hash will be stored in orighash (original hash) then
    #compared to the permuted arrows hash (permhash) at the end.
    for orighash in range(rn.ahash(largest_arrow,dim)+1):
        arrow_config = rn.ainvhash(orighash,narrows,dim)
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
            permhash = rn.ahash(new_arrow_config,dim)
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


