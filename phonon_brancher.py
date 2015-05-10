#This program is meant to find the possible unique configurations of a
#set of colors with optional displacement vectors attached to each.

#import neede modules
import numpy as np
import branch_method as bm
import Inverse_radix_num as irn 
from copy import deepcopy
#Find concentration takes a 1D array, col, and returns a 1D array, 
#Concs, containing the number of times each element in col appears 
#in the array.
def find_concentrations(col):
    Concs=[]

    tcol=col[:]
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
def how_many_arrows(col):
    arrows=0
    
    tcol=col[:]
    for i in range(len(tcol)):
        if len(tcol[i]) == 2:
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

#advance_arrows takes one of the arrows as input and advances it to
#the next possible arrow.
def advance_arrows(arrow):
    if arrow in(range(1,5)):
        arrow = arrow%4 + 1
    else:
        print("Error, arrow not in expected range. Please ensure that the arrows range between 1 and 4.")
    return(arrow)
#initialize the neede inputs generators is a 2D array that contains
#the generators for the group to be used.  col is a 1D array thot
#contians the colors to be used in the enumeration: a singe letter
#means the atom is not being displaced, a letter and a number means
#that the atom is being displaced.
generators=[[0,1,2,3],[1,2,3,0]]
col=[[1,1],[1],[1],[2,1]]

#first we need to sort the colors according to size so that we
#can enumerate the colors that aren't being displaced first.
col=sorted(col, key=len)

#We also need to know the concentrations of the atoms as they 
#occure in the sorted list. This will be stored in concs.
Concs = find_concentrations(col)

#Also find the colors being used for this system. If a color has an
#arrow tag then it is it's own color distinct from those without
#colors.
colors = color_list(col)

#Find out how many arrows there are in the col array
narrows = how_many_arrows(col)
#generate the list of symmentrically distinct configurations for the
#different colors being used. Each if a color has on arrow on it then
#it is treated as it's own color.
#group contains the full group from the generators.
#configs contains the list of unique configurations.
#t contains the time it took to generate the configurations.
(group,configs,t) = bm.brancher(Concs,generators)

#Convert the configurations from the form used in the branch_methods
#routine to the one we are using locally.
for i in range(len(configs)):
    configs[i] = list(irn.invhash(configs[i],Concs,len(col)))

for i in range(len(configs)):
    for j in range(len(configs[i])):
        configs[i][j] = colors[configs[i][j]-1]

#survivors is the array that contains the end result of the permutations
survivor = []

#This array tells us when we've applied all the possible arrow
#directions to the colors that will take arrows
for col in configs:
    function = np.zeros(len(col))
    for i in range(len(col)):
        if len(col[i]) == 1:
            function[i] = 0
        elif len(col[i]) == 2:
            function[i] = 4

    #Set the branch to match the configuration of the colors that we are
    #handed.
    branch=[]
    for i in col:
        branch.append(i)

    #If there are no arrows on the lattice then append it to the survivors
    #list and be done.
    i=0
    if narrows == 0:
        survivor.append(branch[:])
    else:
        visited = 0
        #while there are arrows that still need to be rotated keep doing
        #this. When all the elements of funtcion are zero we are done.
        while (any(function !=0) or visited != 1):
            visited = 0
            if len(col[i]) == 1:
                branch[i]=branch[i]
            elif len(col[i])==2:
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
                        if function[i] ==0 and any(function != 0):
                            if len(col[i]) == 2:
                                function[i] =4
                            i -= 1
                        else:
                            mv = 1

#done
print('done')

for i in survivor:
    print(i)
print(len(survivor))            

