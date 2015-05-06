#This program is meant to find the possible unique configurations of a
#set of colors with optional displacement vectors attached to each.

#import neede modules
import numpy as np
import branch_method as bm
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

def how_many_arrows(col):
    arrows=0
    
    tcol=col[:]
    for i in range(len(tcol)):
        if len(tcol[i]) == 2:
            arrows += 1
    return(arrows)

#initialize the neede inputs generators is a 2D array that contains
#the generators for the group to be used.  col is a 1D array thot
#contians the colors to be used in the enumeration: a singe letter
#means the atom is not being displaced, a letter and a number means
#that the atom is being displaced.
generators=[[0,1,2,3],[1,2,3,0]]
col=[[1,1],[1],[1],[2,1],[3,1],[5,1]]

#first we need to sort the colors according to size so that we
#can enumerate the colors that aren't being displaced first.
col=sorted(col, key=len)

#We also need to know the concentrations of the atoms as they 
#occure in the sorted list. This will be stored in concs.
Concs = find_concentrations(col)

#Find out how many arrows there are in the col array
narrows = how_many_arrows(col)

#survivors is the array that contains the end result of the permutations
survivors = []

#This array tells us when we've applied all the possible arrow
#directions to the colors that will take arrows
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
    branch.append(i[0])

#If there are no arrows on the lattice then append it to the survivors
#list and be done.
i=0
if narrows == 0:
    survivors.append(branch)
else:
    #while there are arrows that still need to be rotated keep doing
    #this. When all the elements of funtcion are zero we are done.
    while any(function !=0):
        if len(col[i]) == 1:
            branch[i]=branch[i]
        elif len(col[i])==2:
            branch[i] = branch[i]*1j
            function[i] -= 1
        else:
            print("ERROR: cannot have an entry with more than two components.")
        #If we reached the last element in our configuration save it
        #to the survivors.
        if i == len(col)-1:
            survivors.append(branch)
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
        elif i == len(col) -1 and len(col[i]) == 2:
            mv = 0
            while mv ==0 and i > 0:
                if function[i] ==0 and any(function != 0):
                    function[i] =4
                    i -= 1
                    if len(col[i]) == 1:
                        i += 1
                        mv =1
                else:
                    mv = 1

for i in survivors:
    print(i)
print(len(survivors))            

