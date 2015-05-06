#This program is meant to find the possible unique configurations of a
#set of colors with optional displacement vectors attached to each.

#import neede modules
import numpy as np
import generator as gn

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
col=[[1],[1,1],[1],[2,1]]

#first we need to sort the colors according to size so that we
#can enumerate the colors that aren't being displaced first.
col=sorted(col, key=len)

#We also need to know the concentrations of the atoms as they 
#occure in the sorted list. This will be stored in concs.
Concs = find_concentrations(col)

#Find out how many arrows there are in the col array
narrows = how_many_arrows(col)

#Now we need to populate the lattice with the atoms we want one at a
#time and store the lattices in the survivors list.  
#branch is an array that represents the lattice and will be populated
#with the colors from the col array.
survivors = []
branch = col[:]
nbranch = []

print(branch)
print(nbranch)

for x in range(len(branch)):
    if len(branch[x]) ==1:
        nbranch.append(branch[x])
    if len(branch[x]) ==2:
        for i in range(4):
            if i==0:
                nbranch.append(branch[x])
            elif i==1:
                print()
            elif i==2:
                print()
            elif i==3:
                print()
            for j in range(narrows):
                print('h')
        
print(nbranch)
