#phonon_enumeration.py
#import needed modules
import branch_method as bm
import phonon_brancher as pb
import arrow_group as ar


print("Performing initial setup:")

#initialize the needed inputs generators is a 2D array that contains
#the generators for the group to be used.  col is a 1D array thot
#contians the colors to be used in the enumeration: a number combined
#with a negative number, i.e. [1,-1], means the atom is not being
#displaced, a number combined with a positive number between 0 and 3,
#i.e. [1,2] means that the atom is being displaced.
#trans is the quotient group, or translations of the lattice
#rots is the group of rotations on the lattice paired with their
#effects on the arrows.
col=[[1,-1],[1,-1],[1,-1],[2,1],[1,-1],[1,-1],[1,-1],[1,-1],[1,-1]]
trans=[[0,1,2,3,4,5,6,7,8],[1,2,0,4,5,3,7,8,6],[2,0,1,5,3,4,8,6,7],[3,4,5,6,7,8,0,1,2],[4,5,3,7,8,6,1,2,0],[5,3,4,8,6,7,2,0,1],[6,7,8,0,1,2,3,4,5],[7,8,6,1,2,0,4,5,3],[8,6,7,2,0,1,5,3,4]]
rots=[[[0,1,2,3,4,5,6,7,8],[0,1,2,3]],[[2,5,8,1,4,7,0,3,6],[1,2,3,0]],[[8,7,6,5,4,3,2,1,0],[2,3,0,1]],[[6,3,0,7,4,1,8,5,2],[3,0,1,2]],[[8,5,2,7,4,1,6,3,0],[1,0,3,2]],[[0,3,6,1,4,7,2,5,8],[3,2,1,0]],[[2,1,0,5,4,3,8,7,6],[0,3,2,1]],[[6,7,8,3,4,5,0,1,2],[2,1,0,3]]]

# trans = [[0,1,2,3,4,5],[1,2,0,5,4,3],[2,0,1,4,5,3],[5,4,3,2,1,0],[4,3,5,0,2,1],[3,5,4,1,0,2]]
# rots = [[[0,1,2,3,4,5],[0,1,2,3]],[[3,4,5,0,1,2],[2,3,0,1]],[[5,4,3,2,1,0],[2,1,0,3]],[[2,1,0,5,4,3],[0,3,2,1]]]
# col = [[3,-1],[1,-1],[2,1],[2,1],[1,-1],[1,-1]]                                  

#we start by sorting the colors to be used so that arraws appear last
#and the highest concentration of the colors is handled first.
col = pb.col_sort(col)

print('Now combining the effect of the arrows with the group action.')
#Find the group using the translations and the rotations of the
#lattice.
agroup = ar.a_group(trans,rots)

print('Done!')

#We also need to know the concentrations of the atoms as they 
#occure in the sorted list. This will be stored in concs.
Concs = pb.find_concentrations(col)

# generate the list of symmentrically distinct configurations for the
# different colors being used. Each if a color has on arrow on it then
# it is treated as it's own color.
#Concs is the list of concentrations of each color.
# group contains the full group from the generators.
#col is the sorted list of the colors.
# configs contains the list of unique configurations.
# t contains the time it took to generate the configurations.
print("Finding unique configurations of colors and arraws.")
(configs,t) = bm.brancher(Concs,agroup,col)

print("Found all unique configurations")
print("They are:")

#print the results.
for i in configs:
    print(i)
print(len(configs))

print("And it took", t, "seconds.")
