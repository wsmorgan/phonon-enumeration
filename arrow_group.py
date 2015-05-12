#arrow_group.py

#a_group take the set of translations of the lattice and the set of
#rotations of the lattice paired with their effect on the arrows of
#the lattice. It makes a direct product of the transformations and the
#rotations and pairs them with the combination of th arraws. If the
#resultant operation is not in the group it then adds it to the group.
def a_group(trans,rots):
    q = 1
    groupi = []
    for i in trans:
        for j in rots:
            c = []            
            c.append([j[0][i[l]] for l in range(0,len(i))])
            c.append(j[1])
            if c not in groupi:
                groupi.append(c)
        
    return(groupi)

