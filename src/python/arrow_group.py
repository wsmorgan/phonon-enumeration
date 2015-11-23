#arrow_group.py

#a_group take the set of translations of the lattice and the set of
#rotations of the lattice paired with their effect on the arrows of
#the lattice. It makes a direct product of the transformations and the
#rotations and pairs them with the combination of th arraws. If the
#resultant operation is not in the group it then adds it to the group.
def a_group(trans,rots):
    """
    This subroutine combines that translations of the lattice with the
    rotatians of the lattice.
    """
    groupi = []
    for i in trans:
        for j in rots:
            c = []            
            c.append([j[0][i[l]] for l in range(0,len(i))])
            c.append(j[1])
            if c not in groupi:
                groupi.append(c)

    groupi = group(groupi)
    
    return(groupi)

def group(gen):
    """
    This subroutine takes the generators of the group then uses them
    to form the entire group.
    """
    # q is a counter to track how many group elements we've found
    q = 1
    # the final output group
    groupi = []
    # loop over the generators
    for i in gen:
        # a second set of generators used to ensure that a group
        # action doesn't act an its self
        gen2 = []      
        # if k is not the ith element of the generators add it to the
        # set gen2
        for k in gen:
            if k != i:
                gen2.append(k)
        # loop over all except the ith generator
        for j in gen2:
            # handle the site group and arrow group independently
            trans = [j[0][i[0][l]] for l in range(0,len(i[0]))]
            rots = [j[1][i[1][l]] for l in range(0,len(i[1]))]
            # c is the combination of the trans and rots, if it's not
            # in groupi then it should be added.
            c = [trans,rots]
            if c not in groupi:
                groupi.append(c)

            # for each new group element, c, we should also see if
            # the ith group element acting on it will make another
            # new group element.
            trans = [c[0][i[0][l]] for l in range(0,len(i[0]))]
            rots = [c[1][i[1][l]] for l in range(0,len(i[1]))]
            # d is the group action after the ith element has acted on
            # c, if it's not in groupi then add it.
            d = [trans, rots]
            if d not in groupi:
                groupi.append(d)
            # until d and c are the same again we want to apply the
            # ith group element to d so that we've found all possible
            # group elements for the combination
            while d != c:
                trans = [d[0][i[0][l]] for l in range(0,len(i[0]))]
                rots = [d[1][i[1][l]] for l in range(0,len(i[1]))]
                d = [trans,rots]
                # increment q because we found a new group element
                q += 1
                # every time d is not in group i add it
                if d not in groupi:
                    groupi.append(d)
    # new is a tracker that is 0 while we are still finding new group
    # elements and 1 when no new ones are found
    new = 0

    #Here each generator is applied to the group over and over until no new group
    #elements are found
    while new == 0:
        # a second set to store any additional group elements
        group2 = []
        for i in gen:
            for h in groupi:
                # apply the ith generator to the hth element of groupi
                # from above to search for new group elements
                trans = [h[0][i[0][l]] for l in range(0,len(i[0]))]
                rots = [h[1][i[1][l]] for l in range(0,len(i[1]))]
                
                d = [trans,rots]
                # if do isn't in groupi or in group2 then add it to
                # group2
                if d not in groupi:
                    if d not in group2:
                        group2.append(d)
        # add anything that is in group2 to groupi
        if len(group2) > 0:
            groupi.extend(group2)
        else:
        # if no new group elements were found then we'v found the
        # whole group and new is set to 1
            new = 1
        
    return(groupi)
