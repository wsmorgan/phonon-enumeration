#generators.py

#group takes a list of generators and applies them to generate the group. For
#arguments it takes gen, the list of generators as an array. It then applies
#the generaters to each other to make the element c, if c is unique it is saved
# as part of the group. the generaters are thec to make d until d = c and each 
#unique d is saved.  It returns the group.
def group(gen):
    q = 1
    groupi = []
    for i in gen:
        gen2 = []      
        for k in gen:
            if k != i:
                gen2.append(k)       
        for j in gen2:
            c = [j[i[l]] for l in range(0,len(i))]
            if c not in groupi:
                groupi.append(c)
            d = [c[i[l]] for l in range(0,len(i))]
            if d not in groupi:
                groupi.append(d)
            while d != c:
                d = [d[i[l]] for l in range(0,len(i))]
                q += 1
                if d not in groupi:
                    groupi.append(d)
    new = 0

#Here each generator is appnied to the group over and over until no new group
#elements are found
    while new == 0:
        group2 = []
        for i in gen:
            for h in groupi:
                d = [h[i[l]] for l in range(0,len(i))]
                if d not in groupi:
                    if d not in group2:
                        group2.append(d)

        if len(group2) > 0:
            groupi.extend(group2)
        else:
            new = 1
        
    return(groupi)

