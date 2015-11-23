#/bin/py
import polyaburnside as burn
import phonon_brancher as pb
import arrow_group as ar
from copy import deepcopy
import os

trans = [[0,1,2,3],[1,0,3,2],[2,3,0,1],[3,2,1,0]]
rots = [[[0,1,2,3],[0,1,2,3]],[[0,1,2,3],[2,1,0,3]],[[0,1,2,3],[0,3,2,1]],[[0,1,2,3],[2,3,0,1]],[[0,3,2,1],[1,0,3,2]],[[0,3,2,1],[3,0,1,2]],[[0,3,2,1],[1,2,3,0]],[[0,3,2,1],[3,2,1,0]]]

group = ar.a_group(trans,rots)

colors = [[[-1,1],[-1,2],[0,3],[0,4]],[[-1,1],[-1,1],[-1,1],[-1,1]],[[0,1],[0,1],[0,1],[0,1]],[[-1,1],[-1,1],[-1,2],[0,3]],[[-1,1],[0,2],[0,2],[0,3]],[[-1,1],[-1,2],[0,3],[0,3]],[[-1,1],[-1,1],[0,2],[0,3]],[[-1,1],[-1,1],[-1,2],[-1,2]],[[-1,1],[-1,1],[0,2],[0,2]],[[-1,1],[-1,1],[-1,1],[-1,2]],[[-1,1],[-1,1],[-1,1],[0,4]],[[-1,1],[0,2],[0,2],[0,2]],[[0,1],[0,1],[0,1],[0,2]],[[0,1],[0,1],[0,2],[0,2]]]

for i in range(len(colors)):
    tcolors = deepcopy(colors)
    tgroup = deepcopy(group)
    direc = "group.in." +str(i+1)
    colf = "colors.in." + str(i+1)
    Concf = "concentrations.in." + str(i+1)
    naf ="number_of_arrows.in." + str(i+1)
    pf = "polya_burnside.out." + str(i+1)
    
    
    (n_arrows,arrow_types) = pb.how_many_arrows(tcolors[i])
    Concs = pb.find_concentrations(tcolors[i])

    print(Concs, arrow_types)
    print(burn.polya(Concs,group,arrowings=arrow_types))

    f = open(pf, 'w+')
    f.write(str(burn.polya(Concs,group,arrowings=arrow_types)))
    f.close()

    if not os.path.exists(direc):
        os.makedirs(direc)
        
    tgroup = deepcopy(group)
    tgroup = [[[j + 1 for j in l] for l in t] for t in tgroup]
    z = 1    
    for t in tgroup:
        transf = "_-" + str(z) + "-site_perm"
        f = open(direc + "/" + transf, 'w+')
        for j in t[0]:
            j = '{0:g}'.format(j)
            f.write(str(j) + ' ')
        f.write('\n')
        z += 1
    f.close()

    tgroup = deepcopy(group)
    tgroup = [[[j + 1 for j in l] for l in t] for t in tgroup]
    
    z = 1
    for t in tgroup:
        rotf = "_-" + str(z) + "-arrow_perm"
        f = open(direc + "/" + rotf, 'w+')
        for j in (t[1]):
            j = '{0:g}'.format(j)
            f.write(str(j) + ' ')
        f.write('\n')
        z += 1
    f.close()

    tcolors = deepcopy(colors)

    f = open(colf, 'w+')
    print(tcolors,i)
    for t in tcolors[i]:
        for j in range(len(t)):
            t[j] = '{0:g}'.format(t[j])
            f.write(str(t[j]) + ' ')
        f.write('\n')
    f.close()

    f = open(Concf, 'w+')
    for j in range(len(Concs)):
        Concs[j] = '{0:g}'.format(Concs[j])
        f.write(str(Concs[j]) + ' ')
    f.close()

    f = open(naf, 'w+')
    arrow_types = '{0:g}'.format(arrow_types)
    f.write(str(arrow_types) + ' ')
    f.close()
    
i = 15

for this_group in group:
    colf = "colors.in." + str(i)
    direc = "group.in." + str(i)
    Concf = "concentrations.in." + str(i)
    naf ="number_of_arrows.in." + str(i)
    pf = "polya_burnside.out." + str(i)

    tgroup = [deepcopy(this_group)]

    f = open(pf, 'w+')
    (n_arrows,arrow_types) = pb.how_many_arrows(colors[0])
    Concs = pb.find_concentrations(colors[0])
    print(Concs, tgroup)
    f.write(str(burn.polya(Concs,tgroup,arrowings=arrow_types)))
    f.close()

    if not os.path.exists(direc):
        os.makedirs(direc)

    tgroup = [deepcopy(this_group)]
    tgroup = [[[j + 1 for j in l] for l in t] for t in tgroup]
    print(tgroup)
    z = 1
    for t in tgroup:
        transf = "_-"+str(z)+"-site_perm"        
        f = open(direc + "/" + transf, 'w+')
        for j in t[0]:
            j = '{0:g}'.format(j)
            f.write(str(j) + ' ')
        f.write('\n')
        z += 1
    f.close()

    tgroup = [deepcopy(this_group)]
    tgroup = [[[j + 1 for j in l] for l in t] for t in tgroup]

    z = 1
    for t in tgroup:
        rotf = "_-"+str(z)+"-arrow_perm"
        f = open(direc + "/" + rotf, 'w+')
        for j in (t[1]):
            j = '{0:g}'.format(j)
            f.write(str(j) + ' ')
        f.write('\n')
        z += 1
    f.close()

    tcolors = deepcopy(colors)
    f = open(colf, 'w+')
    for t in tcolors[0]:
        for j in range(len(t)):
            t[j] = '{0:g}'.format(t[j])
            f.write(str(t[j]) + ' ')
        f.write('\n')
    f.close()

    f = open(Concf, 'w+')
    for j in range(len(Concs)):
        Concs[j] = '{0:g}'.format(Concs[j])
        f.write(str(Concs[j]) + ' ')
    f.close()

    f = open(naf, 'w+')
    arrow_types = '{0:g}'.format(arrow_types)
    f.write(str(arrow_types) + ' ')
    f.close()
    i += 1

tt = 47

trans = [[j - 1 for j in i] for i in[[1, 2, 3, 4, 5, 6, 7, 8], [2, 1, 4, 3, 6, 5, 8, 7], [3, 4, 5, 6, 7, 8, 1, 2], [4, 3, 6, 5, 8, 7, 2, 1], [5, 6, 7, 8, 1, 2, 3, 4], [6, 5, 8, 7, 2, 1, 4, 3], [7, 8, 1, 2, 3, 4, 5, 6], [8, 7, 2, 1, 4, 3, 6, 5]]]
rots = [[[j - 1 for j in i] for i in t] for t in [[[1, 2, 3, 4, 5, 6, 7, 8], [1, 2, 3, 4, 5, 6]], [[1, 2, 3, 4, 5, 6, 7, 8], [4, 2, 3, 1, 5, 6]], [[1, 2, 3, 4, 5, 6, 7, 8], [1, 5, 3, 4, 2, 6]], [[1, 2, 3, 4, 5, 6, 7, 8], [4, 5, 3, 1, 2, 6]], [[1, 2, 7, 8, 5, 6, 3, 4], [1, 2, 6, 4, 5, 3]], [[1, 2, 7, 8, 5, 6, 3, 4], [4, 2, 6, 1, 5, 3]], [[1, 2, 7, 8, 5, 6, 3, 4], [1, 5, 6, 4, 2, 3]], [[1, 2, 7, 8, 5, 6, 3, 4], [4, 5, 6, 1, 2, 3]]]]

group = ar.a_group(trans,rots)

colors= [[[-1, 1], [-1, 1], [-1, 2], [-1, 2], [-1, 2], [0, 3], [0, 4], [0, 4]],[[-1, 1], [-1, 1], [-1, 2], [-1, 2], [0, 2], [0, 3], [0, 4], [0, 4]],[[-1, 1], [-1, 2], [-1, 3], [0, 2], [0, 2], [0, 3], [0, 4], [0, 4]],[[-1, 1], [-1, 1], [-1, 2], [-1, 2], [-1, 2], [-1, 3], [-1, 4], [-1, 4]]]

for i in range(len(colors)):
    tcolors = deepcopy(colors)
    tgroup = deepcopy(group)
    direc = "group.in." +str(tt)
    colf = "colors.in." + str(tt)
    Concf = "concentrations.in." + str(tt)
    naf ="number_of_arrows.in." + str(tt)
    pf = "polya_burnside.out." + str(tt)
    
    
    (n_arrows,arrow_types) = pb.how_many_arrows(tcolors[i])
    Concs = pb.find_concentrations(tcolors[i])

    print(Concs, arrow_types)
    print(burn.polya(Concs,group,arrowings=arrow_types))

    f = open(pf, 'w+')
    f.write(str(burn.polya(Concs,group,arrowings=arrow_types)))
    f.close()

    if not os.path.exists(direc):
        os.makedirs(direc)
        
    tgroup = deepcopy(group)
    tgroup = [[[j + 1 for j in l] for l in t] for t in tgroup]
    z = 1    
    for t in tgroup:
        transf = "_-" + str(z) + "-site_perm"
        f = open(direc + "/" + transf, 'w+')
        for j in t[0]:
            j = '{0:g}'.format(j)
            f.write(str(j) + ' ')
        f.write('\n')
        z += 1
    f.close()

    tgroup = deepcopy(group)
    tgroup = [[[j + 1 for j in l] for l in t] for t in tgroup]
    
    z = 1
    for t in tgroup:
        rotf = "_-" + str(z) + "-arrow_perm"
        f = open(direc + "/" + rotf, 'w+')
        for j in (t[1]):
            j = '{0:g}'.format(j)
            f.write(str(j) + ' ')
        f.write('\n')
        z += 1
    f.close()

    tcolors = deepcopy(colors)

    f = open(colf, 'w+')
    print(tcolors,i)
    for t in tcolors[i]:
        for j in range(len(t)):
            t[j] = '{0:g}'.format(t[j])
            f.write(str(t[j]) + ' ')
        f.write('\n')
    f.close()

    f = open(Concf, 'w+')
    for j in range(len(Concs)):
        Concs[j] = '{0:g}'.format(Concs[j])
        f.write(str(Concs[j]) + ' ')
    f.close()

    f = open(naf, 'w+')
    arrow_types = '{0:g}'.format(arrow_types)
    f.write(str(arrow_types) + ' ')
    f.close()
    tt += 1
