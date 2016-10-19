"""This module contains methods used for visualizing the enuerated structures."""

from phenum.base import testmode
from matplotlib import cm
import matplotlib

import os
if os.name != "nt":
    matplotlib.use("Agg" if testmode else "TkAgg")    

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def HNF_shapes(enum,lattice,show,testmode=False):
    """Plots the shape of each HNF.

    :arg enum: The enum.in style input file.
    :arg lattice: The lattice.in style input file.
    :arg show: If true each HNF is plotted in an interactive window.
    """
    from phenum.io_utils import read_lattice
    from phenum.grouptheory import SmithNormalForm, get_full_HNF
    from phenum.vector_utils import map_enumStr_to_real_space, cartesian2direct
    from operator import mul
    from numpy import array, mgrid, dot
    from itertools import product

    try:
        from functools import reduce
    except ImportError:
        import numpy as np
    
    lattice_data = read_lattice(lattice)
    system = _convert_read_lat_to_system_dat(lattice_data)
    
    with open(enum,"r") as inf:
        for line in inf:
            if "#" in line:
                pass
            else:
                structure = {}
                HNF_name = [int(i) for i in line.strip().split()[:-1]]
                if len(HNF_name) != 6:
                    HNF_name = HNF_name[0:6]
                structure["HNF"] = get_full_HNF(HNF_name)
                (SNF,L,R) = SmithNormalForm(structure["HNF"])
                structure["diag"] = [SNF[0][0],SNF[1][1],SNF[2][2]]
                structure["L"] = L
                structure["n"] = reduce(mul,structure["diag"],1)
                structure["labeling"] = "".join(["0" for i in range(structure["n"])])
                system["nD"] = len(system["dvecs"])
                space_data = map_enumStr_to_real_space(system,structure,True)
                space_data["aBas"] = cartesian2direct(space_data["sLV"],space_data["aBas"],system["eps"])
                # the corners of the polyhedron except the origin.
                x = space_data["sLV"][0]
                y = space_data["sLV"][1]
                z = space_data["sLV"][2]
                xy = (array(x)+array(y)).tolist()
                xz = (array(x)+array(z)).tolist()
                yz = (array(y)+array(z)).tolist()
                xyz = (array(y)+array(x)+array(z)).tolist()

                correct = [[0,0,0],x,y,z,xy,xz,yz,xyz]
                xf = [i[0] for i in correct]
                yf = [i[1] for i in correct]
                zf = [i[2] for i in correct]

                # we need to put the shifted atoms back into cartesian corrdinates.
                atoms = []
                for atom in space_data["aBas"]:
                    atoms.append(dot(atom,space_data["sLV"]).tolist())
                        
                xa = [atom[0] for atom in atoms]
                ya = [atom[1] for atom in atoms]
                za = [atom[2] for atom in atoms]

                fig = plt.figure()
                ax = fig.gca(projection='3d')
                ax.set_aspect('equal')
                ax.set_axis_off()
                ax.scatter(xa,ya,za,zdir='z',c='red',s=1000)
                # ax.scatter(xf,yf,zf,zdir='z',c='red')

                ax.plot([0,x[0]],[0,x[1]],[0,x[2]],'k')
                ax.plot([0,y[0]],[0,y[1]],[0,y[2]],'k')
                ax.plot([0,z[0]],[0,z[1]],[0,z[2]],'k')
                ax.plot([xz[0],x[0]],[xz[1],x[1]],[xz[2],x[2]],'k')
                ax.plot([xz[0],z[0]],[xz[1],z[1]],[xz[2],z[2]],'k')
                ax.plot([xy[0],x[0]],[xy[1],x[1]],[xy[2],x[2]],'k')
                ax.plot([xy[0],y[0]],[xy[1],y[1]],[xy[2],y[2]],'k')
                ax.plot([xz[0],z[0]],[xz[1],z[1]],[xz[2],z[2]],'k')
                ax.plot([yz[0],z[0]],[yz[1],z[1]],[yz[2],z[2]],'k')
                ax.plot([yz[0],y[0]],[yz[1],y[1]],[yz[2],y[2]],'k')
                ax.plot([yz[0],xyz[0]],[yz[1],xyz[1]],[yz[2],xyz[2]],'k')
                ax.plot([xz[0],xyz[0]],[xz[1],xyz[1]],[xz[2],xyz[2]],'k')
                ax.plot([xy[0],xyz[0]],[xy[1],xyz[1]],[xy[2],xyz[2]],'k')

                max_range = array([array(xf).max()-array(xf).min(),array(yf).max()-array(yf).min(),
                                   array(zf).max()-array(zf).min()]).max()
                Xb = 0.5*max_range*mgrid[-1:2:2,-1:2:2,-1:2:2][0].flatten() + 0.5*(array(xf).max()+array(xf).min())
                Yb = 0.5*max_range*mgrid[-1:2:2,-1:2:2,-1:2:2][0].flatten() + 0.5*(array(yf).max()+array(yf).min())
                Zb = 0.5*max_range*mgrid[-1:2:2,-1:2:2,-1:2:2][0].flatten() + 0.5*(array(zf).max()+array(zf).min())

                for xb, yb, zb in zip(Xb,Yb,Zb):
                    ax.plot([xb],[yb],[zb],'w')

                if not testmode: # pragma: no cover
                    fig.savefig("{}.pdf".format("".join([str(i) for i in HNF_name])))
                if show and not testmode: #pragma: no cover
                    plt.show()
                else:
                    plt.close()
    
def HNF_atoms(enum,lattice,show,testmode=False):
    """Plots the atomic positions of the atoms in the cells.

    :arg enum: The enum.in style input file.
    :arg lattice: The lattice.in style input file.
    :arg show: If true each HNF is plotted in an interactive window.
    """
    from phenum.io_utils import read_lattice
    from phenum.grouptheory import SmithNormalForm, get_full_HNF
    from phenum.vector_utils import map_enumStr_to_real_space, cartesian2direct
    from operator import mul
    from numpy import array, mgrid, dot
    from itertools import product
    
    try:
        from functools import reduce
    except ImportError:
        import numpy as np

    lattice_data = read_lattice(lattice)
    system = _convert_read_lat_to_system_dat(lattice_data)
    
    with open(enum,"r") as inf:
        for line in inf:
            if "#" in line:
                pass
            else:
                structure = {}
                HNF_name = [int(i) for i in line.strip().split()[:-1]]
                if len(HNF_name) != 6:
                    HNF_name = HNF_name[0:6]
                structure["HNF"] = get_full_HNF(HNF_name)
                (SNF,L,R) = SmithNormalForm(structure["HNF"])
                structure["diag"] = [SNF[0][0],SNF[1][1],SNF[2][2]]
                structure["L"] = L
                structure["n"] = reduce(mul,structure["diag"],1)
                structure["labeling"] = "".join(["0" for i in range(structure["n"])])
                system["nD"] = len(system["dvecs"])
                space_data = map_enumStr_to_real_space(system,structure,True)
                space_data["aBas"] = cartesian2direct(space_data["sLV"],space_data["aBas"],system["eps"])
                atoms = []
                for atom in space_data["aBas"]:
                    atoms.append(dot(atom,space_data["sLV"]).tolist())
                        
                xf = [atom[0] for atom in atoms]
                yf = [atom[1] for atom in atoms]
                zf = [atom[2] for atom in atoms]

                fig = plt.figure()
                ax = fig.gca(projection='3d')
                ax.set_aspect('equal')
                ax.set_axis_off()
                ax.scatter(xf,yf,zf,zdir='z',c='red',s=1000)

                max_range = array([array(xf).max()-array(xf).min(),array(yf).max()-array(yf).min(),
                                   array(zf).max()-array(zf).min()]).max()
                Xb = 0.5*max_range*mgrid[-1:2:2,-1:2:2,-1:2:2][0].flatten() + 0.5*(array(xf).max()+array(xf).min())
                Yb = 0.5*max_range*mgrid[-1:2:2,-1:2:2,-1:2:2][0].flatten() + 0.5*(array(yf).max()+array(yf).min())
                Zb = 0.5*max_range*mgrid[-1:2:2,-1:2:2,-1:2:2][0].flatten() + 0.5*(array(zf).max()+array(zf).min())

                for xb, yb, zb in zip(Xb,Yb,Zb):
                    ax.plot([xb],[yb],[zb],'w')

                if not testmode: #pragma: no cover
                    fig.savefig("{}.pdf".format("".join([str(i) for i in HNF_name])))
                if show and not testmode:#pragma: no cover
                    plt.show()
                else:
                    plt.close()

def _convert_read_lat_to_system_dat(lattice):
    """This method converts the lattice_data dictionary returned by
    io_utils.read_lattice to the dictionary structure of io_utils.read_enum_out.

    :arg lattice: The io_utils.read_lattice style dictionary.
    """

    system_data = {
        "title": "HNF plotter",
        "bulksurf": lattice["bulk"],
        "plattice": lattice["lat_vecs"],
        "nD": 0,
        "dvecs": lattice["basis_vecs"],
        "k": lattice["nspecies"],
        "eps": 1E-10
    }

    return system_data
