import numpy as np

# The dictionary of all the elements on the periodic table
element_volume ={"H":37.2958,"He":32.1789,"Li":21.2543,"Be":8.49323,"B":7.24205,"C":5.68741,
                 "N":46.6002,"O":22.2802,"F":17.0258,"Ne":21.7346,"Na":23.2596,"Mg":23.3928,
                 "Al":16.6075,"Si":7.8511,"P":9.1459,"S":17.1672,"Cl":35.2074,"Ar":36.3829,
                 "K":71.5278,"Ca":43.4353,"Sc":25.6478,"Ti":18.1565,"V":13.7718,"Cr":11.9439,
                 "Mn":19.3207,"Fe":11.82,"Co":11.1838,"Ni":10.9036,"Cu":11.7615,"Zn":13.311,
                 "Ga":18.4496,"Ge":19.3638,"As":20.4270,"Se":58.6173,"Br":33.3170,"Kr":46.7873,
                 "Rb":87.3384,"Sr":56.1889,"Y":33.0792,"Zr":23.8327,"Nb":17.9685,"Mo":15.6279,
                 "Tc":14.5458,"Ru":13.9206,"Rh":13.718,"Pd":14.716,"Ag":17.1045,"Cd":18.7161,
                 "In":26.6861,"Sn":29.3238,"Sb":27.1733,"Te":62.3227,"I":24.3807,"Xe":59.582,
                 "Cs":110.723,"Ba":63.253,"Hf":23.1748,"Ta":18.1323,"W":15.7772,"Re":14.8694,
                 "Os":14.5485,"Ir":14.1558,"Pt":15.0591,"Au":16.9793,"Hg":27.6914,"Tl":29.2949,
                 "Pd":30.3218,"Bi":31.2849}

def get_lat_param_element(lat_vecs,n_basis_atoms,element):
    """Finds the lattice parameter of an element for the given set of lattice.

    Args:
        lat_vecs (list): The lattice vectors for the system.
        n_basis_atoms (int): The number of atoms in the atomic basis.
        element (str): The symbol for the element.
    
    Returns:
        lat_param (float): The lattice parameter.
    """

    lat_vol = abs(np.linalg.det(lat_vecs)/n_basis_atoms)
    atom_vol = element_volume[element]/lat_vol

    return atom_vol**(1./3.)   

def get_lattice_parameter(elements, concentrations, lat_vecs, n_basis_atoms, default_title,remove_zeros=False):
    """Finds the lattice parameters for the provided atomic species using Vagars law.

    Args:
        elements (list of str): A list of the elements in the system.
        default_title (str): The default system title.
        concentrations (list of int): The concentrations of each element.
        lat_vecs (list): The lattice vectors for the system.
        n_basis_atoms (int): The number of atoms in the atomic basis.
        remove_zeros (bool): True if zeros are to be removed from the elements list.

    Returns:
        lat_param (float): The lattice parameter of the system.

    Raises:
        ValueError: if the number of elements doesn't match the len of the concentration array.
    """

    if elements is None:
        lat_param = 1.0
        title = default_title
    else:
        if len(elements) != len(concentrations):
            raise ValueError("You have provided {0} element names when {1} elements are present "
                "in the system. Please provide the correct number of elements."
                .format(len(elements),len(concentrations)))

        else:
            title = ""
            lat_param = 0
            for i, elem in enumerate(elements):
                lat_param += concentrations[i]*get_lat_param_element(lat_vecs,n_basis_atoms,elem)
                if concentrations[i] > 0:
                    title += " {0} ".format(elem)
                elif not remove_zeros:
                    title += " {0} ".format(elem)
            lat_param = float(lat_param) / sum(concentrations)
            title = "{0} {1}\n".format(title,default_title.strip())
    return lat_param, title
