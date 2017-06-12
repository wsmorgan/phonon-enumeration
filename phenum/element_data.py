# The dictionary of all the elements on the periodic table
all_elements ={"H": 3.75,"He": 3.57,"Li": 3.49,"Be": 2.29,"B": 8.73,"C": 3.57,"N": 4.039,
               "O": 6.83,"Ne": 4.43,"Na": 4.23,"Mg": 3.21,"Al": 4.05,"Si": 5.43,"P": 7.17,
               "S": 10.47,"Cl": 6.24,"Ar": 5.26,"K": 5.23,"Ca": 5.58,"Sc": 3.31,"Ti": 2.95,
               "V": 3.02,"Cr": 2.88,"Mn": 8.89,"Fe": 2.87,"Co": 2.51,"Ni": 3.52,"Cu": 3.61,
               "Zn": 2.66,"Ga": 4.51,"Ge": 5.66,"As": 4.13,"Se": 4.36,"Br": 6.67,"Kr": 5.72,
               "Rb": 5.59,"Sr": 6.08,"Y": 3.65,"Zr": 3.23,"Nb": 3.3,"Mo": 3.15,"Tc": 2.74,
               "Ru": 2.7,"Rh": 3.8,"Pd": 3.89,"Ag": 4.09,"Cd": 2.98,"In": 4.59,"Sn": 5.82,
               "Sb": 4.51,"Te": 4.45,"I": 7.27,"Xe": 6.2,"Cs": 6.05,"Ba": 5.02,"Hf": 3.2,
               "Ta": 3.31,"W": 3.16,"Re": 2.76,"Os": 2.64,"Ir": 3.84,"Pt": 3.92,"Au": 4.08,
               "Hg": 2.99,"Tl": 3.46,"Pb": 4.95,"Bi": 4.75}

def get_lattice_parameter(elements, concentrations, default_title):
    """Finds the lattice parameters for the provided atomic species using Vagars law.

    Args:
        elements (list of str): A list of the elements in the system.
        default_title (str): The default system title.
        concentrations (list of int): The concentrations of each element.

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
                lat_param += concentrations[i]*all_elements[elem]
                if concentrations[i] > 0:
                    title += " {0} ".format(elem)
            lat_param = float(lat_param) / sum(concentrations)
            title = "{0} {1}\n".format(default_title.strip(),title)
    return lat_param, title
