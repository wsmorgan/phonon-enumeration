"""Methods used to turn a labeling into a hash, or radix number"""

#import needed modules
import binomial_calculator as bc
import phonon_brancher as pb
from copy import deepcopy

def Coefficients(number_of_each_color,number_of_slots):
        """This method generates the radix numbers, or total number of
          combinations of each color, for a given aray.

          :arg number_of_each_color: is an integer array integer indicating the
          number of times each color is found in the array.
          :arg number_of_slots: an integer equal to the sum of the
          number_of_each_color array.

        """

        #renames variables for convience and creates an empty array
	Coe = []
	for i in range(0,len(number_of_each_color)):
        #Uses the binomial_coefficient method to build an array of radix numebrs
		Coe.append(bc.binomial_coefficient(number_of_slots,number_of_each_color[i]))
		number_of_slots -= number_of_each_color[i]
	return Coe

def hash(listi,color):
        """Generates a hash for each color in the ladeling returning a array
          of the hash for each color.

          :arg listi: An integer array that contains the labeling.
          :arg color: An integer array of the concentrations of the
          colors.

          For details on this method see: http://msg.byu.edu/papers/enum3.pdf
        """
        m = len(listi)
        li = listi[::-1]
        rm1 = m - li.index(1)
        y = 0
        z = listi[:rm1].count(0)#m - listi[:rm1].count(1)
        li = listi
        for i in range(z):
                p0 = li.index(0) #next(x[0] for x in enumerate(li) if x[1] > j)
                k = li[p0:].count(1)
                li = li[p0+1:]
                y += bc.binomial_coefficient(len(li),k-1)
        return y

def ahash(coloring,dim):
    """Produces a unique number for each possible configuration of
      arrows. This number is used to compare the order that the arrows
      occure in.

      :arg coloring: is a 2D integer array of the full coloring with colors and arrows      
      :arg dim: the number of directions the arrows can point
    """
        
    narrow = 0
    for i in range(len(coloring)):
        narrow = narrow + coloring[i]*dim**i
    return(narrow)

#anum is a unique number that is associated with an array of arrows.
#num_of_arrows in the number of arrows that are in the array.
#dim is the number of directions the arrows can point.
def ainvhash(anum,num_of_arrrows,dim):
        """Turns an arrow hash back into the array of arrow directions.

          :arg anum: the arrow hash number
          :arg num_of_arrows: the number of arrows in the system
          :arg dim: the number of directions the arrows can point
        """
        arrows = [0]*num_of_arrrows
        for i in range(num_of_arrrows):
                base = dim**(num_of_arrrows-1-i)
                arrows[num_of_arrrows-1-i] = anum//base
                anum = anum -base*(anum//base)
        return(arrows)
