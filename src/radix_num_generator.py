#radix_num_generator.py
#import needed modules
import binomial_calculator as bc
import phonon_brancher as pb
from copy import deepcopy

#This method generates the radix numbers for a given aray. It calls
#for number_of_each_color and number_of_slots.
#number_of_each_color is an array of integers indicating the number of
#times each color is found in the array.
#number_of_slots is an integer, it should be equal to the sum of the
#number_of_each_color array.
#It returns an array of radix numbers
def Coefficients(number_of_each_color,number_of_slots):
        #renaims variables for convience and creates an empty array
	Coe = []
	for i in range(0,len(number_of_each_color)):
        #Uses the binomial_coefficient method to build an array of radix numebrs
		Coe.append(bc.binomial_coefficient(number_of_slots,number_of_each_color[i]))
		number_of_slots -= number_of_each_color[i]
	return Coe

#This method generates numerical labels for each color in the provided
#array. It calls for listi and number_of_each_color, and returns an
#array of integers.
#listi is an array of integers it is the array that will be used to
#generate the labels.
#color is and array of integers that indicates how many
#times each color is found in the array
def hash(listi,color):
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

#ahash produces a unique number for each possible configuration. This
#number is used to compare the order that the arrows occure in.
#coloring is the full coloring with colors and arrows that get passed
#in
#arrows is the unique number that is returned by the program.
#dim is the number of directions the arrows can point.
def ahash(coloring,dim):

    narrow = 0
    for i in range(len(coloring)):
        narrow = narrow + coloring[i]*dim**i  
    return(narrow)

#ainvhash takes a number and the number of arrows and returns the
#arrow array associated with it.
#anum is a unique number that is associated with an array of arrows.
#num_of_arrows in the number of arrows that are in the array.
#dim is the number of directions the arrows can point.
def ainvhash(anum,num_of_arrrows,dim):
        arrows = [0]*num_of_arrrows
        for i in range(num_of_arrrows):
                base = dim**(num_of_arrrows-1-i)
                arrows[num_of_arrrows-1-i] = anum//base
                anum = anum -base*(anum//base)
        return(arrows)
