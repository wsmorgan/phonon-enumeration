"""Methods used to turn a labeling into a hash, or radix number and back again."""

#import needed modules
from copy import deepcopy
import numpy as np

def coefficients(number_of_each_color,number_of_slots):
    """This method generates the radix numbers, or total number of
    combinations of each color, for a given aray.

    :arg number_of_each_color: is an integer array integer indicating the
    number of times each color is found in the array.
    :arg number_of_slots: an integer equal to the sum of the
    number_of_each_color array.
    """

    #renames variables for convience and creates an empty array
    coe = []
    for i in range(0,len(number_of_each_color)):
    #Uses the binomial_coefficient method to build an array of radix numebrs
	coe.append(binomial_coefficient(number_of_slots,number_of_each_color[i]))
	number_of_slots -= number_of_each_color[i]
    return coe

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
	y += binomial_coefficient(len(li),k-1)
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

def invhash(branch,colors,n):
    """Turns a hash array into a labeling, i.e., undoes the hash()
    subroutine.
    
    :arg branch: an integer list of the radix number/hash array
    :arg colors: an integer array of the concentration of the colors
    :arg n: the length of the labeling array, i.e., the
    number of sites in the system

    For the details on this method see: http://msg.byu.edu/papers/enum3.pdf
    """

    new = [0 for y in range(n)]
    coluse = 1
    for i in range(0,len(branch)):
        newi = [0 for y in range(n)]
        I = branch[i]
        t = colors[i]
        for l in range(n,0,-1):
            binom = binomial_coefficient(l-1, t-1)
            if binom <= I:
                I -= binom
            else:
                newi[n-l] = coluse
                t -= 1
        coluse += 1
        n -= colors[i]
        si = 0
        for s in range(0,len(new)):
            if new[s] == 0:
                new[s] = newi[si]
                si += 1
    return(new)

#This method finds the factorial of a given number(num).
#The method calls for an integer input num
#The method returns the factorial(represented as fact in the method)
def factorial(num):
    """Finds the factorial of the input integer.
    
    :arg num: an integer
    """
    #If the number provided is zero then the factorial is 1
    if num == 0:
	fact = 1
    #Otherwise set fact to 1 and begin finding the factorial r is
    #used to find each num-n for n=0 to n=num each value in r is
    #then used to compute the factorial
    else:
	fact = 1
	r = range(1,num+1)
	for i in r:
	    fact *= i
    return fact

#This method finds a binomial coefficient that is it calculates n
#choose r The method calls for two integer values n and r and the
#method factorial The method returns the binomial factorial as result
def binomial_coefficient(n,r):
    """Finds the binomial coefficient, n choose r, for a given set of
    integers.
    
    :arg n: an integer
    :arg r: an integer
    """
    #If r is less than zero then the binomial coefficient is zero
    if r < 0:
	result = 0
    #Otherwise use the factorial method to calculate the binomial
    #coefficient n!/r!/(n-r)!
    else:
	result = factorial(n) / factorial(r) / factorial(n-r)
    return result
