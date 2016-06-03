"""Numerical methods for mathematical functions neede for the program"""

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
