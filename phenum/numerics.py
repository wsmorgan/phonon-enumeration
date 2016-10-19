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
        r = list(range(1,num+1))
        for i in r:
            fact *= i
    return fact

def binomial_coefficient(n, k):

    if k < 0 or k > n:
        return 0
    if k==0 and n == 0:
        return 1
    t = 1
    if k < n-k:
        for i in range(n, n-k, -1):
            t = t*i//(n-i+1)
    else:
        for i in range(n, k, -1):
            t = t*i//(n-i+1)

    return t
