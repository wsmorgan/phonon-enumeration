"""Methods for turning a hash into a labeling"""

#import needed modules
import binomial_calculator as bc
import numpy as np  

def invhash(branch,colors,n):
        """Turns a hash array into a labeling, i.e., undoes the hash function
          from radix_number_generator.py.

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
                        binom = bc.binomial_coefficient(l-1, t-1)
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
