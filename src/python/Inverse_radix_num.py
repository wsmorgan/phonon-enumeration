#Inverse_radix_num.py

#import needed modules
import binomial_calculator as bc
import numpy as np  

#This module takes a number and turns it into an array of numbers
#based off of the number and other values passed to it.
#branch is an integer list of the input radix number that needs to
#have the hashing undone.
#Colors is an array containing the concentrations of each color to be
#used.
#n is the length of the new array.
#Reurns an array filled with colors.
def invhash(branch,colors,n):
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
