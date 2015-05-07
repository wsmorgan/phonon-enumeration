#Inverse_radix_num.py

#import needed modules
import binomial_calculator as bc
import numpy as np 
import radix_num_generator as rng 

#This module takes a number and turns it into an array of numbers based off of the number and
#other values passed to it. Coefficients is an array of radix numbers to be used in the algorithm.
#Coloros is an array indicating how many of each color will be used in the array. Radix_number is the 
#integer number that will be converted into an array. Number_of_slots is the length of the new array.
#Reurns an array filled with colors. 
def invid(Coefficient,colors,radix_number,number_of_slots):
#Set the inputs to easily used variables and creat an emoty array of length m.
        new = [0 for y in range(number_of_slots)]#np.zeros(number_of_slots)
#rev is used as a counter and coluse is a counter for the numeral to represent each color
        coluse = 1
        for rev in range(0,len(Coefficient)):
#For each go through the loop create a temporary array of length m
                newi = [0 for f in range(number_of_slots)] #np.zeros(number_of_slots)
#The variables t, I, y are defiend for later use in the algorithm and l is a duplicate of m 
                t = colors[rev]
                I = radix_number%Coefficient[rev]
                radix_number = radix_number//Coefficient[rev]
                for l in range(number_of_slots,0,-1):
                        #Performs loops for each index in the new array to determine if that index needs to be filled
                        binom = bc.binomial_coefficient(l-1, t-1)
                        if binom <= I:
                                I -= binom
                        else:
                                newi[number_of_slots-l] = coluse
                                t -= 1
                coluse += 1
                #Changes the value of m so a new array can be generated from only the empty parts of newi
                number_of_slots -= colors[rev]
                si = 0
#This loop fills in the blanks in the new array for the colors from array newi for each color
                for s in range(0,len(new)):
                        if new[s] == 0:
                                new[s] = newi[si]
                                si += 1
        return(new)

def binomial(I,n,col):
	new = [0 for y in range(n)]#np.zeros(n)
	for l in range(n,0,-1):
#Performs loops for each index in the new array to determine if that index needs to be filled
		if bc.binomial_coefficient(l-1, col-1) <= I:
			I = I - bc.binomial_coefficient(l-1, col-1)
		else:
			new[n-l] = 1	
			col -= 1
	return new

def invhash(branch,colors,n):
        new = [0 for y in range(n)]#np.zeros(n)
        coluse = 1
        for i in range(0,len(branch)):
                newi = [0 for y in range(n)]#np.zeros(
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

def invhashbin(y,m,a):
        I = y
        l = m
        t = a
        result = [0]*m
        while l > 0:
                binom = bc.binomial_coefficient(l-1,t-1) 
                if binom <= I:
                        result[m-l] = 0
                        I -= binom
                else:
                        result[m-l] = 0
                        t += 1
                l -= 1
        return(result)

#print(binomial(3,10,2))
#print(invhash([0,3,0],[4,2,2],8))
