import matplotlib.pyplot as plt
import numpy as np
import math
import sys

def calculate_pi_naught(profile_number):
    # Read in radii and brunt-vasala freqs
    r, bvf, flags = np.genfromtxt(
        fname="../profile{}.data".format(profile_number),
        unpack=True, usecols=[4,97,13], skip_header=6)
    
    r = 10**r

    sum = 0
    for i in range(1,len(r)):
        if flags[i] == 0:
            left  = bvf[i-1] / r[i-1]
            right = bvf[i]   / r[i]
            midpoint = (left + right) / 2
            
            dr = abs(r[i] - r[i-1])
            
            sum = sum + midpoint * dr

    pi_0 = 2 * (math.pi**2) / sum
    return pi_0

profiles = np.loadtxt(fname="Data/filteredProfileNumbers.dat")

pi_1 = []
models = []
for i in range(len(profiles)):
    
    print(" {}% Complete...         ".format(int(i/len(profiles)*100)),end="\r")
    
    models.append(profiles[i])
    pi = calculate_pi_naught(int(profiles[i])) / (2**0.5)
    pi_1.append(pi)

with open("Data/DPi_1.dat", "w") as file:
    file.write("Model_Number\tDelta_Pi_1\n")
    for i in range(len(pi_1)):
        file.write(str(i+1) + "\t\t" + str(pi_1[i]) + "\n")
