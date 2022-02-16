#!/usr/bin/env python

#--------------------------------------------------
# This script is useful for creating a file specifying a simple matrix of
#   competition coefficients. 
# Its output is a file with num_species rows and columns, and whitespace
#   separating the elements.
#
# If a number is given on the command line, it is used for all the off-
#   diagonal elements.
# Otherwise, the off-diagonal elements are chosen randomly (uniform on [0, 1)).
#-------------------------------------------------- 

import sys
import random
rand = random.Random()

num_sp = 2               # number of species

try:
	a_ij = sys.argv[1]
except IndexError:
	a_ij = None

fp = open("alpha.in", "w")

if a_ij:
	for i in range(num_sp):
		for j in range(num_sp):
			if i == j:
				fp.write("1.00  ")
			else:
				fp.write("%s  " % a_ij)
		fp.write("\n")

else:
	for i in range(num_sp):
		for j in range(num_sp):
			if i == j:
				fp.write("1.00  ")
			else:
				fp.write("%0.2f  " % rand.uniform(0., 1.))
		fp.write("\n")

fp.close()
