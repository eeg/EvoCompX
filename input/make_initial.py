#!/usr/bin/env python

#--------------------------------------------------
# This script is useful for creating files specifying the initial condition.
# Its output is two files, one for abundance and one for mean phenotype.
# Each of those output files has one row per spatial cell and one column per species.
#-------------------------------------------------- 

UNDEF = -9999    # this must match the value in compete.c

def list_of_lists(n):
	'''return a list containing n empty lists'''
	a = []
	for i in range(n):
		a.append([])
	return a

def get_optimum(x):
	'''return the optimum phenotype for cell x'''
	opt_slope = 0.5
	return x * opt_slope


num_cells = 100          # landscape size
num_sp = 2               # number of species

# one element for each species; start numbering cells with 0
start = [0, 99]           # first cell in range where the species is
stop = [0, 99]            # last cell in range where the species is
abun = [4.9, 5.1]         # abundance in each cell in those ranges
zoffset = [-0.1, 0.1]     # offset from optimum phenotype
  # or instead, could give fixed values for initial zbar
# (if there are many species, could use a loop to generate the above lists)

num_fp = open("num.in", "w")
zbar_fp = open("zbar.in", "w")

num = list_of_lists(num_sp)
zbar = list_of_lists(num_sp)

for sp in range(num_sp):
	for cell in range(num_cells):

		if cell < start[sp] or cell > stop[sp]:
			num[sp] = num[sp] + [0]
			zbar[sp] = zbar[sp] + [UNDEF]
		else:
			num[sp] = num[sp] + [abun[sp]]
			zbar[sp] = zbar[sp] + [ get_optimum(cell) + zoffset[sp] ]

for cell in range(num_cells):
	for sp in range(num_sp):
		num_fp.write("%f\t" % num[sp][cell])
		zbar_fp.write("%f\t" % zbar[sp][cell])
	num_fp.write("\n")
	zbar_fp.write("\n")
