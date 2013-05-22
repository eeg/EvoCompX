#!/usr/bin/env python

#--------------------------------------------------
# This script is useful for creating files specifying the initial condition.
# Its output is two files, one for abundance and one for mean phenotype.
# Each output file has one row per spatial cell and one column per species.
#-------------------------------------------------- 

# TODO: pull params from same input file that EvoCompX uses

UNDEF_PHEN = -9999    # this must match the value in input.h

def list_of_lists(n):
	'''return a list containing n empty lists'''
	a = []
	for i in range(n):
		a.append([])
	return a

def get_best_abar(x):
    '''return the optimum mean breeding value for cell x'''
    B = 0.5           # slope of the optimum phenotype function
    C = 1             # slope of the environment function
    bbar = 0.1        # mean plasticity
    return (B - bbar) * C * x


num_cells = 100          # landscape size
num_sp = 2               # number of species

# one element for each species; start numbering cells with 0
start = [0, 90]          # first cell in range where the species is
stop  = [9, 99]          # last cell in range where the species is
abun  = [4.0, 6.0]       # abundance in each cell in those ranges
abar_off = [-0.5, 0.3]  # offset from optimum mean breeding value

num_fp = open("num.in", "w")
abar_fp = open("abar.in", "w")

num  = list_of_lists(num_sp)
abar = list_of_lists(num_sp)

for sp in range(num_sp):
	for cell in range(num_cells):

		if cell < start[sp] or cell > stop[sp]:
			num[sp] = num[sp] + [0]
			abar[sp] = abar[sp] + [UNDEF_PHEN]
		else:
			num[sp] = num[sp] + [abun[sp]]
			abar[sp] = abar[sp] + [ get_best_abar(cell) + abar_off[sp] ]

for cell in range(num_cells):
	for sp in range(num_sp):
		num_fp.write("%f\t" % num[sp][cell])
		abar_fp.write("%f\t" % abar[sp][cell])
	num_fp.write("\n")
	abar_fp.write("\n")
