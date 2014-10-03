# fitting a triplet of library residues into the backbone

import os
import re
import random
from routines import *
import matplotlib.pyplot as plt
import pylab
from mpl_toolkits.mplot3d import Axes3D
from collections import deque
from auction_algo import *


#name_file = '/1CPZ.mr.txt'
name_file = '/1G6J.mr.txt' # ubiquitin

nodes, distance_matrix = read_restraints(name_file)

def dispersion_score(u,v,distance_matrix):
	dispersion_score = 0.0
	
	set_u, set_v = set(), set()
	
	for node in distance_matrix[u]:
		if node!=u:
			set_u.add(node)
	for node in distance_matrix[v]:
		if node!=v:
			set_v.add(node)
	set_commun = set_u.intersection(set_v)
	
	if len(set_commun)==0:
		try:
			distance_ = distance_matrix[u][v]
			dispersion_score = (dispersion_score+1.)/float(len(set_u)*len(set_v))
			return dispersion_score
		except KeyError:
			return dispersion_score
	else:
		sum_local = 0.5
		for node_1 in set_commun:
			for node_2 in set_commun:
				try:
					distance_ = distance_matrix[node_1][node_2]
					sum_local+=0.5
				except KeyError:
					sum_local+=1.
		dispersion_score = sum_local/float(len(set_u)*len(set_v))
		return dispersion_score

u='HA_2'
print u, distance_matrix[u]

