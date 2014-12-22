
# read NMR restraints from mr.txt file
import os
import re
import math
from routines import *
#from spring_modif import *
import matplotlib.pyplot as plt
import pylab
from mpl_toolkits.mplot3d import Axes3D
import numpy as np


def d_v(s,t,u,v,distance_matrix):
	'''distance function on the nodes of C_uv'''
	try:
		distance_st = distance_matrix[s][t]
		return 0.5
	except KeyError:
		set_s = set(distance_matrix[s].keys())
		set_t = set(distance_matrix[t].keys())
		set_intersect = set_s.intersection(set_t)
		if len(set_intersect) > 2: return 0
		else: return 1

def dispersion(u,v,distance_matrix):
	disp = 0
	set_u = set(distance_matrix[u].keys())
	set_v = set(distance_matrix[v].keys())
	set_commun = set_u.intersection(set_v)
	for s in set_commun:
		for t in set_commun:
			if s!=t:
				disp+=d_v(s,t,u,v,distance_matrix)
	if len(set_commun) > 0:
		return disp/float(len(set_commun))
	else:
		return 0.0

def rec_disp(u,emb,distance_matrix,disp):
	neighbors = []
	set_commun = dict()
	set_u = set(distance_matrix[u].keys())
	for node in distance_matrix[u].keys():
		if node !=u: 
			neighbors.append(node)
			set_node = set(distance_matrix[node].keys())
			set_unode = set_u.intersection(set_node)
			set_commun[node] = set_unode
	
	for node in neighbors: disp[u][node] = 1.0	# reinitialising the dispersion dictionary
	dict_value = dict()
	
	for nb_iteration in xrange(3):
		for v in neighbors:
			sum_nominator = 0.0
			set_uv = set_commun[v]
			for s in set_uv:
				for t in set_uv:
					if s!=t:
						sum_nominator+=2.*d_v(s,t,u,v,distance_matrix)
			for w in neighbors:
				sum_nominator+=disp[u][w]**2.
			try:
				score_v = sum_nominator/float(emb[u][v])
				dict_value[v] = score_v
			except ZeroDivisionError:
				dict_value[v] = 0.0
		for v in neighbors:
			disp[u][v] = dict_value[v]
	max = 0.0
	for v in disp[u].keys():
		if disp[u][v] > max: max = disp[u][v]
	if max == 0.0: 
		pass
	else:
		for v in disp[u].keys():
			disp[u][v] = disp[u][v]/max
	return disp

path_local = os.getcwd()
name_file = '/1CPZ.mr.txt'
#name_file = '/2KA0.mr.txt'
#name_file = '/2MU2.mr.txt'

nodes, distance_matrix = read_restraints(name_file)

print len(nodes), ' nodes'

for node in nodes:
	distance_matrix[node][node] = 0.0

#nb_add = 0		
#nb_correct = 0
core = set()
core = set(nodes)

for i in xrange(20):				
	nodes_removed = set()
	for node in core:
		nb_connect = 0
		for node_ in distance_matrix[node].keys():
			if node_ in core and distance_matrix[node][node_] <  float('inf') and distance_matrix[node][node_] > 0.0:
				nb_connect +=1
		if nb_connect < 4:
			nodes_removed.add(node)
	for node_ in nodes_removed:
		try:
			core.remove(node_)
		except KeyError:
			pass
	print i, len(core)


#for node in core:
#	for node_ in distance_matrix[node].keys():
#		if node_ in core and distance_matrix[node][node_] < float('inf') and distance_matrix[node][node_]>0.0:
#			print node, node_, distance_matrix[node][node_]
#	print '_'*10

print "NEW TEST: removing points from core set"

def extraction_order(node_core, core, distance_matrix):
	'''return extraction order'''
	new_core = set(core)
	new_core.remove(node_core)
	old_score = len(new_core)+1
	new_score = len(new_core)
	while old_score!=new_score:
		nodes_removed = set()
		old_score = len(new_core)
		for node in new_core:
			nb_connect=0
			for node_ in distance_matrix[node].keys():
				if node_ in new_core and 0.0<distance_matrix[node][node_]<float('inf'):
					nb_connect+=1
			if nb_connect<4:
				nodes_removed.add(node)
		for node_ in nodes_removed:
			try:
				new_core.remove(node_)
				new_score = len(new_core)
			except KeyError:
				pass
	return len(new_core)

construction_next = list()

for i in xrange(266):			
	max = 0
	for node_core in core:
		nb_extract = extraction_order(node_core, core, distance_matrix)
		if max<nb_extract:
			max=nb_extract
			node_removed = node_core
	if len(core) == max+1:
		core.remove(node_removed)
		print node_removed, max
	else:
		pass
