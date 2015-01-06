import os
import re
import math
import random
import itertools
from routines import *
import matplotlib.pyplot as plt
import pylab
from mpl_toolkits.mplot3d import Axes3D
import numpy as np

name_file='/1CPZ.mr.txt'
nodes, distance_matrix=read_restraints(name_file)
print len(nodes), 'nodes'

for node in nodes:
	distance_matrix[node][node]=0.0

core=set()
core=set(nodes)

for i in xrange(20):
	nodes_removed=set()
	for node in core:
		nb_connect=0
		for node_ in distance_matrix[node].keys():
			if node_ in core and float('inf')>distance_matrix[node][node_]>0.0:
				nb_connect+=1
		if nb_connect<4:
			nodes_removed.add(node)
	for node_ in nodes_removed:
		try:
			core.remove(node_)
		except KeyError:
			pass
	print len(core)

neighbors=dict()
for node_core in core:
	neighbors[node_core]=set()
	for node_ in distance_matrix[node_core]:
		if node_ in core and 0.0<distance_matrix[node_core][node_]<float('inf'):
			neighbors[node_core].add(node_)

def extraction_order(node_core,core,distance_matrix):
	new_core=set(core)
	new_core.remove(node_core)
	old_score=len(new_core)+1
	new_score=len(new_core)
	while old_score!=new_score:
		nodes_removed=set()
		old_score=len(new_core)
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
				new_score=len(new_core)
			except KeyError:
				pass
	return len(new_core)

list_extract=[]
while 1:
	max=0
	old_len=len(core)
	for node_core in core:
		nb_extract=extraction_order(node_core,core,distance_matrix)
		if max<nb_extract:
			max=nb_extract
			node_removed=node_core
	if old_len==max+1:
		core.remove(node_removed)
		list_extract.append(node_removed)
		print node_removed, len(core)
	else:
		break
print len(core)
#print list_extract
# net_core

net_core=set(core)
net_matrix=dict()
for node_core in net_core:
	net_matrix[node_core]=dict()
for node_1 in net_core:
	for node_2 in net_core:
		try:
			net_matrix[node_1][node_2]=distance_matrix[node_1][node_2]
		except KeyError:
			if node_1==node_2:
				net_matrix[node_1][node_2]=0.0
			else:
				net_matrix[node_1][node_2]=float('inf')
#net_matrix=fw_algo(net_core,net_matrix)
#dict_models=iterative_procedure(net_core,net_matrix)

def check_five(list_points, distance_matrix):
	nb_edges=0
	missing=[]
	for node_1 in list_points:
		for node_2 in list_points:
			if 0.0<distance_matrix[node_1][node_2]<float('inf'):
				nb_edges+=1
			elif distance_matrix[node_1][node_2]==float('inf'):
				missing.append([node_1,node_2])	
	return nb_edges/2, missing

def small_construction(node, list_points,net_matrix):
	nb_check,missing=check_five(list_points,net_matrix)
	if nb_check>5:
		dict_coord=constructing_4_points(list(list_points),net_matrix)
		x_,y_,z_=fifth_point_cramer(node,list(list_points),net_matrix,dict_coord)
		dict_coord[node]=Point(node,x_,y_,z_)
		sum=0.0
		nb=0.0
		for node_ in list_points:
			print node_, dict_coord[node_].x, dict_coord[node_].y, dict_coord[node_].z
		print list(list_points)+[node]
		for node_1 in list(list_points)+[node]:
			for node_2 in list(list_points)+[node]:
				distance_=math.sqrt((dict_coord[node_1].x-dict_coord[node_2].x)**2+ \
					(dict_coord[node_1].y-dict_coord[node_2].y)**2+\
					(dict_coord[node_1].z-dict_coord[node_2].z)**2)
				#sum+=(distance-net_matrix[node_1][node_2])**2
				#nb+=1.
				print node_1,node_2,'%.4f'%distance_, net_matrix[node_1][node_2]
	else:
		pass
		

for node in net_core:
	list_points = neighbors[node]
	#print node
	nb_ = 0
	for list_4_points in itertools.combinations(list_points,4):
		#print list_4_points
		nb_+=1
	#print node, nb_
	if nb_==1: 
		#print node, neighbors[node]
		nb_check,missing=check_five(list_points,net_matrix)
		#print nb_check
		if nb_check==6:
			small_construction(node,list_points,net_matrix)
		#for node_1 in list_points:
		#	for node_2 in list_points:
		#		try:
		#			if 0.0<net_matrix[node_1][node_2]<float('inf'):
		#				print node_1, node_2, net_matrix[node_1][node_2]
		#		except KeyError: pass

