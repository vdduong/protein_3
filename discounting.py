
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
import itertools


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

def check_five(list_points,distance_matrix):
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
	eps=0.1
	if nb_check>5:
		for list_4_points in itertools.combinations(list_points,4):
			for node_ in list_4_points:
				list_new=[] # creating simple permutations of 4 nodes
				list_new.append(node_)
				for item in list_4_points:
					if item!=node_:
						list_new.append(item)
				try:
					dict_coord=constructing_4_points(list_new,net_matrix)
					x_,y_,z_=fifth_point_cramer(node,list_new,net_matrix,dict_coord)
					#dict_coord[node]=Point(node,x_,y_,z_)
					for i in xrange(10):
						appro=raphson(x_,y_,z_,node,list_new,dict_coord,net_matrix,eps)
						x_,y_,z_,eps=x_-appro[0],y_-appro[1],z_-appro[2],eps-appro[3]
					if 5.5>=net_matrix[node][list_new[0]]+eps >0.0:
						distance_target = net_matrix[node][list_new[0]]+eps
						distance_initial = net_matrix[node][list_new[0]]
						
						print node, list_new[0], net_matrix[node][list_new[0]],'->','%.2f'%(net_matrix[node][list_new[0]]+eps)
					else: pass
				except ZeroDivisionError: pass
				except ValueError: pass		
	else:
		pass

#--------------------------

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

#list_extract=[]
#while 1:
#	max=0
#	old_len=len(core)
#	for node_core in core:
#		nb_extract=extraction_order(node_core,core,distance_matrix)
#		if max<nb_extract:
#			max=nb_extract
#			node_removed=node_core
#	if old_len==max+1:
#		core.remove(node_removed)
#		list_extract.append(node_removed)
#		print node_removed, len(core)
#	else:
#		break
list_append=['HB3_20', 'HE21_34', 'H_45', 'H_44', 'HB2_60', 'HA_44', 'QG1_41', 'H_53', \
			'HB3_48', 'QG2_49', 'H_68', 'HA_29', 'H_61', 'H_60', 'H_63', 'QD1_68', 'QD1_61',\
			'HA_49', 'HA_20', 'H_37', 'HG2_19', 'HA_26', 'HB3_63', 'HA_24', 'HB2_29', 'HB_68',\
			'HA_21', 'QE_44', 'HA2_27', 'HG2_24', 'HD2_24', 'HA_22', 'HB2_48', 'HA_59', 'HA_55',\
			'HA_51', 'HB2_59', 'QD_44', 'H_21', 'HA_9', 'H_25', 'H_24', 'H_27', 'H_26', 'H_29',\
			'H_28', 'HA_2', 'HA_4', 'QB_1', 'HG_61', 'HB3_14', 'QG_20', 'H_55', 'H_56', 'H_57',\
			'H_50', 'H_52', 'H_33', 'H_59', 'QD1_54', 'HG12_54', 'QG1_31', 'HB3_37', 'HB3_35',\
			'H_31', 'HB3_30', 'QD_4', 'QG2_58', 'HB_54', 'HB2_13', 'QG2_54', 'H_20', 'HA_48',\
			'HA_30', 'HA_13', 'HA_16', 'HA_14', 'HA_19', 'HA_33', 'HA_41', 'HA_34', 'HA_31',\
			'HA_35', 'HE22_34', 'HA_1', 'HD22_48', 'QG2_33', 'QG2_31', 'HB3_59', 'HB2_32',\
			'HB2_35', 'HB2_37', 'HG2_50', 'HA_25', 'H_19', 'H_10', 'H_11', 'HA_32', 'H_14',\
			'H_17', 'QB_26', 'HB_25', 'QG1_6', 'HB3_19', 'H_6', 'H_7', 'H_4', 'H_2', 'H_3',\
			'HB2_50', 'H_9', 'QD_63', 'HB_58', 'HG2_34', 'H_36', 'QB_57', 'H_34', 'H_35',\
			'H_32', 'H_30', 'H_38', 'H_39', 'HG3_34', 'HD2_13', 'QD2_61', 'HD21_48', 'QD1_25',\
			'HA_68', 'HA_61', 'HA_60', 'HA_65', 'HB_49', 'HG12_68']
net_core=set()
for node in list_append:
	net_core.add(node)
	
neighbors=dict()
for node_core in net_core:
	neighbors[node_core]=set()
	for node_ in distance_matrix[node_core]:
		if node_ in net_core and 0.0<distance_matrix[node_core][node_]<float('inf'):
			neighbors[node_core].add(node_)
			
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
nb_=0
for node in net_core:
	list_neighbors=neighbors[node]
	for list_neighbors_ in itertools.combinations(list_neighbors,4):
		nb_check,missing=check_five(list_neighbors_,net_matrix)
		if nb_check==6:
			#print node, list_neighbors_
			#list_neighbors_=list(list_neighbors_)
			#list_neighbors_.append(node)
			#for list_ in itertools.combinations(list_neighbors_,4):
			#	for item in list_neighbors_:
			#		if item in list_: pass
			#		else: node_ext=item
			#	small_construction(node_ext,list_,net_matrix)
			#print '==================================='
			pass
		elif nb_check==4:
			nb_+=1
			#print node, list_neighbors_, missing[0]

print nb_
