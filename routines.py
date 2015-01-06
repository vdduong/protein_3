import numpy as np
import re
import math
import os
import random

class Point:
	def __init__(self, name='',x=0.0,y=0.0,z=0.0):
		self.name=name
		self.x=x
		self.y=y
		self.z=z
	
def kabsh(P,Q):
	"""KABSH algorithm to align two point clouds
	return the RMSD between two structures and the aligned structures"""
	C=np.dot(np.transpose(P),Q)
	V,S,W=np.linalg.svd(C)
	d=(np.linalg.det(V)*np.linalg.det(W))<0.0
	if (d):
		S[-1]=-S[-1]
		V[:,-1]=-V[:,-1]
	U=np.dot(V,W)
	P=np.dot(P,U)
	
	return rmsd(P,Q),P

def rmsd(V,W):
	"""return the RMSD from two sets of vectors V and W"""
	D=len(V[0])
	N=len(V)
	rmsd=0.0
	for v,w in zip(V,W):
		rmsd+=sum((v[i]-w[i])**2 for i in range(D))
	return np.sqrt(rmsd/N)

def centroid(X):
	"""return the centroid from a vector set X"""
	C=sum(X)/len(X)
	return C

def fw_algo(vertices, dict_distance):
	"""running the FORD algorithm to generate triangular inequality distance matrix"""
	d = dict(dict_distance)
	for k in vertices:
		for i in vertices:
			for j in vertices:
				d[i][j]=min(d[i][j], d[i][k]+d[k][j])
	return d

def constructing_4_points(list_4_points, dict_distance):
	dict_coord=dict()
	
	atom_1=list_4_points[0]
	x1,y1,z1=0.0,0.0,0.0
	dict_coord[atom_1]=Point(atom_1,x1,y1,z1)
	
	atom_2=list_4_points[1]
	x2=float('%.4f'%dict_distance[atom_1][atom_2])
	y2,z2=0.0,0.0
	dict_coord[atom_2]=Point(atom_2,x2,y2,z2)
	
	atom_3=list_4_points[2]
	d31=dict_distance[atom_3][atom_1]
	d32=dict_distance[atom_3][atom_2]
	x3=float('%.4f'%((d31**2.-d32**2.)/(2.*x2) + x2/2.))
	y3=float('%.4f'%((d31**2.-x3**2.)**0.5))
	z3=0.0
	dict_coord[atom_3]=Point(atom_3,x3,y3,z3)
	
	atom_4=list_4_points[3]
	d41=dict_distance[atom_4][atom_1]
	d42=dict_distance[atom_4][atom_2]
	d43=dict_distance[atom_4][atom_3]
	x4=float('%.4f'%((d41**2.-d42**2.)/(2.*x2) + x2/2.))
	y4=float('%.4f'%((d42**2.-d43**2. - (x4-x2)**2. + (x4-x3)**2.)/(2.*y3) + y3/2.))
	z4=float('%.4f'%((d41**2.-x4**2.-y4**2.)**0.5))
	dict_coord[atom_4]=Point(atom_4,x4,y4,z4)
	
	return dict_coord

def fifth_point(atom_5,list_4_points,dict_distance,dict_coord):
	atom_1,atom_2,atom_3,atom_4 = list_4_points[0],list_4_points[1],list_4_points[2],list_4_points[3]
	x1,y1,z1=dict_coord[atom_1].x,dict_coord[atom_1].y,dict_coord[atom_1].z
	x2,y2,z2=dict_coord[atom_2].x,dict_coord[atom_2].y,dict_coord[atom_2].z
	x3,y3,z3=dict_coord[atom_3].x,dict_coord[atom_3].y,dict_coord[atom_3].z
	x4,y4,z4=dict_coord[atom_4].x,dict_coord[atom_4].y,dict_coord[atom_4].z
	
	d54=dict_distance[atom_5][atom_4]
	d53=dict_distance[atom_5][atom_3]
	d52=dict_distance[atom_5][atom_2]
	d51=dict_distance[atom_5][atom_1]
	
	a=np.array([[2.*(x3-x4),2.*(y3-y4),2.*(z3-z4)],\
				[2.*(x2-x4),2.*(y2-y4),2.*(z2-z4)],\
				[2.*(x1-x4),2.*(y1-y4),2.*(z1-z4)]])
	
	b=np.array([d54**2.-d53**2.+x3**2.+y3**2.+z3**2.-x4**2.-y4**2.-z4**2.,\
				d54**2.-d52**2.+x2**2.+y2**2.+z2**2.-x4**2.-y4**2.-z4**2.,\
				d54**2.-d51**2.+x1**2.+y1**2.+z1**2.-x4**2-y4**2.-z4**2.])
	solve=np.linalg.solve(a,b)
	x5,y5,z5=solve[0],solve[1],solve[2]
	x5=float('%.4f'%x5)
	y5=float('%.4f'%y5)
	z5=float('%.4f'%z5)
	
	return x5,y5,z5

def fifth_point_cramer(atom_5,list_4_points,dict_distance,dict_coord):
	atom_1,atom_2,atom_3,atom_4 = list_4_points[0],list_4_points[1],list_4_points[2],list_4_points[3]
	x1,y1,z1=dict_coord[atom_1].x,dict_coord[atom_1].y,dict_coord[atom_1].z
	x2,y2,z2=dict_coord[atom_2].x,dict_coord[atom_2].y,dict_coord[atom_2].z
	x3,y3,z3=dict_coord[atom_3].x,dict_coord[atom_3].y,dict_coord[atom_3].z
	x4,y4,z4=dict_coord[atom_4].x,dict_coord[atom_4].y,dict_coord[atom_4].z
	
	d54=dict_distance[atom_5][atom_4]
	d53=dict_distance[atom_5][atom_3]
	d52=dict_distance[atom_5][atom_2]
	d51=dict_distance[atom_5][atom_1]
	
	X=np.array([[(d51**2-d52**2)-(x1**2-x2**2)-(y1**2-y2**2)-(z1**2-z2**2),2.*(y2-y1),2.*(z2-z1)],\
				[(d51**2-d53**2)-(x1**2-x3**2)-(y1**2-y3**2)-(z1**2-z3**2),2.*(y3-y1),2.*(z3-z1)],\
				[(d51**2-d54**2)-(x1**2-x4**2)-(y1**2-y4**2)-(z1**2-z4**2),2.*(y4-y1),2.*(z4-z1)]])
	Y=np.array([[2.*(x2-x1),(d51**2-d52**2)-(x1**2-x2**2)-(y1**2-y2**2)-(z1**2-z2**2),2*(z2-z1)],\
				[2.*(x3-x1),(d51**2-d53**2)-(x1**2-x3**2)-(y1**2-y3**2)-(z1**2-z3**2),2.*(z3-z1)],\
				[2.*(x4-x1),(d51**2-d54**2)-(x1**2-x4**2)-(y1**2-y4**2)-(z1**2-z4**2),2.*(z4-z1)]])
	Z=np.array([[2.*(x2-x1),2.*(y2-y1),(d51**2-d52**2)-(x1**2-x2**2)-(y1**2-y2**2)-(z1**2-z2**2)],\
				[2.*(x3-x1),2.*(y3-y1),(d51**2-d53**2)-(x1**2-x3**2)-(y1**2-y3**2)-(z1**2-z3**2)],\
				[2.*(x4-x1),2.*(y4-y1),(d51**2-d54**2)-(x1**2-x4**2)-(y1**2-y4**2)-(z1**2-z4**2)]])
	commun=np.array([[2.*(x2-x1),2.*(y2-y1),2.*(z2-z1)],\
					[2.*(x3-x1),2.*(y3-y1),2.*(z3-z1)],\
					[2.*(x4-x1),2.*(y4-y1),2.*(z4-z1)]])
	if np.linalg.det(commun)!=0.0:
		x5=np.linalg.det(X)/np.linalg.det(commun)
		y5=np.linalg.det(Y)/np.linalg.det(commun)
		z5=np.linalg.det(Z)/np.linalg.det(commun)
		return x5,y5,z5
	else:
		raise 'ZeroDivisionError'

def jacobien(x_,y_,z_,list_4_points,dict_coord,distance_matrix):
	Jac=[]
	atom_1=list_4_points[0]
	x1,y1,z1=dict_coord[atom_1].x,dict_coord[atom_1].y,dict_coord[atom_1].z
	Jac.append(np.array([2.*(x_-x1),2.*(y_-y1),2.*(z_-z1)]))
	
	atom_2=list_4_points[1]
	x2,y2,z2=dict_coord[atom_2].x,dict_coord[atom_2].y,dict_coord[atom_2].z
	Jac.append(np.array([2.*(x_-x2),2.*(y_-y2),2.*(z_-z2)]))

	atom_3=list_4_points[2]
	x3,y3,z3=dict_coord[atom_3].x,dict_coord[atom_3].y,dict_coord[atom_3].z
	Jac.append(np.array([2.*(x_-x3),2.*(y_-y3),2.*(z_-z3)]))
	
	Jac=np.array(Jac)
	return Jac

def iterative_procedure(vertices, distance_matrix):
	list_c=list(vertices)
	list_vertices=list(vertices)
	dict_models=dict()
	nb_model_build=100
	nb_try=0
	while nb_try<nb_model_build:
		try:
			dict_coord=dict()
			random.shuffle(list_c)
			list_4_points=list_c[:4]
			dict_coord=constructing_4_points(list_4_points,distance_matrix)
			for i in xrange(4,len(list_vertices)):
				next_point=list_c[i]
				x_n,y_n,z_n=fifth_point(next_point,list_4_points,distance_matrix,dict_coord)
				dict_coord[next_point]=Point(next_point,x_n,y_n,z_n)
			
			V=[]
			for key in list_vertices:
				x,y,z=dict_coord[key].x,dict_coord[key].y,dict_coord[key].z
				coord=[x,y,z]
				V.append(np.array(coord))
			V=np.array(V)
			V_c=centroid(V)
			V-=V_c
			dict_models[nb_try]=V
			nb_try+=1
		except ValueError: pass
		except ZeroDivisionError: pass
	return dict_models

def breath_first_search(neighbors, start):
	blacknodes=[]
	graynodes=[start]
	while graynodes:
		current=graynodes.pop()
		for neighbor in neighbors[current]:
			if not neighbor in blacknodes+graynodes:
				graynodes.insert(0,neighbor)
		blacknodes.append(current)
	return blacknodes

def connected_components(nodes,distance_matrix):
	nodes=set(nodes)
	max_group=set()
	while nodes:
		start=nodes.pop()
		bfs_group=breath_first_search(distance_matrix,start)
		nodes=nodes.difference(set(bfs_group))
		if len(bfs_group)>len(max_group): max_group=set(bfs_group)
		else:pass
	if len(max_group)>0.9*len(nodes): return max_group
	else:
		print 'Error: network not dense enough'
		return None

def read_restraints(name_file):
	"""return connected component(s) of nodes and the raw (pre-defined) distance matrix"""
	path_local=os.getcwd()
	nmr_restraints=file(path_local+name_file)
	
	for nb_line in xrange(25):
		nmr_restraints.readline()
	
	nodes=set()
	distance_matrix=dict()
	
	for line in nmr_restraints:
		if line!='\n':
			line=line.rstrip('\n\r')
			data_line=re.split(' +',line)
			try:
				int(data_line[0])
				node_1=data_line[2]+'_'+data_line[0]
				node_2=data_line[5]+'_'+data_line[3]
				if (node_1 in nodes) and (node_2 not in nodes):
					distance_matrix[node_1][node_2]=float(data_line[6])
					nodes.add(node_2)
					distance_matrix[node_2]=dict()
					distance_matrix[node_2][node_1]=distance_matrix[node_1][node_2]
				elif (node_2 in nodes) and (node_1 not in nodes):
					distance_matrix[node_2][node_1]=float(data_line[6])
					nodes.add(node_1)
					distance_matrix[node_1]=dict()
					distance_matrix[node_1][node_2]=distance_matrix[node_2][node_1]
				elif (node_1 not in nodes) and (node_2 not in nodes):
					nodes.add(node_1)
					nodes.add(node_2)
					distance_matrix[node_1],distance_matrix[node_2]=dict(),dict()
					distance_matrix[node_1][node_2]=float(data_line[6])
					distance_matrix[node_2][node_1]=float(data_line[6])
				else:
					distance_matrix[node_1][node_2]=float(data_line[6])
					distance_matrix[node_2][node_1]=float(data_line[6])
			except ValueError:
				node_1=data_line[3]+'_'+data_line[1]
				node_2=data_line[6]+'_'+data_line[4]
				if (node_1 in nodes) and (node_2 not in nodes):
					distance_matrix[node_1][node_2]=float(data_line[7])
					nodes.add(node_2)
					distance_matrix[node_2]=dict()
					distance_matrix[node_2][node_1]=float(data_line[7])
				elif (node_2 in nodes) and (node_1 not in nodes):
					distance_matrix[node_2][node_1]=float(data_line[7])
					nodes.add(node_1)
					distance_matrix[node_1]=dict()
					distance_matrix[node_1][node_2]=distance_matrix[node_2][node_1]
				elif (node_2 not in nodes) and (node_1 not in nodes):
					nodes.add(node_1)
					nodes.add(node_2)
					distance_matrix[node_1], distance_matrix[node_2]=dict(),dict()
					distance_matrix[node_1][node_2]=float(data_line[7])
					distance_matrix[node_2][node_1]=float(data_line[7])
				else:
					distance_matrix[node_1][node_2]=float(data_line[7])
					distance_matrix[node_2][node_1]=float(data_line[7])
		else: break
	nodes=connected_components(nodes,distance_matrix) #figure out the biggest connected component
	new_distance_matrix=dict()
	for node_1 in distance_matrix:
		if node_1 in nodes:
			new_distance_matrix[node_1]=dict()
			for node_2 in distance_matrix[node_1]:
				if node_2 in nodes:
					new_distance_matrix[node_1][node_2]=distance_matrix[node_1][node_2]
	return nodes, new_distance_matrix

			
	
