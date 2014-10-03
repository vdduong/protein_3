
# read NMR restraints from mr.txt file
import os
import re
import math
from routines import *
from spring_modif import *
import matplotlib.pyplot as plt
import pylab
from mpl_toolkits.mplot3d import Axes3D
import numpy as np

path_local = os.getcwd()

ref_model = dict()
nodes_ref_model = []
file = open('/Users/james/Documents/protein/H_coord_1cpz')
for line in file:
	line = line.rstrip('\n\r')
	line = re.split(' +', line)
	ref_model[line[1]] = Point(line[1], float(line[2]), float(line[3]), float(line[4]))
	nodes_ref_model.append(line[1])

nmr_restraints = open(path_local + '/1CPZ.new_restraints.txt')

nodes = set()
distance_matrix = {}
initial_restraints = {}

for line in nmr_restraints:
	if line!='\n':
		line = line.rstrip('\n\r')
		data_line = re.split(' +', line)
		node_1 = data_line[0]
		node_2 = data_line[1]
		if (node_1 in nodes) and (node_2 not in nodes):
			distance_matrix[node_1][node_2] = float(data_line[2])
			initial_restraints[node_1][node_2] = float(data_line[2])
			nodes.add(node_2)
			distance_matrix[node_2] = {}
			distance_matrix[node_2][node_1] = distance_matrix[node_1][node_2]
			
			initial_restraints[node_2] = {}
			initial_restraints[node_2][node_1] = initial_restraints[node_1][node_2]
			
		elif (node_2 in nodes) and (node_1 not in nodes):
			distance_matrix[node_2][node_1] = float(data_line[2])
			initial_restraints[node_2][node_1] = float(data_line[2])
			nodes.add(node_1)
			distance_matrix[node_1] = {}
			distance_matrix[node_1][node_2] = distance_matrix[node_2][node_1]
			initial_restraints[node_1] = {}
			initial_restraints[node_1][node_2] = initial_restraints[node_2][node_1]
			
		elif (node_1 not in nodes) and (node_2 not in nodes):
			nodes.add(node_1)
			nodes.add(node_2)
			distance_matrix[node_1], distance_matrix[node_2] = {},{}
			distance_matrix[node_1][node_2] = float(data_line[2])
			distance_matrix[node_2][node_1] = float(data_line[2])
			
			initial_restraints[node_1], initial_restraints[node_2] = {},{}
			initial_restraints[node_1][node_2] = float(data_line[2])
			initial_restraints[node_2][node_1] = float(data_line[2])
		else:
			distance_matrix[node_1][node_2] = float(data_line[2])
			distance_matrix[node_2][node_1] = float(data_line[2])
			initial_restraints[node_1][node_2] = float(data_line[2])
			initial_restraints[node_2][node_1] = float(data_line[2])
	else:
		break

nodes = connected_components(nodes, distance_matrix)
print len(nodes), ' nodes'

nodes_del = set(distance_matrix.keys()).difference(nodes)

nodes_del_1 = set()
for node in nodes:
	if 13<= int(re.split('_',node)[1]) <=24:
		pass
	else:
		nodes_del_1.add(node)
		nodes_del.add(node)

for node in nodes_del_1:
	nodes.remove(node)
		
for node_1 in nodes_del:
		del distance_matrix[node_1]
		del initial_restraints[node_1]

for node in nodes:
	distance_matrix[node][node] = 0.0
	initial_restraints[node][node]=0.0

distance_to_be_del = []
for node_1 in distance_matrix:
	for node_2 in distance_matrix[node_1]:
		if node_2 not in nodes:
			distance_to_be_del.append([node_1, node_2])

for item in distance_to_be_del:
	del distance_matrix[item[0]][item[1]]

for node_1 in nodes:
	for node_2 in nodes:
		try:
			distance_ = distance_matrix[node_1][node_2]
		except KeyError:
			distance_matrix[node_1][node_2] = float('inf')
			distance_matrix[node_2][node_1] = float('inf')


for node_1 in nodes:
	for node_2 in nodes:
		try:
			distance_ = initial_restraints[node_1][node_2]
		except KeyError:
			initial_restraints[node_1][node_2] = float('inf')
			initial_restraints[node_2][node_1] = float('inf')

distance_matrix = fw_algo(nodes, distance_matrix)	

#for node_1 in distance_matrix:
#	for node_2 in distance_matrix[node_1]:
#		distance_ = (ref_model[node_1].x - ref_model[node_2].x)**2 +\
#			(ref_model[node_1].y - ref_model[node_2].y)**2 + (ref_model[node_1].z - ref_model[node_2].z)**2
		#if abs(distance_matrix[node_1][node_2] - math.sqrt(distance_)) > 0.:
		#	distance_corr = distance_matrix[node_1][node_2]*0.65
		#	if abs(distance_corr - math.sqrt(distance_)) > 1.:
		#		print node_1, node_2, distance_matrix[node_1][node_2], math.sqrt(distance_)
		#		nb_c +=1
		#distance_matrix[node_1][node_2] = math.sqrt(distance_)


dict_models = iterative_procedure(nodes, distance_matrix)

list_nodes = list(nodes)


######### 
## Average model computation begins here



average_model = dict()
for node in list_nodes:
	average_model[node] = Point(node)
nb_average = 0.0

list_rmsd, list_harmonic = [],[]

for key_model in dict_models.keys():
	model = dict_models[key_model]
	d = distance_matrix_generation_from_array(nodes, model)
	sum_ = 0.0
	time_ = 0.0
	for node_1 in list(nodes):
		for node_2 in list(nodes):
			sum_ += (distance_matrix[node_1][node_2]-d[node_1][node_2])**2
			time_ +=1.

	harmonic_score = sum_/time_
	rmsd_ = rmsd_individual(list_nodes, model)
	if harmonic_score < 60:
		list_rmsd.append(rmsd_)
		list_harmonic.append(harmonic_score)	
	#print harmonic_score
	if harmonic_score <= 20.:
		#print "Model ", key_model
		
		#print 'model score %.4f'%(harmonic_score)	
		if nb_average == 0.0: # first model to be added ?
			for index_node, node in enumerate(list_nodes):
				average_model[node].x += model[index_node][0]
				average_model[node].y += model[index_node][1]
				average_model[node].z += model[index_node][2]
			nb_average+=1.
			model_1 = model # assign model_1

		else:
			model_inverse = []
			for index_node, node in enumerate(list_nodes):
				coord_inverse = [model[index_node][0], -model[index_node][1], model[index_node][2]]
				model_inverse.append(np.array(coord_inverse))
			model_inverse = np.array(model_inverse)
			rmsd_k, model = kabsh(model, model_1)
			rmsd_k_inverse, model_inverse = kabsh(model_inverse, model_1)
			if rmsd_k == min(rmsd_k, rmsd_k_inverse):
				for index_node, node in enumerate(list_nodes):
					average_model[node].x += model[index_node][0]
					average_model[node].y += model[index_node][1]
					average_model[node].z += model[index_node][2]
				nb_average +=1.0
			else:
				for index_node, node in enumerate(list_nodes):
					average_model[node].x += model_inverse[index_node][0]
					average_model[node].y += model_inverse[index_node][1]
					average_model[node].z += model_inverse[index_node][2]
				nb_average +=1.0


for node in list_nodes:
	average_model[node].x = average_model[node].x / nb_average
	average_model[node].y = average_model[node].y / nb_average
	average_model[node].z = average_model[node].z / nb_average

print len(average_model), 'nodes'

#fig=plt.figure()
#ax_view=fig.add_subplot(1,1,1)
#ax_view.scatter([rmsd for rmsd in list_rmsd], \
#	[harmonic for harmonic in list_harmonic], color='green',s=30.)

dict_x,dict_y,dict_z = dict(),dict(),dict()
list_bb = []
P,P_inverse,Q = [],[],[]
for index_node, node in enumerate(list_nodes):
	if re.split('_', node)[0]=='H':
		list_bb.append(node)
		dict_x[node] = average_model[node].x
		dict_y[node] = average_model[node].y
		dict_z[node] = average_model[node].z
		coord_exp = [dict_x[node], dict_y[node], dict_z[node]]
		coord_exp_inverse = [dict_x[node], -dict_y[node], dict_z[node]]
		P.append(coord_exp)
		P_inverse.append(coord_exp_inverse)
		coord_ref = [ref_model[node].x, ref_model[node].y, ref_model[node].z]
		Q.append(coord_ref)


P = np.array(P)
P_inverse = np.array(P_inverse)
Q = np.array(Q)

P_c =centroid(P)
P_inverse_c = centroid(P_inverse)
Q_c = centroid(Q)
P-=P_c
P_inverse-=P_inverse_c
Q-=Q_c

rmsd_k, P = kabsh(P,Q)
rmsd_k_inverse,P_inverse = kabsh(P_inverse,Q)
inversion = False # to know if the structure should be reversed or not ?

dict_x, dict_y, dict_z = dict(), dict(), dict()

if rmsd_k < rmsd_k_inverse:		
	print rmsd_k
	for index_node, node in enumerate(list_bb):
		dict_x[node] = P[index_node][0]
		dict_y[node] = P[index_node][1]
		dict_z[node] = P[index_node][2]
	rmsd_ = rmsd_k
else:
	print rmsd_k_inverse
	inversion = True
	for index_node, node in enumerate(list_bb):
		dict_x[node] = P_inverse[index_node][0]
		dict_y[node] = P_inverse[index_node][1]
		dict_z[node] = P_inverse[index_node][2]
	rmsd_ = rmsd_k_inverse
list_ordering = []
for node in list_bb:
	ind = int(re.split('_', node)[1])
	list_ordering.append(ind)
list_ordering.sort()
for index_node, node in enumerate(list_ordering):
	list_ordering[index_node] = 'H_' + str(node)
		
X,Y,Z = [],[],[]
X_ref,Y_ref,Z_ref = [],[],[]
		
for node in list_ordering:
	X.append(float(dict_x[node]))
	Y.append(float(dict_y[node]))
	Z.append(float(dict_z[node]))
	X_ref.append(ref_model[node].x)
	Y_ref.append(ref_model[node].y)
	Z_ref.append(ref_model[node].z)
		
fig = plt.figure()		
ax = Axes3D(fig)
ax.plot(X,Y,Z,c='r',marker='o')
ax.plot(X_ref,Y_ref,Z_ref,c='b', marker='o')
plt.show()

def distance_between_points(proton_a, proton_b, dict_coord):
	distance = (dict_coord[proton_a].x - dict_coord[proton_b].x)**2 + \
		(dict_coord[proton_a].y - dict_coord[proton_b].y)**2 + \
		(dict_coord[proton_a].z - dict_coord[proton_b].z)**2
	return math.sqrt(distance)

def mobility_radius(node, dict_coord, initial_restraints):
	D_0, D_1 = [], []
	nodes = dict_coord.keys()
	for node_ in nodes:
		if node_!=node and initial_restraints[node][node_]<float('inf') and \
			abs(distance_between_points(node_, node, dict_coord) - initial_restraints[node_][node])<=1.:
			D_1.append(abs(initial_restraints[node][node_] - distance_between_points(node_, node, dict_coord)))
	try:
		d_1 = min(D_1)
		r_node = d_1
		return r_node
	except ValueError:
		return None
		
def move(node, dict_coord, initial_restraints):
	r_node = mobility_radius(node, dict_coord, initial_restraints)
	F_x, F_y, F_z = 0.0, 0.0, 0.0
	for node_ in nodes:
		if node_!=node and initial_restraints[node][node_]<float('inf'):
			x_move = dict_coord[node].x - dict_coord[node_].x
			y_move = dict_coord[node].y - dict_coord[node_].y
			z_move = dict_coord[node].z - dict_coord[node_].z
			norm_ = math.sqrt(x_move**2 + y_move**2 + z_move**2)
			F_x +=x_move/norm_
			F_y +=y_move/norm_
			F_z +=z_move/norm_
	norm_total = math.sqrt(F_x**2 + F_y**2 + F_z**2)
	dict_coord[node].x-=F_x*r_node/norm_total
	dict_coord[node].y-=F_y*r_node/norm_total
	dict_coord[node].z-=F_z*r_node/norm_total
	return dict_coord		

def move_new(node, dict_coord, initial_restraints):
	list_nodes_iteration = []
	for node_ in nodes:
		if node_!=node and initial_restraints[node][node_]<float('inf'):
			list_nodes_iteration.append(node_)
	for index_1 in xrange(0, len(list_nodes_iteration)-1):
		for index_2 in xrange(index_1+1, len(list_nodes_iteration)):
			node_1 = list_nodes_iteration[index_1]
			node_2 = list_nodes_iteration[index_2]
			distance_1 = distance_between_points(node_1, node, dict_coord)
			distance_2 = distance_between_points(node_2, node, dict_coord)
			eps_1 = abs(distance_1 - initial_restraints[node][node_1])
			eps_2 = abs(distance_2 - initial_restraints[node][node_2])
			#print eps_1, eps_2
			if eps_1>0.5 and eps_2>0.5:
				x_eps_1 = initial_restraints[node][node_1]*(dict_coord[node].x-dict_coord[node_1].x)/distance_1 + dict_coord[node_1].x
				y_eps_1 = initial_restraints[node][node_1]*(dict_coord[node].y-dict_coord[node_1].y)/distance_1 + dict_coord[node_1].y
				z_eps_1 = initial_restraints[node][node_1]*(dict_coord[node].z-dict_coord[node_1].z)/distance_1 + dict_coord[node_1].z
				
				x_eps_2 = initial_restraints[node][node_2]*(dict_coord[node].x-dict_coord[node_2].x)/distance_2 + dict_coord[node_2].x
				y_eps_2 = initial_restraints[node][node_2]*(dict_coord[node].y-dict_coord[node_2].y)/distance_2 + dict_coord[node_2].y
				z_eps_2 = initial_restraints[node][node_2]*(dict_coord[node].z-dict_coord[node_2].z)/distance_2 + dict_coord[node_2].z
			

				a = np.array([[x_eps_2-x_eps_1, y_eps_2-y_eps_1, z_eps_2-z_eps_1],\
						[y_eps_2-y_eps_1, -(x_eps_2-x_eps_1),0.0],\
						[0.0, -(z_eps_2-z_eps_1),y_eps_2-y_eps_1]])
				b = np.array([dict_coord[node].x*(x_eps_2-x_eps_1)+dict_coord[node].y*(y_eps_2-y_eps_1)+dict_coord[node].z*(z_eps_2-z_eps_1),\
						x_eps_1*(y_eps_2-y_eps_1)-y_eps_1*(x_eps_2-x_eps_1),\
						z_eps_1*(y_eps_2-y_eps_1)-y_eps_1*(z_eps_2-z_eps_1)])
			
				try:
					solve = np.linalg.solve(a,b)
					x_h, y_h, z_h = solve[0], solve[1], solve[2]
					
					x_ = 2.*(x_h - dict_coord[node].x) + dict_coord[node].x
					y_ = 2.*(y_h - dict_coord[node].y) + dict_coord[node].y
					z_ = 2.*(z_h - dict_coord[node].z) + dict_coord[node].z
								
					dict_coord[node] = Point(node, x_, y_, z_)
				except np.linalg.linalg.LinAlgError as err:
					if 'Singular matrix' in err.message:
						pass
					else: raise
			else:
				pass			
	return dict_coord

for i in xrange(150):
#while rmsd_ >=2.:		
	for node in average_model:
		#r_node = mobility_radius(node, average_model, initial_restraints)
		#print node, r_node
		#if r_node!=None:
			#average_model = move(node, average_model, initial_restraints)
		average_model = move_new(node, average_model, initial_restraints)
	
	dict_x,dict_y,dict_z = dict(),dict(),dict()
	list_bb = []
	P,P_inverse,Q = [],[],[]
	for index_node, node in enumerate(list_nodes):
		if re.split('_', node)[0]=='H':
			list_bb.append(node)
			dict_x[node] = average_model[node].x
			dict_y[node] = average_model[node].y
			dict_z[node] = average_model[node].z
			coord_exp = [dict_x[node], dict_y[node], dict_z[node]]
			coord_exp_inverse = [dict_x[node], -dict_y[node], dict_z[node]]
			P.append(coord_exp)
			P_inverse.append(coord_exp_inverse)
			coord_ref = [ref_model[node].x, ref_model[node].y, ref_model[node].z]
			Q.append(coord_ref)
	P = np.array(P)
	P_inverse = np.array(P_inverse)
	Q = np.array(Q)
	
	P_c =centroid(P)
	P_inverse_c = centroid(P_inverse)
	Q_c = centroid(Q)
	P-=P_c
	P_inverse-=P_inverse_c
	Q-=Q_c
	
	rmsd_k, P = kabsh(P,Q)
	rmsd_k_inverse,P_inverse = kabsh(P_inverse,Q)
	inversion = False # to know if the structure should be reversed or not ?
	dict_x,dict_y,dict_z = dict(),dict(),dict()
	if rmsd_k < rmsd_k_inverse:		
		print rmsd_k
		for index_node, node in enumerate(list_bb):
			dict_x[node] = P[index_node][0]
			dict_y[node] = P[index_node][1]
			dict_z[node] = P[index_node][2]
		rmsd_ = rmsd_k
	else:
		print rmsd_k_inverse
		inversion = True
		for index_node, node in enumerate(list_bb):
			dict_x[node] = P_inverse[index_node][0]
			dict_y[node] = P_inverse[index_node][1]
			dict_z[node] = P_inverse[index_node][2]
		rmsd_ = rmsd_k_inverse
	
#d = distance_matrix_generation_from_coord(nodes, average_model)
#nodes = average_model.keys()
#for node in nodes:
#	for node_ in nodes:
#		if node!=node_ and initial_restraints[node][node_]<float('inf'):
#			print node, node_, d[node][node_], initial_restraints[node][node_]
	
	#list_ordering = []
	#for node in list_bb:
	#	ind = int(re.split('_', node)[1])
	#	list_ordering.append(ind)
	#list_ordering.sort()
	#for index_node, node in enumerate(list_ordering):
	#	list_ordering[index_node] = 'H_' + str(node)
		
	#X,Y,Z = [],[],[]
	#X_ref,Y_ref,Z_ref = [],[],[]
		
	#for node in list_ordering:
	#	X.append(float(dict_x[node]))
	#	Y.append(float(dict_y[node]))
	#	Z.append(float(dict_z[node]))
	#	X_ref.append(ref_model[node].x)
	#	Y_ref.append(ref_model[node].y)
	#	Z_ref.append(ref_model[node].z)
	#fig = plt.figure()		
	#ax = Axes3D(fig)
	#ax.plot(X,Y,Z,c='r',marker='o')
	#ax.plot(X_ref,Y_ref,Z_ref,c='b', marker='o')
	#plt.show()

