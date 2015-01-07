# test

class Point:
	def __init__(self, name='',x=0.0,y=0.0,z=0.0):
		self.name=name
		self.x=x
		self.y=y
		self.z=z
		
import numpy as np

def raphson(x_,y_,z_,node_,list_4_points,dict_coord,distance_matrix,eps=0.0):
	Jac=[]
	
	atom_1=list_4_points[0]
	x1,y1,z1=dict_coord[atom_1].x,dict_coord[atom_1].y,dict_coord[atom_1].z
	Jac.append(np.array([2.*(x_-x1),2.*(y_-y1),2.*(z_-z1),-2.*(eps)]))
	
	atom_2=list_4_points[1]
	x2,y2,z2=dict_coord[atom_2].x,dict_coord[atom_2].y,dict_coord[atom_2].z
	Jac.append(np.array([2.*(x_-x2),2.*(y_-y2),2.*(z_-z2),-2.*(eps)]))

	atom_3=list_4_points[2]
	x3,y3,z3=dict_coord[atom_3].x,dict_coord[atom_3].y,dict_coord[atom_3].z
	Jac.append(np.array([2.*(x_-x3),2.*(y_-y3),2.*(z_-z3),-2.*(eps+distance_matrix[node_][atom_3])]))
	
	atom_4=list_4_points[3]
	x4,y4,z4=dict_coord[atom_4].x,dict_coord[atom_4].y,dict_coord[atom_4].z
	Jac.append(np.array([2.*(x_-x4),2.*(y_-y4),2.*(z_-z4),-2.*(eps)]))
	
	Jac_inverse=np.linalg.inv(Jac)
	
	value_1=(x_-x1)**2+(y_-y1)**2+(z_-z1)**2-(distance_matrix[node_][atom_1])**2
	value_2=(x_-x2)**2+(y_-y2)**2+(z_-z2)**2-(distance_matrix[node_][atom_2])**2
	value_3=(x_-x3)**2+(y_-y3)**2+(z_-z3)**2-(distance_matrix[node_][atom_3]+eps)**2
	value_4=(x_-x4)**2+(y_-y4)**2+(z_-z4)**2-(distance_matrix[node_][atom_4])**2
	
	value_=[value_1,value_2,value_3,value_4]
	value_=np.array(value_)
	
	return Jac_inverse.dot(value_.T)
	
distance_matrix=dict()
list_4_points=['A','B','C','D']
node='E'
for node_ in ['A','B','C','D','E']:
	distance_matrix[node_]=dict()

dict_coord=dict()
dict_coord['A']=Point('A',0.,0.,0.)
dict_coord['B']=Point('B',3.83,0.,0.)
dict_coord['C']=Point('C',1.9814,2.5382,0.)
dict_coord['D']=Point('D',0.9159,1.4469,2.5601)

distance_matrix['A']['E']=2.44
distance_matrix['B']['E']=5.26
distance_matrix['C']['E']=5.31 #4.26
distance_matrix['D']['E']=3.7
distance_matrix['E']['A']=2.44
distance_matrix['E']['B']=5.26
distance_matrix['E']['C']=5.31 #4.26
distance_matrix['E']['D']=3.7 #2.16

eps=0.
x_,y_,z_=-0.92,0.35,2.23
for i in xrange(500):
	print x_,y_,z_,eps
	appro = raphson(x_,y_,z_,node,list_4_points,dict_coord,distance_matrix,eps)
	print appro
	x_,y_,z_,eps= x_-appro[0],y_-appro[1],z_-appro[2],eps-appro[3]
	dict_coord[node]=Point(node,x_,y_,z_)
	print '_'*10

for node_ in list_4_points:
	print node, node_, distance_matrix[node][node_]


















def small_construction(node, list_points,net_matrix):
	nb_check,missing=check_five(list_points,net_matrix)
	if nb_check>5:
		# 6
		# group de 5 nodes, alternatively call coordinates of these nodes and 
		# choose the best distance corrected
		pass

	elif nb_check==5:
		edge_missed=missing[0]
		print edge_missed
		edge_missed
		
	else:
		pass
		
