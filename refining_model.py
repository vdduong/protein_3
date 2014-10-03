# refining the model issued from initial computation

# given the dict_coord containing the coordinates of the configuration
# given the initial distance restraints from NMR

import math

def distance_between_points(proton_a, proton_b, dict_coord):
	distance = (dict_coord[proton_a].x - dict_coord[proton_b].x)**2 + \
		(dict_coord[proton_a].y - dict_coord[proton_b].y)**2 + \
		(dict_coord[proton_a].z - dict_coord[proton_b].z)**2
	return math.sqrt(distance)

def well_placed_pair(proton_a, proton_b, dict_coord, initial_restraints, threshold=6.0):
	'''return if proton b is well-placed in respect to proton a'''
	#threshold = 6.0 # angstroms
	distance_ = distance_between_points(proton_a, proton_b, dict_coord)
	if initial_restraints[proton_a][proton_b]<float('inf') and distance_<= threshold:
		return True
	elif initial_restraints[proton_a][proton_b] == float('inf') and distance_ > threshold:
		return True
	else:
		return False

def well_placed(proton_a, dict_coord, initial_restraints, threshold=6.0):
	if 1:
		return True
	else:
		return False

def move(node, dict_coord, initial_restraints, threshold=6.0):
	r_node = # radius of mobility at threshold t of node i
	F_x, F_y, F_z = 0.0, 0.0, 0.0
	nodes = dict_coord.keys()
	for node_ in nodes:
		if not well_placed(node_, node, dict_coord, initial_restraints, threshold=6.0):
			if initial_restraints[node][node_]<float('inf'):
				F_x = F_x - distance_between_points(node_, node, dict_coord)/initial_restraints[node][node_]
				F_y = F_y - distance_between_points(node_, node, dict_coord)/initial_restraints[node][node_]
				F_z = F_z - distance_between_points(node_, node, dict_coord)/initial_restraints[node][node_]
			else:
				F_x = F_x + distance_between_points(node_, node, dict_coord)/initial_restraints[node][node_]
				F_y = F_y + distance_between_points(node_, node, dict_coord)/initial_restraints[node][node_]
				F_z = F_z + distance_between_points(node_, node, dict_coord)/initial_restraints[node][node_]
	x_ = dict_coord[node].x + F_x
	y_ = dict_coord[node].y + F_y
	z_ = dict_coord[node].z + F_z
	return x_, y_, z_
				
def correct(dict_coord, initial_restraints, threshold=6.0):
	nodes = dict_coord.keys()
	for node in nodes:
		if not well_placed(node, dict_coord, initial_restraints, threshold=6.0):
			x_, y_, z_ = move(node, dict_coord, initial_restraints, threshold=6.0)
			dict_coord[node] = Point(node, x_, y_, z_)
		else:
			pass
	return dict_coord


def mobility_radius(node, dict_coord, initial_restraints, threshold=6.0):
	D_0, D_1 = [], []
	nodes = dict_coord.keys()
	for node_ in nodes:
		if node_!=node and initial_restraints[node][node_]<float('inf') and \
			distance_between_points(node_, node, dict_coord)<=threshold:
			D_1.append(distance_between_points(node_, node, dict_coord))
		elif node_!=node and initial_restraints[node][node]==float('inf') and \
			distance_between_points(node_, node, dict_coord)>threshold:
			D_0.append(distance_between_points(node_, node, dict_coord))
	d_0 = min(D_0)
	d_1 = max(D_1)
	r_node = min(d_0-threshold, threshold-d_1)
	return r_node

def move_new(node, dict_coord, initial_restraints):
	list_nodes_iteration = []
	for node_ in nodes:
		if node_!=node and initial_restraints[node][node_]<float('inf'):
			list_nodes_iteration.append(node_)
	for index_1 in xrange(0, len(list_nodes_iteration)-1):
		for index_2 in xrange(index_1, len(list_nodes_iteration)):
			node_1 = list_nodes_iteration[index_1]
			node_2 = list_nodes_iteration[index_2]
			distance_1 = distance_between_points(node_1, node, dict_coord)
			distance_2 = distance_between_points(node_2, node, dict_coord)
			eps_1 = abs(distance_1 - initial_restraints[node][node_1])
			eps_2 = abs(distance_2 - initial_restraints[node][node_2])
			
			x_eps_1 = initial_restraints[node][node_1]*(dict_coord[node].x-dict_coord[node_1].x)/distance_1 + dict_coord[node_1].x
			y_eps_1 = initial_restraints[node][node_1]*(dict_coord[node].y-dict_coord[node_1].y)/distance_1 + dict_coord[node_1].y
			z_eps_1 = initial_restraints[node][node_1]*(dict_coord[node].z-dict_coord[node_1].z)/distance_1 + dict_coord[node_1].z
			
			x_eps_2 = initial_restraints[node][node_2]*(dict_coord[node].x-dict_coord[node_2].x)/distance_1 + dict_coord[node_2].x
			y_eps_2 = initial_restraints[node][node_2]*(dict_coord[node].y-dict_coord[node_2].y)/distance_1 + dict_coord[node_2].y
			z_eps_2 = initial_restraints[node][node_2]*(dict_coord[node].z-dict_coord[node_2].z)/distance_1 + dict_coord[node_2].z
			
			a = np.array([[x_eps_2-x_eps_1, y_eps_2-y_eps_1, z_eps_2-z_eps_1],\
						[y_eps_2-y_eps_1, -(x_eps_2-x_eps_1),0.0],\
						[0.0, -(z_eps_2-z_eps_1),y_eps_2-y_eps_1]])
			b = np.array([dict_coord[node].x*(x_eps_2-x_eps_1)+dict_coord[node].y*(y_eps_2-y_eps_1)+dict_coord[node].z*(z_eps_2-z_eps_1),\
						x_eps_1*(y_eps_2-y_eps_1)-y_eps_1*(x_eps_2-x_eps_1),\
						z_eps_1*(y_eps_2-y_eps_1)-y_eps_1*(z_eps_2-z_eps_1)])
			
			try:
				solve = np.linalg.solve(a,b)
				x_, y_, z_ = solve[0], solve[1], solve[2]
				dict_coord[node] = Point(node, x_, y_, z_)
			except np.linalg.linalg.LinAlgError as err:
				if 'Singular matrix' in err.message:
					pass
				else: raise
					
	return dict_coord





