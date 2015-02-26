# matching the atom set into the sequence

# dict_assignment
dict_assignment = dict()
for atom in set_atom:
	dict_assignment[atom] = dict()

def proba_matching(cs_exp,cs_theo,tol):
	proba = math.exp(-0.5*(cs_exp - cs_theo)**2./tol**2.)
	return float('%.4f'%proba)

def assignment(set_atom, res, dict_assignment):
	for atom in set_atom:
		for atom_ref in dict_res_ref[res]: # if "H" in atom_ref ?
			if abs(shift_exp[atom] - dict_res_ref[res][atom_ref].cs) <= tol_H:
				dict_assignment[atom][atom_ref] = proba_matching(shift_exp[atom], \
								dict_res_ref[res][atom_ref].cs, tol_H)
	return dict_assignment

def assignment_enhancement(set_atom, res, dict_assignment):
	for atom_1 in set_atom:
		set_assignment_1 = dict_assignment[atom_1].keys()
		for atom_2 in set_atom:
			set_assignment_2 = dict_assignment[atom_2].keys()
			if atom_1!=atom_2:
				try:
					distance_ = distance_matrix[atom_1][atom_2]
					for atom_ref_1 in set_assignment_1:
						for atom_ref_2 in set_assignment_2:
							if distance_ref(atom_ref_1,atom_ref_2) <= 5.5:
								dict_assignment[atom_1][atom_ref_1]+=0.05
								dict_assignment[atom_2][atom_ref_2]+=0.05
				except KeyError:
					pass

# auction algorithm
