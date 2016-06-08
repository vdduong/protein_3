import numpy as np
import scipy.stats
import random
import math
from routines import *

def T(x, T0, T1, k=1.):
	# apply an affine transformation to x

	y = x*T0.T
	y[:,0]+=T1[0,0]
	y[:,1]+=T1[1,0]
	y[:,2]+=T1[2,0]
	y[:,3]+=T1[3,0]
	return y*k

def random_unit_vector():
	v0 = np.random.random(3)
	length_ = math.sqrt(sum(v0[i]**2 for i in xrange(3)))	
	for i in xrange(3):
		v0[i]/=length_
	return v0

def rotation_matrix_cp(axis, theta):
    """
    Return the rotation matrix associated with counterclockwise rotation about
    the given axis by theta radians.
    """
    axis = np.asarray(axis)
    theta = np.asarray(theta)
    axis = axis/math.sqrt(np.dot(axis, axis))
    a = math.cos(theta/2.0)
    b, c, d = -axis*math.sin(theta/2.0)
    aa, bb, cc, dd = a*a, b*b, c*c, d*d
    bc, ad, ac, ab, bd, cd = b*c, a*d, a*c, a*b, b*d, c*d
    return np.matrix([[aa+bb-cc-dd, 2*(bc+ad), 2*(bd-ac), 0.],
                     [2*(bc-ad), aa+cc-bb-dd, 2*(cd+ab), 0.],
                     [2*(bd+ac), 2*(cd-ab), aa+dd-bb-cc, 0.],
                     [0., 0., 0., 1.]])

def translate(X, Y, errorfct):
	mx, my = [], []
	for i in xrange(3):
		mx.append([X[:,i].mean()])
		my.append([Y[:,i].mean()])
		
	mx.append([0.])
	my.append([0.])	# chemical shifts are not translatable ?
	
	mx = np.matrix(mx)
	my = np.matrix(my)
	translation = mx-my
	I = np.matrix(np.identity(4))
	Yp = T(Y, I, translation)
	return errorfct(X, Yp), translation


def randrot(X, Y, errorfct):
	theta = random.uniform(0., 2.*np.pi)
	rotation = rotation_matrix_cp(random_unit_vector(), theta).T
	Z = np.matrix(np.zeros((4,1)))
	Yp = T(Y, rotation, Z)
	return errorfct(X, Yp), rotation

def ptSSE(pt, X):
	# the sum of the squares of the distances between each point in X and the nearest point in Y
	difference = pt - X
	xcol = np.ravel(difference[:, 0])
	ycol = np.ravel(difference[:, 1])
	zcol = np.ravel(difference[:, 2])
	cs_col = np.ravel(difference[:, 3])
	tol_matching = 2.5 # ppm
	cs_difference = sum(1-math.exp(-0.5*cs_col[i]**2/tol_matching**2) for i in xrange(len(cs_col)))
	sqr_difference = xcol**2. + ycol**2. + zcol**2. #+ 2.*cs_difference
	distance = np.min(sqr_difference)
	nearest_pt = np.argmin(sqr_difference)
	return distance

def NSSE(X, Y):
	# nearest sum squared error
	err = 0.
	for x in X:
		err+=ptSSE(x, Y)
	return err

def fit_clouds(X, Y, M, N, errorfct, threshold=1e-5):
	T0 = list()
	T1 = list()
	errors = list()
	errors.append(errorfct(X, Y))
	print errors[-1]
	
	Yp = Y.copy()
	for iter in range(M):
		err, translation = translate(X, Yp, errorfct)
		if err < threshold:
			break
		elif err<errors[-1]:
			errors.append(err)
			print errors[-1]
			T1.append(translation)
			I = np.matrix(np.identity(4))
			Yp = T(Yp, I, T1[-1])
		
		rot = [randrot(X, Yp, errorfct) for i in range(N)]
		rot.sort()
		err, rotation = rot[0]
		if err < threshold:
			break
		elif err < errors[-1]:
			errors.append(err)
			print errors[-1]
			T0.append(rotation)
			Z = np.matrix(np.zeros((4,1)))
			Yp = T(Yp, T0[-1], Z)
	
	return Yp, errors[-1]

def closest_point(X,Y, error):
	# print out the closest point in Y for each in X
	nb = 0
	nb_x = 0
	list_corres = []
	
	radius_lim = 1.5*(math.sqrt(error/float(len(X))))
	
	for x in X:
		min_ = float('inf')
		difference = x - Y
		xcol = np.ravel(difference[:, 0])
		ycol = np.ravel(difference[:, 1])
		zcol = np.ravel(difference[:, 2])
		cs_col = np.ravel(difference[:, 3])
		tol_matching = 2.5 # ppm
		cs_difference = sum(math.exp(-0.5*cs_col[i]**2/tol_matching**2) for i in xrange(len(cs_col)))
		sqr_difference = xcol**2. + ycol**2. + zcol**2. #+ 2.*cs_difference
		
		nb_corr_ = 0
		for nb_row in xrange(len(Y)):
			if sqr_difference[nb_row]<min_:
				min_ = sqr_difference[nb_row]
				nb_corr_ = nb_row				
			if math.sqrt(sum((x[i]-Y[nb_row][i])**2 for i in xrange(3)))<=radius_lim:
				nb+=1
				list_corres.append((nb_x, nb_row, math.sqrt(sum((x[i]-Y[nb_row][i])**2 for i in xrange(3)))))
		nb_x+=1
	return list_corres
