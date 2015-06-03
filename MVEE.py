# Python implementation of Nima Moshtagh's 
# Minimum Volume enclosing ellipsoids MATLAB code
# 
# By Haruki Hirasawa
# April 2015

import numpy as np
import numpy.linalg as linalg

# Returns A matrix for an ellipse and the center c
# such that (x-c)^T * A * (x-c) = 1
# defines the ellipse.
# P = d by N numpy array of N points in R^d
def MVEE(P,tol = 0.01):
	d, N = P.shape
	
	Q = np.zeros((d+1,N))
	Q[0:d,:] = P
	Q[d,:] = np.ones((1,N))

	count = 1
	err = 1
	u = (1.0/N)*np.ones(N)
	u = u.transpose()

	while err > tol:
		X = np.dot(np.dot(Q,np.diag(u)),np.transpose(Q))
		try:
			invX = linalg.inv(X)
			M = np.diag(np.dot(np.dot(np.transpose(Q),invX),Q))
			maxM = np.amax(M)
			j = np.argmax(M)
			step_size = (maxM - d - 1)/((d+1)*(maxM - 1))
			new_u = (1-step_size)*u
			new_u[j] = new_u[j] +step_size
			count += 1
			err = linalg.norm(new_u-u)
			u = new_u
		except: 
			print "Error"
			return None
	MatU = np.diag(u)

	u.transpose()
	c = np.dot(P,u)
	try: 
		A = (1./d)*linalg.inv(np.dot(np.dot(P,MatU),np.transpose(P))-np.multiply.outer(c,c))
	except:
		print "Error"
		return None
	return [A,c]
