"""
	Eigen decomposition of 3x3 symmetric matrices.
	
	Author: Alexis Mignon (c)
	e-mail: alexis.mignon@gmail.com

"""
import numpy as np
cimport numpy as np

cdef extern :
	int eigen_sym_3D(double* A, double* w, double* U)

def eigh(np.ndarray[np.float_t, ndim = 2] A):
	"""
		Eigen decomposition of 3x3 symmetric matrices.
		
		arguments:
		- A: a 3x3 array of np.float
		
		returns: 
		- w: array of eigen-values sorted in decreasing order
		- u: array of eigen-vectors where u[:,i] is the eigen-vector
		     corresponding to w[i]
	"""
	cdef :
	
		np.ndarray[np.float_t, ndim = 2] A_ = np.zeros((3,3),'float')
		np.ndarray[np.float_t, ndim = 1] w = np.zeros(3,'float')
		np.ndarray[np.float_t, ndim = 2] U = np.zeros((3,3),'float')
		double *dA
		double *dw = <double*>w.data
		double *dU = <double*>U.data
		int res
		
	if (<object>A).flags.c_contiguous and (<object>A).dtype == "float64" :
		A_ = A
	else :
		A_ = A.copy()
		
	dA = <double*>A_.data
	
	res = eigen_sym_3D(dA, dw, dU)
	if res != 0 : raise ValueError("An error occured, probably due to some bad properties of the matrix A")
	return w,U
