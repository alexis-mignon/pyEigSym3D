import numpy as np
from scipy.linalg import eigh
from eigsym3d import eigh as eigh3D
import time

def test_speed():
    ntest = 1000
    X = np.random.randn(5,3)
    C = np.dot(X.T,X)
    start = time.time()
    for i in range(ntest):
        w,u = eigh(C)
    time0 = (time.time() - start)/ntest
    print "Elapsed time for scipy.linalg.eigh: %fs"%(time0)
    
    start = time.time()
    for i in range(ntest):
        w,u = eigh3D(C)
    time1 = (time.time() - start)/ntest
    print "Elapsed time for eigsym3d.eigh: %fs"%(time1)
    print "eigsym3d was", time0/time1, "faster"
    
if __name__ == '__main__' :
    test_speed()
