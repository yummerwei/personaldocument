import numpy as np
a = np.loadtxt('Chessian.dat')
M=np.mat(a)
c=np.linalg.eig(M)
d=c[1]
e=np.dot(np.dot(np.linalg.inv(d),a),d)
f=np.diagonal(e)
np.savetxt('diaghessian.txt',f,newline='\n')
