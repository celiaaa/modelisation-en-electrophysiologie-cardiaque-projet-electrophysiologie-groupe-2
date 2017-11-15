# -*- coding: utf-8 -*-


import numpy as np


def u0(x): #u(0,x)=u0(x) on ]0,1[

    r = np.cos(np.pi*x)
    # r = 0*x
    # for i in range(0,len(x)):
    #     if 0.4<x[i]<0.6:
    #         r[i] = 10.

    return r

def alpha(x):

    return 0.01


## Tri Diagonal Matrix Algorithm(a.k.a Thomas algorithm) solver
def TDMAsolver(a, b, c, d):
    '''
    TDMA solver, a b c d can be NumPy array type or Python list type.
    refer to http://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm
    and to http://www.cfd-online.com/Wiki/Tridiagonal_matrix_algorithm_-_TDMA_(Thomas_algorithm)
    '''
    nf = len(d) # number of equations
    ac, bc, cc, dc = map(np.array, (a, b, c, d)) # copy arrays
    for it in xrange(1, nf):
        mc = ac[it-1]/bc[it-1]
        bc[it] = bc[it] - mc*cc[it-1] 
        dc[it] = dc[it] - mc*dc[it-1]
        	    
    xc = bc
    xc[-1] = dc[-1]/bc[-1]

    for il in xrange(nf-2, -1, -1):
        xc[il] = (dc[il]-cc[il]*xc[il+1])/bc[il]

    return xc


def init_A(n,dx): #return an array size (n+1)x(n+1) (0:n,0:n)

    A = np.zeros((n+1,n+1))
    for i in range(1,n):
        A[i,i] = alpha((i+0.5)*dx)+alpha((i-0.5)*dx)
        A[i,i+1] = -1*alpha((i+0.5)*dx)
        A[i,i-1] = -1*alpha((i-0.5)*dx)
    A[0,0] = 2*alpha(0.5*dx)
    A[n,n] = 2*alpha((n-0.5)*dx)
    A[0,1] = -2*alpha(0.5*dx)
    A[n,n-1] = -2*alpha((n-0.5)*dx)
    
    return A


def init_abc(A,n): #return 3 vectors of size n+1 wich contains the diagonals of A. a and c (lower and upper diagonals) are completed with 0.

    a = np.zeros((n+1,1))
    b = np.zeros((n+1,1))
    c = np.zeros((n+1,1))
    for i in range(1,n):
        a[i] = A[i,i-1]
        b[i] = A[i,i]
        c[i] = A[i,i+1]
    b[0] = A[0,0]
    b[n] = A[n,n]

    return a,b,c

