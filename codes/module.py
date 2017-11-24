# -*- coding: utf-8 -*-


import numpy as np
import scipy.sparse as sp


def alpha(x):
    "Coefficient alpha(x) de l'équation"
    return 1.*(1.+0.*x) # pour répondre à la question en prenant alpha=0.01

def f(t,x):
    resultat=0# on commence avec 0
    return resultat

def u_0(x):
    resultat=np.cos(np.pi*x) # pour commencer
    return resultat

def u_ex(t,x):
    return np.cos(np.pi*x)*np.exp(-(np.pi)**2*t)



def init_Fn(t,a,pas_x,dimension):
    resultat=np.zeros((1,dimension+1),dtype='d')
    for i in np.arange(0,dimension+1,1):
        resultat[0][i]=f(t,a+i*pas_x)
    return resultat


# -------fonction pour le schéma d'euler implicite

def M_impl(dx,dt,n,a):
    """M_impl(dx,dt,n)
    
    Renvoie la matrice M = Id + (dt/dx^2)*A ou A = matrice de discrétisation de
    -d_xx(u), sur un intervalle [a,b] découpé en n morceaux de tailles
    dx. Il y a donc n+1 points. La matrice est creuse (tridiagonale) de
    taille n+1 x n+1. Avec les conditions limites d_x(u)=0
"""
    # La matrice est M = Id + dt*A
    return sp.identity(n+1,format="csc") + dt*A(dx,n,a)


def A(dx,n,a):
    """A(dx,n,a)
    
    Renvoie la matrice A de la discrétisation par DF de -d_xx(u), sur un
    intervalle [a,b] découpé en n morceaux de tailles dx. Il y a donc
    n+1 points. La matrice est creuse (tridiagonale) de taille n+1 x
    n+1. Avec les conditions limites d_x(u)=0
"""
    # La matrice est
    #  2*a_{1/2}   -2*a_{1/2}
    #  -a_{1/2}  a_{1/2}+a_{3/2} -a_{3/2}
    #              -a_{3/2}         .         .
    #                               .         .      .
    #                                       -2*a_{n-1/2}  2a_{n-1/2}
    x_ip12 = a + dx*(0.5+np.arange(n)) # x_{i+1/2} i = 0 à n-1
    a_ip12 = alpha(x_ip12)
    diag_0 = np.hstack([a_ip12,0]) + np.hstack([0,a_ip12])
    diag_0[0] = diag_0[0] + a_ip12[0]
    diag_0[n] = diag_0[n] + a_ip12[n-1]
    diag_m1 = -a_ip12
    diag_m1[n-1] = diag_m1[n-1] - a_ip12[n-1]
    diag_p1 = -a_ip12
    diag_p1[0] = diag_p1[0] - a_ip12[0]
    return 1./dx**2 * sp.diags([diag_0,diag_m1,diag_p1],[0,-1,+1],format="csc") # dia ou csc ou csr


# ------------ Foncitons pour la shcéma de Crank-Nicolson

def M_CN(dx,dt,n,a):
    """M_CN(dx,dt,n)
    
    Renvoie les matrices M = Id + (dt/2dx^2)*A et J = Id - (dt/2dx^2)*A ou A = matrice de discrétisation de
    -d_xx(u), sur un intervalle [a,b] découpé en n morceaux de tailles
    dx. Il y a donc n+1 points. La matrice est creuse (tridiagonale) de
    taille n+1 x n+1. Avec les conditions limites d_x(u)=0
"""
    # Les matrices sont M = Id + (dt/2)*A et J = Id - (dt/2)*A 
    return sp.identity(n+1,format="csc") + 0.5*dt*A(dx,n,a) , sp.identity(n+1,format="csc") - 0.5*dt*A(dx,n,a)




# ------------ Pour la résolution du système linéaire Mx^k = x^k-1

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


# def fact_LU(A):
#     print 'LU factorisation:'
#     plu = sp.lu_factor(A)
#     # if np.all(np.dot(P,np.dot(L,U))==A):
#     #     print 'Done correctly'
#     # else:
#     #     print 'Not correct'
    
#     return plu
