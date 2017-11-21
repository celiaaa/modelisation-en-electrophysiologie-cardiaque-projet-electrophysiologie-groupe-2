#! /usr/bin/env python
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

def M_expl(dx,dt,n,a):
    """M_expl(dx,dt,n)
    
    Renvoie la matrice M = Id - dt*A ou A = matrice de discrétisation de
    -d_xx(u), sur un intervalle [a,b] découpé en n morceaux de tailles
    dx. Il y a donc n+1 points. La matrice est creuse (tridiagonale) de
    taille n+1 x n+1. Avec les conditions limites d_x(u)=0
"""
    # La matrice est M = Id - dt*A
    return sp.identity(n+1) - dt*A(dx,n,a)
    
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
    return 1./dx**2 * sp.diags([diag_0,diag_m1,diag_p1],[0,-1,+1],format="dia") # dia ou csc ou csr
    
