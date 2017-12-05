# -*- coding: utf-8 -*-

import module as m

import scipy as sp
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from mpl_toolkits.mplot3d import Axes3D
from math import ceil,floor
from scipy.sparse.linalg import splu


# Initialisation des paramètres
a,b = 0. , 1. # On modèlise sur le segment (0,1)
t = 5.
D = 1.e-3

N=50
dx = (b-a)/float(N)
dt = dx*dx
P = int(t/dt)

T = np.linspace(0.,t,P+1)
X = np.linspace(a,b,N+1)

y0 = np.zeros((4,N+1))
y1 = np.zeros((4,N+1))

A = D*m.A(dx,N,a)

M,J = m.M_CN(dx,dt,N,a)         # Matrice du schéma Euler Implicite
B = splu(M)                     # Factorisation LU

# Conditions initiales
y0[0,:] = 0.  # u
y0[1,:] = 1.  # v
y0[2,:] = 1.  # w
y0[3,:] = 0.  # s 

# Stockage pour l'affiche
U = np.zeros((P+1,N+1))
U[0,:] = y0[0,:]
S = np.zeros((3,P+1))

for n in np.arange(1,P+1,1):
    y1 = y0 + dt*m.G(y0,(n-1)*dt)

    y1[0,:] = B.solve(J.dot(y1[0,:]))
    y0 = y1
    
    U[n,:] = y0[0,:]
    S[0,n] = y1[1,1]
    S[1,n] = y1[2,1]
    S[2,n] = y1[3,1]
    
    
    # Affichage des résultats

U = m.U_mv(U)
sx,st = sp.meshgrid(X,T)

fig = plt.figure()
ax = fig.gca(projection='3d')
ax.set_xlabel('espace', fontsize=12)
ax.set_ylabel('temps',fontsize=12)
ax.set_zlabel('potentiel',fontsize=14)

norm = colors.Normalize(np.max(U[n,:]),np.min(U[n,:]))
surf = ax.plot_surface(sx,st,U ,cstride=1,linewidth=0,cmap='jet')

cb = fig.colorbar(surf,ax=ax)

# plt.plot(T,S[2 ,:])
plt.show()
