# -*- coding: utf-8 -*-

import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from mpl_toolkits.mplot3d import Axes3D
from math import ceil,floor
from scipy.sparse.linalg import splu

import module as mod

#Initialisation parametres
a,b = 0.,1.
T=0.5

N = input("Quelle valeur de N ?")
dx = (b-a)/N
dt=0.5*dx*dx
P=int(T/dt)

T = np.linspace(0,T,P+1)
X = np.linspace(a,b,N+1)
u = np.zeros((N+1,P+1))

u[:,0] = mod.u_0(X)
M = mod.M_impl(dx,dt,N,a)
B = splu(M)                     # Factorisation LU


for n in np.arange(1,P+1,1):
    u[:,n] = B.solve(u[:,n-1])




# Calcul de l'erreur

print "Erreur de la méthode d'Euler implicite : ",mod.erreur(u,T,X)

    # Affichage des résultats

# sx,st = sp.meshgrid(X,T)

# fig = plt.figure()
# ax = fig.gca(projection='3d')
# ax.set_xlabel('X', fontsize = 16)
# ax.set_ylabel('temps', fontsize = 16)
# ax.set_zlabel('temperature', fontsize = 16)

# norm = colors.Normalize(0,10)
# surf = ax.plot_surface(sx,st,np.transpose(u),cstride=1,linewidth=0,cmap='jet')
# cb = fig.colorbar(surf,ax = ax)
# plt.show()
