# -*- coding: utf-8 -*-

import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from mpl_toolkits.mplot3d import Axes3D
from math import ceil,floor
import module as mod

#Initialisation parametres
CFL = 0.5
t_fin = 1. 
nx = 30
X = np.linspace(0,1,nx+1)
dx = X[1]-X[0]

dt = CFL*dx*mod.alpha(0)
nt = int(floor(t_fin/dt))
T =  np.linspace(0,t_fin,nt)

u = np.zeros((np.size(X),np.size(T)))

u[:,0] = mod.u0(X) #On initialise u pour t=0

A = mod.init_A(nx,dx) #initialisation de A

M = np.eye(nx+1) + (dt/dx/dx)*A #Matrice M à inverser pour résoudre

# On calcule la première itération avec euler explicite
a,b,c = mod.init_abc(M,nx) #récupère les diago de M dans des vecteurs
d = u[:,0]
u[:,1] = np.transpose(mod.TDMAsolver(a,b,c,d))

for k in range(1,nt-1):
    
    u[:,k+1] = u[:,k-1] - 2*(dt/dx/dx)*A.dot(u[:,k]) 


sx,st = sp.meshgrid(X,T)

fig = plt.figure()
ax = fig.gca(projection='3d')
ax.set_xlabel('X', fontsize = 16)
ax.set_ylabel('temps', fontsize = 16)
ax.set_zlabel('temperature', fontsize = 16)

norm = colors.Normalize(0,10)
surf = ax.plot_surface(sx,st,np.transpose(u),cstride=1,linewidth=0,cmap='jet')
cb = fig.colorbar(surf,ax = ax)
plt.show()
