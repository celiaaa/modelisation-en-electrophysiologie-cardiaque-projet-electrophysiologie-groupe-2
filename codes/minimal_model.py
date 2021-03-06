# -*- coding: utf-8 -*-

import model as m

import scipy as sp
import numpy as np

import matplotlib.animation as animation
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.colors as colors
from mpl_toolkits.mplot3d import Axes3D
from math import ceil,floor
from scipy.sparse.linalg import splu
from scipy.sparse import csc_matrix , diags
from time import sleep
# Initialisation des paramètres
Omega = [0.,1.,0.,1.] # On modèlise sur le carré (0,1)x(0,1)
t = 100.
D = 1.e-3

N=20
n2=(N+1)*(1+N)
h = (Omega[1]-Omega[0])/float(N)
dt = h*h

P = int(t/dt)
T = np.linspace(0.,t,P+1)

X,Y = np.mgrid[Omega[0]:Omega[1]+h:h, Omega[2]:Omega[3]+h:h]

y0 = np.zeros((4,n2))
y1 = np.zeros((4,n2))


M,J = m.M_CN(h,dt,N,Omega)         # Matrice du schéma CN
B = splu(M)                     # Factorisation LU

# Conditions initiales
y0[0,:] = 0.  # u
y0[1,:] = 1.  # v
y0[2,:] = 1.  # w
y0[3,:] = 0.  # s 

# Stockage pour l'affiche
U = np.zeros((N,N))
S = np.zeros((3,P+1))
# S[0,0] = y0[1,40]
# S[1,0] = y0[2,40]
# S[2,0] = y0[3,40]

for n in np.arange(1,P+1,1):
    y1 = y0 + dt*m.G(y0,(n-1)*dt,X,Y)

    y1[0,:] = B.solve(J.dot(y1[0,:]))
    y0 = y1
    
    
    # S[0,n] = y1[1,40]
    # S[1,n] = y1[2,40]
    # S[2,n] = y1[3,40]
    
    # Affichage des résultats
    if n%10==0:
        U = np.reshape(y0[0,:],(N+1,N+1))
        U = m.U_mv(U)
        
        fig = plt.figure()
        # ax = fig.gca(projection='3d')
        # ax.set_xlabel('X', fontsize = 16)
        # ax.set_ylabel('Y', fontsize = 16)
        # ax.set_zlabel('Potentiel', fontsize = 16)
        
        # surf = ax.plot_surface(X,Y,U,cstride=1,linewidth=0,cmap='jet')
        ax = plt.subplot(111)
        ax, _ = mpl.colorbar.make_axes(plt.gca())
        norm=mpl.colors.Normalize(vmin=-85., vmax=65.)
        
        cax  = plt.matshow(U,vmin=-85.,vmax=60.,cmap='jet',norm=norm)
        
        
        # plt.contourf(X,Y,U,50,vmin=-85.,vmax=60.)
        
        plt.show()
        
        


# plt.plot(T,S[0,:])
# plt.show()

