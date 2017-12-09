#! /usr/bin/env python
# -*- coding: utf-8 -*-
import numpy as np
import schemas as g
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from scipy.integrate import odeint
import sys
from mpl_toolkits.mplot3d import Axes3D
import math
reload (sys)
sys.setdefaultencoding('utf8')

a,b=0.,1.
T=0.3

N = input("Quelle valeur pour N : ")

dx=(b-a)/N

dt= 0.5*dx*dx
P=int(T/dt)

t = np.linspace(0,T,P+1)
x = np.linspace(a,b,N+1)
u = np.zeros((N+1,P+1))

u[:,0] = g.u_0(x)

Jn = g.M_expl(dx,dt,N,a)
#print(Jn)
for n in np.arange(1,P+1,1):
    u[:,n] = Jn*u[:,n-1] + dt*g.f(t[n-1],x)
    


    
xx,tt = np.meshgrid(x,t)

z = g.u_ex(xx,tt)
diff = 0*u
# Affichage des résulats

fig = plt.figure()
ax = fig.gca(projection='3d')


ax.set_xlabel('X', fontsize = 12)
ax.set_ylabel('temps', fontsize = 12)
ax.set_zlabel('u', fontsize = 12)

# --------Affiche la solution exacte
#ax.set_title('Solution exacte',fontsize = 16)
#surf = ax.plot_surface(xx,tt,z,linewidth=0, cmap='jet')

# --------Affiche la solution approchée
ax.set_title('Solution approchée', fontsize = 16)
surf = ax.plot_surface(xx,tt,np.transpose(u),linewidth=0, cmap='jet')


norm = colors.Normalize(-1,1)
cb = fig.colorbar(surf,ax=ax)

plt.show()



        
    
