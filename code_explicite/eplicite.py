#! /usr/bin/env python
# -*- coding: utf-8 -*-
import numpy as np
import schemas as g
import matplotlib.pyplot as plt
from scipy.integrate import odeint
import sys
from mpl_toolkits.mplot3d import Axes3D
import math
reload (sys)
sys.setdefaultencoding('utf8')

a,b=0.,1.
T=0.5

N=30

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
    


    
xx,tt = np.meshgrid(t,x)
fig = plt.figure()
ax = fig.add_subplot(111,projection='3d')

z = g.u_ex(xx,tt)
diff = 0*u
for i in np.arange(0,P+1,1):
    diff[:,i] = np.abs(z[:,i]-u[:,i])

surf = ax.plot_surface(xx,tt,diff,linewidth=0)

# ax.set_xlabel('X', fontsize = 12)
# ax.set_ylabel('temps', fontsize = 12)
# ax.set_zlabel('u', fontsize = 12)
# ax.set_title('Solution exacte',fontsize = 16)

# surf = ax.plot_surface(xx,tt,z,linewidth=0)

# ax = fig.add_subplot(211,projection='3d')
# ax.set_xlabel('X', fontsize = 12)
# ax.set_ylabel('temps', fontsize = 12)
# ax.set_zlabel('u', fontsize = 12)
# ax.set_title('Solution approchée', fontsize = 16)
# surf = ax.plot_surface(xx,tt,u, linewidth=0)
plt.show()



        
    
