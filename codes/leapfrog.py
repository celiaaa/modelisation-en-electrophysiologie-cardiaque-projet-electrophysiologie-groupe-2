#! /usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from mpl_toolkits.mplot3d import Axes3D
from math import ceil

c = 1.
xmin = 0.
xmax = 1.
N = 100
tmax = 10.
time_to_save = 20

dx = (xmax-xmin)/N
dt = dx*dx/2./c
a = 2.*dt/dx/dx
nb_iter = int(ceil(tmax/dt))

print 'dx / dt', dx , '  /  ',dt
print 'T = ',tmax, ' | ' , nb_iter, 'itérations'

X = np.linspace(xmin,xmax,N+1)
Y = np.linspace(0.,tmax,nb_iter+1)
SX, ST = np.meshgrid(X, Y)

sol = np.zeros((N+1,nb_iter))

sol[N/2:0] = 100.
print sol.shape
sol[0,1] = sol[0,0] + a*(sol[1,0] - sol[0,0])
sol[N,1] = sol[N,0] + a*(sol[N,0] - sol[N-1,0])
for i in range(1,N-1):
    sol[i,1] = sol[i,0] + (dt/dx/dx)*(sol[i+1,0] - 2*sol[i,0] + sol[i-1][0])



for k in range(1,nb_iter-1):

    sol[0,k+1] = sol[0,k-1] + 2*a*(sol[1,k]-sol[0,k])
    sol[N,k+1] = sol[N,k-1] + 2*a*(sol[N-1,k] - sol[N,k])
    sol[1:N-1,k+1] = sol[1:N-1,k-1]*(1-a)/(1+a) + (a/(1+a))*(sol[0:N-2,k]+sol[2:N,k])


fig = plt.figure()
ax = fig.gca(projection='3d')
surf = ax.plot_surface(SX,ST,np.transpose(sol))


plt.show()
