#! /usr/bin/env python
# -*- coding: utf-8 -*-
import numpy as np
import fonction as g
import matplotlib.pyplot as plt
from scipy.integrate import odeint
import sys
from mpl_toolkits.mplot3d import Axes3D
import math
reload (sys)
sys.setdefaultencoding('utf8')
a=1
b=4
T=1
N=40
x = np.linspace(a,b,N+1)
dx=x[1]-x[0]
dt= 0.5*dx*dx
P=int(T/dt)
S=np.zeros((P+1,N+1))
t = np.linspace(0,T,P+1)
uo = g.init_U(a,N,dx)
S[:][0]=uo
Jn = g.initialise_Jn(dx,dt,N,a)
#print(Jn)
for n in np.arange(1,P+1,1):
    Fn = g.init_Fn(t[n],a,dx,N)
    A=np.zeros((1,N+1),dtype='d')
    A=g.prod_mat_creuse(Jn,uo,N)
    U1=A+Fn
    print("n=",n)
    print(U1)
    uo=U1
    S[:][n]=U1


SX,ST=np.meshgrid(x,t)
fig = plt.figure()
ax = fig.add_subplot(111,projection='3d')

surf = ax.plot_surface(SX,ST,S, linewidth=0)
plt.show()



        
    
