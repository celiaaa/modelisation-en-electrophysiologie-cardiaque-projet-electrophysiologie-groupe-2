#! /usr/bin/env python
# -*- coding: utf-8 -*-
import numpy as np
import fonction as g
import matplotlib.pyplot as plt
from scipy.integrate import odeint
import sys  
reload (sys)
sys.setdefaultencoding('utf8')
a=1
b=4
T=1
N=5
P=10
x = np.linspace(a,b,N+1)
t = np.linspace(0,T,P+1)
print t[P]
dx = x[1]-x[0]
dt = t[1]
uo = g.init_U(a,N,dx)
Jn = g.initialise_Jn(dx,dt,N,a)
for n in np.arange(1,P+1,1):
    Fn = g.init_Fn(t[n],a,dx,N)
    A=np.zeros((1,N+1),dtype='d')
    A=g.prod_mat_creuse(Jn,uo,N)
    U1=A+Fn
    print(" pour n =",n," on a:")
    print(U1)
    uo=U1
    

        
    
