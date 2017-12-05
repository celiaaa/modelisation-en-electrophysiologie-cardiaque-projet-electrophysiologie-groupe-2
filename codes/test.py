# -*- coding: utf-8 -*-

import model as m

import scipy as sp
import numpy as np
import matplotlib.pyplot as plt

t=500.
P = 13680

T = np.linspace(0.,t,P+1)
dt = T[1]-T[0]
print "dt : ", dt
y = np.zeros((4,P+1))

y[0,0] = 0.
y[1,0] = 1.
y[2,0] = 1.
y[3,0] = 0.

for n in np.arange(0,P,1):
    y[:,n+1] = y[:,n] + dt*m.G(y[:,n],n*dt,0.)

y[0,:] = m.U_mv(y[0,:])

plt.plot(T,y[0,:])
# plt.plot([0.90,0.90],[-87,80.],c='r',linewidth=0.5,label='Debut de la stimulation')
# plt.legend(loc='upper left')
plt.show()
