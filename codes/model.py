# -*- coding: utf-8 -*-


import numpy as np
import scipy.sparse as sp
import poisson as p

def alpha(x):
    "Coefficient alpha(x) de l'équation"
    return 1.*(1.+0.*x) # pour répondre à la question en prenant alpha=0.01

def A(dx,n,a):
    """A(dx,n,a)
    
    Renvoie la matrice A de la discrétisation par DF de -d_xx(u), sur un
    intervalle [a,b] découpé en n morceaux de tailles dx. Il y a donc
    n+1 points. La matrice est creuse (tridiagonale) de taille n+1 x
    n+1. Avec les conditions limites d_x(u)=0
"""
    # La matrice est
    #  2*a_{1/2}   -2*a_{1/2}
    #  -a_{1/2}  a_{1/2}+a_{3/2} -a_{3/2}
    #              -a_{3/2}         .         .
    #                               .         .      .
    #                                       -2*a_{n-1/2}  2a_{n-1/2}
    x_ip12 = a + dx*(0.5+np.arange(n)) # x_{i+1/2} i = 0 à n-1
    a_ip12 = alpha(x_ip12)
    diag_0 = np.hstack([a_ip12,0]) + np.hstack([0,a_ip12])
    diag_0[0] = diag_0[0] + a_ip12[0]
    diag_0[n] = diag_0[n] + a_ip12[n-1]
    diag_m1 = -a_ip12
    diag_m1[n-1] = diag_m1[n-1] - a_ip12[n-1]
    diag_p1 = -a_ip12
    diag_p1[0] = diag_p1[0] - a_ip12[0]
    return 1./dx**2 * sp.diags([diag_0,diag_m1,diag_p1],[0,-1,+1],format="csc") # dia ou csc ou csr

def M_CN(dx,dt,n,Omega):
    """M_CN(dx,dt,n)
    
    Renvoie les matrices M = Id + (dt/2dx^2)*A et J = Id - (dt/2dx^2)*A ou A = matrice de discrétisation de
    -d_xx(u), sur un intervalle [a,b] découpé en n morceaux de tailles
    dx. Il y a donc n+1 points. La matrice est creuse (tridiagonale) de
    taille n+1 x n+1. Avec les conditions limites d_x(u)=0
"""
    # Les matrices sont M = Id + (dt/2)*A et J = Id - (dt/2)*A 
    return sp.identity((n+1)*(n+1),format="csc") + 0.5*dt*p.matrix_neumann2D(Omega,n,n) , sp.identity((n+1)*(n+1),format="csc") - 0.5*dt*p.matrix_neumann2D(Omega,n,n)

def heaviside(x):
    return 1*(x>0)
def vinf(x):
    return 1*(x-thetamv<0)
def taumv(u):
    return (1-heaviside(u-thetamv))*taumv1+heaviside(u-thetamv)*taumv2
def taumw(u):
    return taumw1+(taumw2-taumw1)*(1+np.tanh(kmw*(u-umw)))/2
def tauso(u):
    return tauso1+(tauso2-tauso1)*(1+np.tanh(kso*(u-uso)))/2
def taus(u):
    return (1-heaviside(u-thetaw))*taus1+heaviside(u-thetaw)*taus2
def tauo(u):
    return (1-heaviside(u-thetao))*tauo1+heaviside(u-thetao)*tauo2
def winf(u):
    return (1-heaviside(u-thetao))*(1-u/tauwinf)+heaviside(u-thetao)*wainf
def Jfi(u,v):
    return -v*heaviside(u-thetav)*(u-thetav)*(uu-u)/taufi
def Jso(u):
    return (u-uo)*(1-heaviside(u-thetaw))/tauo(u) + heaviside(u-thetaw)/tauso(u)
def Jsi(u,w,s):
    return -heaviside(u-thetaw)*w*s/tausi
def I_app(t):
    return 1.*(4.5<t<5.7)#*(x<=0.15)*(x>=0.)*(y<=0.15)*(y>=0.)

def I(u,v,w,s,t,x,y):
    return -(Jfi(u,v)+Jso(u)+Jsi(u,w,s)) + I_app(t)
def g_v(u,v,w,s):
    return (1-heaviside(u-thetav))*(vinf(u)-v)/taumv(u) - heaviside(u-thetav)*v/taupv
def g_w(u,v,w,s):
    return (1-heaviside(u-thetaw))*(winf(u)-w)/taumw(u) - heaviside(u-thetaw)*w/taupw
def g_s(u,v,w,s):
    return ((1+ np.tanh(ks*(u-us)))/2.-s)/taus(u)

def G(y,t,x,z):
    u = y[0,:]
    v = y[1,:]
    w = y[2,:]
    s = y[3,:]
    return np.stack((I(u,v,w,s,t,x,z), g_v(u,v,w,s), g_w(u,v,w,s), g_s(u,v,w,s)))

uo=0
uu=1.55
thetav=0.3
thetaw=0.13
thetamv=0.006
thetao=0.006
taumv1=60
taumv2=1150
taupv=1.4506
taumw1=60
taumw2=15
kmw=65
umw=0.03
taupw=200
taufi=0.11
tauo1=400
tauo2=6
tauso1=30.0181
tauso2=0.9957
kso=2.0458
uso=0.65
taus1=2.7342
taus2=16
ks=2.0994
us=0.9087
tausi=1.8875
tauwinf=0.07
wainf=0.94

def U_mv(u):
    return 85.7*u - 84.
