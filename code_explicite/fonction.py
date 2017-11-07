#! /usr/bin/env python
# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
import sys  
reload (sys)
sys.setdefaultencoding('utf8')

def alpha(x):
    resultat=1  # pour répondre à la question en prenant alpha=1
    return resultat

def initialise_Jn(pas_x,pas_t,dimension,a): # a est la borne inférieur de l'intervalle
    resultat=np.zeros((dimension+1,dimension+1),dtype='d')
    for i in np.arange(1,dimension,1):
        resultat[i][i]=1-((pas_t)/pas_x*pas_x)*(alpha(a+(i-0.5)*pas_x)+ alpha(a+(i+0.5)*pas_x))
        resultat[i][i-1]=(pas_t)/((pas_x)*(pas_x)) * alpha(a+(i-0.5)*pas_x)
        resultat[i][i+1]= (pas_t)/((pas_x)*(pas_x)) * alpha(a+(i+0.5)*pas_x)
    resultat[0][0]=1-2*((pas_t)/(pas_x*pas_x))*alpha(a+pas_x/2)
    resultat[0][1]=2*((pas_t)/(pas_x*pas_x))*alpha(a+pas_x/2)
    resultat[dimension][dimension-1]=2*(pas_t/pas_x*pas_x)*alpha(a+(dimension-0.5)*pas_x)
    resultat[dimension][dimension]=1-2*(pas_t/pas_x*pas_x)*alpha(a+(dimension-0.5)*pas_x)
    return resultat

def f(t,x):
    resultat=0 # on commence avec 0
    return resultat

def init_Fn(t,a,pas_x,dimension):
    resultat=np.zeros((1,dimension+1),dtype='d')
    for i in np.arange(0,dimension+1,1):
        resultat[0][i]=f(t,a+i*pas_x)
    return resultat


def init_uo(x):
    resultat=1 # pour commencer
    return resultat

def init_U(a,dimension,pas_x):
    resultat=np.zeros((1,dimension+1),dtype='d')
    for i in np.arange(0,dimension+1,1):
        resultat[0][i]=init_uo(a+i*pas_x)
    return resultat

def prod_mat_creuse(A,x,dimension):
    resultat=np.zeros((1,dimension+1),dtype='d')
    n=dimension
    resultat[0][0]=A[0][0]*x[0][0]+A[0][1]*x[0][1]
    resultat[0][n]=A[n-1][n-1]*x[0][n-1]+A[n-1][n-1]*x[0][n-1]
    for i in np.arange(1,n,1):
        resultat[0][i]=A[i][i-1]*x[0][i-1]+A[i][i]*x[0][i]+A[i][i+1]*x[0][i+1]
    return resultat   
