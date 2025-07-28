# -*- coding: utf-8 -*-
"""
Created on Mon Jul 21 13:33:50 2025

@author: Jimi_
"""

import numpy as np
from numpy import linalg as lg
from Interpolacion import proyeccion,Closer
def Zoom_matrix(lam):
    return np.diag([lam,lam])



def Zoom(Picture,lam,method=proyeccion):
    W=Picture.shape
    Blank=np.zeros([int(lam*W[0]),int(lam*W[1]),W[2]])
    for i in range(W[0]):
        for j in range(W[1]):
            v=(1.0/lam)*np.array([i,j])
            Blank[i,j]=method(Picture,v)
    return Blank

def Zoom_Area(Picture,Area,lam,method=Closer):
    S0=Area[0]
    S1=Area[1]
    ZoomArea=Picture[S0[0]:S0[1],S1[0]:S1[1]]
    return Zoom(ZoomArea,lam,method)
    
def Transpose(Picture):
    W=Picture.shape
    [N,K,c]=[W[0],W[1],W[2]]
    Blank=np.zeros((K,N,c))
    for i in range(N):
        for j in range(K):
            Blank[j,i]=Picture[i,j]
    return Blank
def Mirror_Transpose(Picture):
    W=Picture.shape
    [N,K,c]=[W[0],W[1],W[2]]
    Blank=np.zeros((K,N,c))
    for i in range(N):
        for j in range(K):
            Blank[j,i]=Picture[N-1-i,K-1-j]
    return Blank
def Mirror_H(Picture):
    W=Picture.shape
    [N,K,c]=[W[0],W[1],W[2]]
    Blank=np.zeros((N,K,c))
    for i in range(N):
        for j in range(K):
            Blank[i,j]=Picture[i,K-1-j]
    return Blank
def Mirror_V(Picture):
    W=Picture.shape
    [N,K,c]=[W[0],W[1],W[2]]
    Blank=np.zeros((N,K,c))
    for i in range(N):
        for j in range(K):
            Blank[i,j]=Picture[N-1-i,j]
    return Blank