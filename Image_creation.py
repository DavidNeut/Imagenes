# -*- coding: utf-8 -*-
"""
Created on Mon Jul 21 13:04:41 2025

@author: Jimi_
"""

import numpy as np
from numpy import linalg as lg
import matplotlib.pyplot as plt
import skimage
from skimage import io
from numpy import random as rd

def random_image(Size):
    N0=Size[0]
    N1=Size[1]
    Blank=np.zeros([N0,N1,3])
    for i in range(N0):
        for j in range(N1):
            for k in range(3):
                Blank[i,j,k]=rd.uniform()
    return Blank


def random_Walk(N,K,steps):
    Blank=np.zeros((N,K,3))
    Visits=np.zeros((N,K))
    Blank[int(0.5*N),int(0.5*K)]=np.ones(3)
    P=np.array([int(0.5*N),int(0.5*K)])
    R=np.array([1,0])
    L=np.array([-1,0])
    U=np.array([0,-1])
    D=np.array([0,1])
    Dic={0:L,1:R,2:U,3:D}
    for s in range(steps):
        choose=[0,1,2,3]
        if P[0]==0:
            choose.remove(0)
        if P[0]==N-1:
            choose.remove(1)
        if P[1]==0:
            choose.remove(2)
        if P[1]==K-1:
            choose.remove(3)
        i=rd.randint(0,len(choose))
        P=P+Dic[choose[i]]
        #print(choose, i, P[0],P[1])
        
        
        #Blank[int(P[0]),int(P[1])]=Color[i]
        Visits[int(P[0]),int(P[1])]=Visits[int(P[0]),int(P[1])]+1
    return [Blank,Visits]

