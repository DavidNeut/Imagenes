# -*- coding: utf-8 -*-
"""
Created on Mon Jul 21 12:43:30 2025

@author: DavidNeut
"""
import numpy as np
from numpy import linalg as lg
import matplotlib.pyplot as plt
import skimage
import numpy as np
from skimage import io


"""
Parte entera
"""
def proyeccion(Picture,u):
    """
    Parameters
    ----------
    Picture :Imagen a tratar
    u : coordenadas a interpolar

    Returns
    -------
    list
        Los valores de los coloress de la proyeccion (esquina)
    """
    x=int(u[0])
    y=int(u[1])
    W=Picture.shape
    if x>=0 and x<=int(W[0])-1:
        if y>=0 and y<=int(W[1])-1:
            return [True, Picture[x,y]]
    return [False,None]

#################################################

"""
Nodo más cercano
"""

def Closer(Picture,v):
    """
    Metodo de interpolacion: El nodo más cercano

    Parameters
    ----------
    Picture : La imagen
    v : Coordenada a encontrar el valor de la interpolción
    

    Returns
    -------
    list
        El valor correspondiente de su nodo más cercano
    """
    W=Picture.shape
    x=int(v[0])
    y=int(v[1])
    
    if x>=0.0 and x<=W[0]-1:
        if v[0]-x>x+1-v[0]:
            x=x+1
        if y>=0 and y<=W[1]-1:
            if v[1]-y>y+1-v[1]:
                y=y+1 
            return [True,Picture[x,y]]
    return [False, None]
            
###############################################
"""
Interpolacion Bilineal
"""
def Bilineal(Picture,u):    
    x=int(u[0])
    y=int(u[1])
    W=Picture.shape
    if x>=0 and x<=int(W[0])-2:
        if y>=0 and y<=int(W[1])-2:         
            
            v=(u[0]-x)*Picture[x,y]+(x+1-u[0])*Picture[x+1,y]
            w=(u[0]-x)*Picture[x,y+1]+(x+1-u[0])*Picture[x+1,y+1]
            return [True,(u[1]-y)*v+(y+1-u[1])*w]
    if x==int(W[0])-1:
        if y<=int(W[1])-2:
            return [True, (u[1]-y)*Picture[x,y]+(y+1-u[1])*Picture[x,y+1]]
        else:
            return [True, Picture[x,y]]
    if y==int(W[1])-1:
        if x<=int(W[0])-2:
            return [True, (u[0]-x)*Picture[x,y]+(x+1-u[0])*Picture[x+1,y]]
        else:
            return [True, Picture[x,y]]   
        
    return [False,None]               
            

###############################################

"""
Funciones para quue funcione la interpolación bicubica
Usando la derivada
"""
def submatrix(A,coord,lenght):
    W=A.shape
    if len(W)==2:
        Z=np.zeros(lenght)
    else:
        Z=np.zeros([lenght[0],lenght[1],W[2]])
    for i in range(lenght[0]):
        for j in range(lenght[1]):
            if i<W[0]-1-coord[0] and i+coord[0]>=0:
                if j<W[1]-1-coord[1] and j+coord[1]>=0:
                    Z[i,j]=A[coord[0]+i,coord[1]+j]
                else:
                    if j>=W[1]-1-coord[1]:
                        Z[i,j]=A[i+coord[0],W[1]-1]
                    else:
                        Z[i,j]=A[i+coord[0],0]
            else:
                if j<W[1]-1-coord[1] and j+coord[1]>=0:
                    Z[i,j]=A[W[0]-1,coord[1]+j]
                elif i>=W[0]-1-coord[0]:
                    if j>=W[1]-1-coord[1]:
                        Z[i,j]=A[W[0]-1,W[1]-1]
                    else:
                        Z[i,j]=A[W[0]-1,0]
                else:
                    if j>=W[1]-1-coord[1]:
                        Z[i,j]=A[0,W[1]-1]
                    else:
                        Z[i,j]=A[0,0]
                    
    return Z


def Fun(A,coord):
    return submatrix(A,coord,[2,2])
def derx(A,coord):
    D=submatrix(A,coord,[2,2])
    D_R=submatrix(A,[coord[0]+1,coord[1]],[2,2])
    D_L=submatrix(A,[coord[0]-1,coord[1]],[2,2])
    return 0.5*((D_L+D_R)-(D+D))
def dery(A,coord):
    D=submatrix(A,coord,[2,2])
    D_U=submatrix(A,[coord[0],coord[1]+1],[2,2])
    D_D=submatrix(A,[coord[0],coord[1]-1],[2,2])
    return 0.5*((D_U+D_D)-(D+D))

def derxy(A,coord):
    DX=derx(A,coord)
    DX_U=derx(A,[coord[0],coord[1]+1])
    DX_D=derx(A,[coord[0],coord[1]-1])
    DY=dery(A,coord)
    DY_R=dery(A,[coord[0]+1,coord[1]])
    DY_L=dery(A,[coord[0]-1,coord[1]])
    return 0.25*((DX_U+DX_D)-(DX+DX)+(DY_R+DY_L)-(DY+DY))

def Matrix_bicubic(A,coord):
    W=A.shape
    n=W[2]
    F=np.zeros([4,4,n])
    Fu=Fun(A,coord)
    Dx=derx(A, coord)
    Dy=dery(A, coord)
    Dxy=derxy(A, coord)
    for i in range(2):
        for j in range(2):
            F[i,j]=Fu[i,j]
            F[i+2,j]=Dx[i,j]
            F[i,j+2]=Dy[i,j]
            F[i+2,j+2]=Dxy[i,j]
    Q=[]
    for k in range(n):
        Q.append(np.dot(np.dot(Right_M,F[:,:,k]),Left_M))
        
    return Q
    
def polinomio_cubico(Q,u):
    n=len(Q)
    x=u[0]-int(u[0])
    y=u[1]-int(u[1])
    R=np.zeros(n)
    
    for k in range(n):
        for i in range(4):
            for j in range(4):
                R[k]=R[k]+pow(x,i)*pow(y,j)*Q[k][i,j]
        if R[k]<0.0:
            R[k]=0.0
        elif R[k]>1:
            R[k]=1
    
            
    return R
   
def bicubico(Picture,u):
    W=Picture.shape
    coord=[int(u[0]),int(u[1])]
    Q=Matrix_bicubic(Picture, coord)
    if int(u[0])<0 or int(u[1])<0:
        return [False, None]
    elif int(u[0])>W[0]-1 or int(u[1])>W[1]-1:
        return [False,None]
    return [True,polinomio_cubico(Q,u)]




###########################################################
"""
Proceso para interpolar las imagenes con cualquier metodo.
"""

def Corners(Picture,A):
    """
    Abarca los valores de las nuevas 4 esquinas
    de la imagen tranformada por la matriz A
    """
    W=Picture.shape
    S1=np.array((W[0],0))
    S2=np.array((0,W[1]))
    S3=np.array((W[0],W[1]))
    return [np.zeros(2),np.dot(A, S1),np.dot(A, S2),np.dot(A, S3)]
def New_size(Picture,A):
    """

    Parameters
    ----------
    Picture :Imagen a transformar
    A :Matriz de 2x2

    Returns
    -------
    list
        S es el nuevo tamaño de la imagen
        center es el nuevo origen de la imagen

    """
    C=Corners(Picture, A)
    X=[]
    Y=[]
    for i in range(4):
        X.append(C[i][0])
        Y.append(C[i][1])
    mx=min(X)
    my=min(Y)
    Mx=max(X)
    My=max(Y)
    S=[Mx-mx,My-my]
    center=np.array([-mx,-my])
    return [S,center]

def Inversa_f(v,A,center=np.zeros(2)):
    """
    Parameters
    ----------
    v : vector o coordenadas a transformar
    A : Matriz de 2x2.
    center : Es un arreglo de tamaño 2. Las coordenadas de la imagen del origen del plano cartesiano
        DESCRIPTION. The default is np.zeros(2).

    Returns
    -------
    TYPE
        La coordenada de la imagen de la transformación A.

    """
    w=v-center
    
    return np.dot(A,w)

def New_picture(Picture,A,method=proyeccion):
    """
    Parameters
    ----------
    Picture : Imagen a interpolar
    A : Matriz de 2x2, la transformación de la imagen
    method : El método de interpolación que se quiera usar.
        DESCRIPTION. The default is proyeccion.

    Returns
    -------
    Blank :La imagen ya aplicada la transformación lineal.

    """
    if abs(lg.det(A))<0.001:
        print("No se puede realizar")
        return None
    W=Picture.shape
    S=np.array([[W[0],0],[0,W[1]]])
    In=lg.inv(A)
    [S,Cen]=New_size(Picture, A)
    N0=int(S[0])
    N1=int(S[1])
    Blank=np.zeros((N0,N1,W[2]))
    for i in range(N0):
        for j in range(N1):
            w=np.array([i,j])
            v=Inversa_f(w,In,Cen)
            Proy=method(Picture, v)
            if Proy[0]:
                Blank[i,j]=Proy[1]
            
    return Blank
