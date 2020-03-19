#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 18 09:46:52 2020
This scipt is used to generate SPT inp file(s) for ABAQUS Simulations with 
specific parameters (thickness, materials parameters)
It reads an origin file ("SPT.inp") and modify the file accordingly

Multiple inp files with random repartition over a defined parameter space can be generated


@author: mowgli
"""

import numpy as np
import pandas as pd
from pyDOE import lhs

## functions to write inp files from input parameters 

def StrLength(y,length):
    if len(y)<length:
        y = ' '*(length-len(y))+y
    else:
        y = ' '+y[0:length-1]
    return y

def INP_NODES(NODES,t):
    NewNODES = ''
    while len(NODES)!=0:     #line per line: change of the y value of the node 
        y = float(NODES[22:35])*t/500
        y = StrLength(str(y),13)
        NewNODES = NewNODES + NODES[0:22]+y+'\n'
        NODES = NODES[36:]
    return NewNODES[:-1]

def INP_DIES(DIES,t):
    y = str(1.2+(float(DIES[22:35])-1.2)*t/500)
    y = StrLength(str(y),13)
    NewDIES = DIES[0:22] + y + DIES[35:185]
    y = StrLength(str(1.2+(float(DIES[185:198])-1.2)*t/500),13)
    NewDIES = NewDIES + y + DIES[198:220] + StrLength(str(t/1000),13) + DIES[233:248] + StrLength(str(1.2+t/1000),13)
    NewDIES = NewDIES + DIES[261:949] + StrLength(str(t/1000),13) + DIES[962:1108]+ StrLength(str(t/1000),13) 
    NewDIES = NewDIES + DIES[1121:1142] + StrLength(str(t/1000),13) + DIES[1155:]
    return NewDIES

def INP_Material(MAT,E,nu,K,E0,n):
    NewMAT = MAT[0:50] + str(E) + MAT[57:59] + str(nu) + MAT[63:73]
    x=np.power(np.linspace(0,np.sqrt(2),100),2)
    y=K*np.power(E0+x,n)
    for i in np.array([y,x]).T:
        NewMAT = NewMAT + StrLength(str(i[0]),8) + ',' + StrLength(str(i[1]),12) + '\n'
    return NewMAT[:-1] + MAT[80:]

def INP_Generation(INP,E,nu,K,E0,n,t,i,Mat): # i: simulation number (Job Name)
    NewINP = INP[0:22] + Mat + '_SPT_{:04d}'.format(i) + INP[25:47]   # initial section of the inp file 
    NewINP = NewINP + '\n** E = '+str(E)+'\n** nu = '+str(nu)+'\n** K = '+str(K)
    NewINP = NewINP + '\n** e0 = '+str(E0)+'\n** n = '+str(n)+'\n** t = '+str(t) + INP[47:245]
    NewINP = NewINP + INP_NODES(INP[245:28360],t)   # nodes location as a function of t
    NewINP = NewINP + INP[28360:47929]
    NewINP = NewINP + INP_DIES(INP[47929:49190],t)
    NewINP = NewINP + INP_Material(INP[49190:49271],E,nu,K,E0,n)
    NewINP = NewINP + INP[49271:]
    return NewINP

def WriteINP(NewInp,fold,Mat,i):
    #NewInp = NewInp.encode('ascii')
    with open(fold+'/'+Mat+'_SPT_{:04d}.inp'.format(i),'w') as file:
        file.write(NewInp)
    return "File writen"


# Principal code
    
# Open the reference inp file 
with open('SPT.inp', 'r') as file:
    INP = file.read()

# Writing directory 
fold = 'Generated_inp'

# Material Name
Mat = input("Material Name : ")

# Ask for multiple or simple inp files
Choice = input("Multiple files (m) or Single file (s)? - (m/s) : ")

if Choice=='s':
    E = float(input("Modulus (GPa) : "))
    nu = float(input("Poisson ratio : "))
    K = float(input("K : "))
    E0 = float(input("E0 : "))
    n = float(input("n : "))
    t = float(input("Thickness (um) : "))
    i = int(input("Index : "))

    NewInp = INP_Generation(INP,E*1000,nu,K,E0,n,t,i,Mat)
    print(WriteINP(NewInp,fold,Mat,i))
    
elif Choice=='m':
    # Number of files
    N = int(input("Number of inp files to generate : "))
    Erange = np.array(input('Range of Modulus (GPa) : ').split()).astype(float)
    vrange = np.array(input('Range of Poisson ratio : ').split()).astype(float)
    Krange = np.array(input('Range of K : ').split()).astype(float)
    E0range = np.array(input('Range of E0 : ').split()).astype(float)
    nrange = np.array(input('Range of n : ').split()).astype(float)
    trange = np.array(input('Range of Thickness (um) : ').split()).astype(float)
    
    # Latin Hypercube sampling
    x = lhs(6,N)
    # LHC sample at corresponding values (from range)
    X=np.array(x).T
    X[0] = X[0] * (Erange[1]-Erange[0]) + Erange[0]
    X[1] = X[1] * (vrange[1]-vrange[0]) + vrange[0]
    X[2] = X[2] * (Krange[1]-Krange[0]) + Krange[0]
    X[3] = X[3] * (E0range[1]-E0range[0]) + E0range[0]
    X[4] = X[4] * (nrange[1]-nrange[0]) + nrange[0]
    X[5] = (X[5] * (trange[1]-trange[0]) + trange[0]).astype(int)
    X = pd.DataFrame(X.T,columns=['E','nu','K','E0','n','t'])
    
    # Save spt_design
    filename = fold+'/'+Mat+'_SPT_design.csv'
    with open(filename,mode='w') as file:
        X.to_csv(file,index=False)
        
    for i in np.arange(N)+1:
        [E, nu, K, E0, n ,t] = (X.iloc[i-1]).values
        NewInp = INP_Generation(INP,E*1000,nu,K,E0,n,int(t),i,Mat)
        WriteINP(NewInp,fold,Mat,i)
        
    print("All inp files written")

else:
    print("Wrong input : [s/m] - Respect letter case")