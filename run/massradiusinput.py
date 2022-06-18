#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 17 18:14:12 2022

@author: davidr
"""


filename='water' #Give Descriptive file name
minmass=0.1  #Minimal mass
maxmass=4  #Maximum mass
mstep=0.1 #Number of data points

fcore=.325  #fraction mass for core, mantle, water. 1-fcore-fmantle-fwater=fatmosphere
fmantle=.675
fwater=0


f=open('input'+filename+'.txt','w')
f.write('Mass, fCore, fMantle, fWater\n')

mass=minmass
while mass<=maxmass:
    f.write(str(mass)+' '+ str(mass*fcore) +' '+ str(mass*fmantle) +' '+ str(mass*fwater) +'\n')
    mass=round(mass+mstep,12)
    
f.close()