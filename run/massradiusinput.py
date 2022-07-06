#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 17 18:14:12 2022

@author: davidr
"""


filename='core' #Give Descriptive file name
minmass=0.1  #Minimal mass
maxmass=4  #Maximum mass
mstep=0.1 #Step in mass

fcore=1.0  #fraction mass for core, mantle, water. 1-fcore-fmantle-fwater=fatmosphere
fmantle=0.0
fwater=0.0


f=open('input'+filename+'.txt','w')
f.write('Mass, fCore, fMantle, fWater\n')

mass=minmass
while mass<=maxmass:
    f.write(str(mass)+' '+ str(fcore) +' '+ str(fmantle) +' '+ str(fwater) +'\n')
    mass=round(mass+mstep,12)
    
f.close()