#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May  7 21:53:27 2020

@author: davidr
"""

#Plot the interior structure from Magrathea

#Required Packages
from matplotlib import pyplot as plt
import re


inf=open('../result/StructureEarth.txt','r')  #Takes output of input_mode=0
lines=inf.readlines()
#parameters from file radius, pressure, mass, density, temperature, phase
magrad=[]
magpress=[]
magmass=[]
magdens=[]
magtemp=[]
magphase=[]
for line in lines[1:-1]: #Skips first line and last empty line
    div=re.split(r'\s{2,}',line)
    magrad.append(float(div[1]))
    magpress.append(float(div[2]))
    magmass.append(float(div[3]))
    magdens.append(float(div[4]))
    magtemp.append(float(div[5]))
    magphase.append(div[6])
inf.close()

plt.figure(figsize=(6,10))

#Desity vs Radius Plot
plt.subplot(3,1,1)
plt.grid(color='lightgrey',linewidth=0.5,zorder=1)
plt.plot(magrad,magdens,label='R='+str(round(magrad[-1],4))) #Record radius in legend
plt.ylabel('Density', fontsize=14)
plt.legend(framealpha=1,fontsize=10)

#Pressure vs Radius
plt.subplot(3,1,2)
plt.grid(color='lightgrey',linewidth=0.5,zorder=1)
plt.plot(magrad,magpress)
plt.ylabel('Pressure', fontsize=14)

#Temperature vs Radius
plt.subplot(3,1,3)
plt.grid(color='lightgrey',linewidth=0.5,zorder=1)
plt.plot(magrad,magtemp)

plt.xlabel('Radius', fontsize=14)
plt.ylabel('Temperature', fontsize=14)

#Print Planet's radius
print(magrad[-1])

plt.savefig('interior.pdf',bbox_inches='tight') #Save figure as pdf
