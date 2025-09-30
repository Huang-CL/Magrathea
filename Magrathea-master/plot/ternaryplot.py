#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 22 16:21:27 2022

@author: davidr
"""

#Plot a Ternary-Diagram of planet radii for given planet mass

#Required Packages
import ternary  #Install from https://github.com/marcharper/python-ternary
from matplotlib import pyplot as plt
import re
#from astropy import constants as const  #Use if want Earth-radii instead of km
import cmasher as cmr


plt.rcParams.update({'font.size': 12})

#5151 planets of one-earth mass with integer differences in mass percentage in each layer
inf=open('../result/OneEarthTern.txt','r')  #Output from input_mode=6

lines=inf.readlines()
coreid=[]
mantleid=[]
iceid=[]
radius=[]

mass=1.0

for line in lines[1:]:
    ##Get percent core, mantle, and ice
    div=re.split('\s{2,}',line)
    coreid.append(int(round(float(div[0])/mass*100,2)))
    mantleid.append(int(round(float(div[1])/mass*100,2)))
    iceid.append(int(round(float(div[2])/mass*100,2)))
    prevrad=float(div[7])
    radius.append(float(div[7]))  #*const.R_earth.si.value/1000 to change from Earth-Radii to km

#python ternary uses a dictionary ((bottom,right,left):float)
index=list(zip(coreid,iceid,mantleid))  #Rotation of Zeng & Seager 08
data=dict((index[i],radius[i]) for i in range(len(index))) 
#core to mantle ratio increases along the bottom axes
#the water percentage increases perpendicular to the bottomaxis

#Set-up Figure
scale = 100  #Must be the length of side and number of possible indexes must be s*(s+1)/2
figure,tax=ternary.figure(scale=scale)
figure.set_size_inches(6, 5.5)

# Draw Gridlines
tax.boundary(linewidth=1)
tax.gridlines(color="black", multiple=5.0)

#HeatMap
mycmap = plt.cm.get_cmap('cmr.rainforest')
cb_kwargs={"shrink" : 0.8,"orientation":"horizontal","pad" : 0.02}
tax.heatmap(data,cmap=mycmap,cb_kwargs=cb_kwargs,style='triangular') 
plt.text(45,-30, r'R$_\oplus$',fontsize=16)

# Set Axis labels and Title
axes_colors = {'b': 'k', 'l': 'k', 'r': 'k'}
plt.text(-22.5,-2,'Mantle',fontsize=16)
plt.text(103,-7,'Core',fontsize=16)
plt.text(50,92,'Water',fontsize=16)

# Set ticks
tax.ticks(axis='lbr', linewidth=1, multiple=10, offset=0.025,axes_colors=axes_colors)
tax.boundary(linewidth=1.0, axes_colors=axes_colors)

# Remove default Matplotlib Axes
tax.clear_matplotlib_ticks()
tax.get_axes().axis('off')

#Add descriptive text
plt.text(1,80, str(mass)+r' M$_\oplus$',fontsize=16)

tax.get_axes().set_aspect(1)
tax._redraw_labels()
#tax.show()
tax.savefig('oneearthternary.pdf')
tax.close()
inf.close()

