'''
A temporary python plot script that draw the Magrathea interior structure profile saved in "../result/Structure.txt".
If multiple results exist in the file, this script shows the last profile.
'''

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from io import StringIO
from matplotlib import rc
import math

rc('text',usetex=True)          
mpl.rcParams['font.size']=16.5
mpl.rcParams["savefig.directory"] = ""

# convert to latex format scientific notation
# def latex_float(f):
#     float_str = "{0:.2g}".format(f)
#     if "e" in float_str:
#         base, exponent = float_str.split("e")
#         return r"${0} \times 10^{{{1}}}$".format(base, int(exponent))
#     else:
#         return float_str
RE=6.3781E8
ME=5.972E27

with open('../result/Structure.txt', encoding='utf8') as f:
  text = f.read().strip()
sections = text.split('\n\n')
datatype=[('Index', '<i8'), ('R', '<f8'), ('P', '<f8'), ('M', '<f8'), ('Density', '<f8'), ('T', '<f8'), ('Phase', 'S17')]
data = []

for i in range(len(sections)):
  data.append(np.genfromtxt(StringIO(sections[i]),delimiter='\t ', skip_header=1, dtype=datatype))

Mphase = []
phases = [data[-1][0]['Phase'].decode("utf-8")]

for i in range(1,len(data[-1])):
  if data[-1][i]['Phase'] != data[-1][i-1]['Phase']:
    phases.append(data[-1][i]['Phase'].decode("utf-8"))
    Mphase.append(data[-1][i-1]['M'])
  
my_dpi=96

fig,((ax1,ax2,ax3))=plt.subplots(3,1, sharex=True)

lns1=ax1.plot(data[-1][:]['M'],data[-1][:]['T'],'k-', label = 'Temperature')
ax1.set_ylabel(r'$T$ (K)',labelpad=0,fontsize=20)
ax1.tick_params('x',pad=8)

for i in range(len(Mphase)):
  ax1.axvline(x=Mphase[i])

ax2.plot(data[-1][:]['M'],data[-1][:]['R'],'k-',label = r'R')
ax2.set_ylabel(r'$R$ ($R_\oplus$)',labelpad=0 ,fontsize=20)

for i in range(len(Mphase)):
  ax2.axvline(x=Mphase[i])
  
ax3.plot(data[-1][:]['M'],data[-1][:]['Density'],'k-', label = 'Density')
ax3.set_ylabel(r'$\rho$ ($\mathrm{g~cm^{-3}}$)',labelpad=0,fontsize=20)
ax3.set_xlabel(r'Mass ($M_\oplus$)',labelpad=0,fontsize=20)
ax3.tick_params('x',pad=8)

for i in range(len(Mphase)):
  ax3.axvline(x=Mphase[i])
  ax3.text(Mphase[i]-0.04,12,phases[i],rotation=90,size=14)
ax3.text(Mphase[-1]+0.02,12,phases[-1],rotation=90,size=14)

for ax in fig.get_axes():
    ax.label_outer()

plt.tight_layout()
plt.subplots_adjust(hspace=.0)
plt.show()
#plt.savefig("profile_panels.pdf", dpi=my_dpi)
