import numpy as np
from os import system
system('f2py --fcompiler=gfortran -c ./response1d_combined.f'+
 ' ./simpleshear.f ./combinedlinear.f -m combstress')
from combstress import rcombined
import matplotlib.pyplot as plt
plt.close('all')
# List of labels
labels = ['Isotropic','Kinematic','Combined']
styles = ['k-.','k--','k-']
# Array containing the material properties
props = np.array([[10000., 0.3, 50.,100.,0],
                  [10000., 0.3, 50.,0.,100],
                  [10000., 0.3, 50.,50.,50]])
# Strain increments from 0 to 0.9
strains0 =np.linspace(0,0.9,51)
# Cycle
strains = np.r_[strains0,strains0[-2::-1],-strains0[1:],
  -strains0[-2::-1],strains0[1:]]
# cycle through tests:
plt.figure(1)
for lbl in range(3):
  # Initial state variables
  statev = np.zeros(7)
  # Initial stress and stress array 
  stress,stressarr = 0.,[0.]
  for i in range(strains.size-1):
    dstrain = strains[i+1]-strains[i]
    stress,statev = rcombined(stress,dstrain,1.,1.,statev,props[lbl])
    stressarr += [stress]
  plt.plot(strains,np.array(stressarr),styles[lbl],label=labels[lbl])
plt.grid(True)
plt.xlabel(r"$\varepsilon$",fontsize=16.)
plt.ylabel(r"$\sigma$",fontsize=16.)
plt.title('One Dimensional Incremental Plasticity')
plt.legend(loc='upper left')
plt.show()