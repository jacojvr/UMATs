import numpy as np
from os import system
system('f2py --fcompiler=gfortran -c ./response1d_isotropic.f'+
  ' ./simpleshear.f ./linearharden.f -m isostress')
from isostress import risotropic
import matplotlib.pyplot as plt
plt.close('all')
# Array containing the material properties
props = np.array([10000., 0.3, 50.,100.])
# Strain increments from 0 to 0.9
strains0 =np.linspace(0,0.9,51)
# Cycle
strains = np.r_[strains0,strains0[-2::-1],-strains0[1:],
  -strains0[-2::-1],strains0[1:]]
# Initial state variables
statev = np.zeros(1)
# Initial stress and stress array 
stress,stressarr = 0.,[0.]
for i in range(strains.size-1):
 dstrain = strains[i+1]-strains[i]
 stress,statev = risotropic(stress,dstrain,1.,1.,statev,props)
 stressarr += [stress]
# Create Figure
plt.figure(1)
plt.plot(strains,np.array(stressarr),'k-s',ms=2)
plt.grid(True)
plt.xlabel(r"$\varepsilon$",fontsize=14.)
plt.ylabel(r"$\sigma$",fontsize=14.)
plt.title('One Dimensional Isotropic Incremental Plasticity')
plt.show()