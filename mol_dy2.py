# Overview file

#import python classes
import numpy as np
import random as rn
import math
import matplotlib.pyplot as plt
import pylab
from mpl_toolkits.mplot3d import Axes3D



#import self produced classes
import forcemodule as fm
import init_sys


# independent parameters
dt = 0.0004
N=864
lpnum = 2000
density = 0.85
temp = 0.8
Ttarg = 300
Kb =  1


# Loading initial conditions
mom = init_sys.init_mom(N, temp) 
pos, l = init_sys.init_pos(N, density) 
forces = init_sys.init_forc(N)
pot = init_sys.init_pot(N)






# Iteration Verlet method

forces, pot = fm.calc_forces(pos,forces,pot,l,[N])
formersummom = 0
for lp in range(lpnum):
  mom = mom + forces*0.5*dt
  pos = (pos + mom*dt) % l           # % l means modulo of l, hence it adds/subtracts n*l untill 0<pos<l
  forces, pot = fm.calc_forces(pos,forces,pot,l,[N])
  mom = mom + forces*0.5*dt
  Ken = np.sum(mom*mom*0.5, axis=1)
  toten = sum(Ken) + sum(pot)
  if lp % 100 == 0:
    lamda = np.sqrt((Ttarg*3.*(N-1.)*Kb)/(np.sum(Ken)*2.))
#  print type(lamda), len(lamda), type(forces), type(mom), type(0.5*dt)
    mom = mom*lamda
    print lp, 'lp',lamda, 'lambda', np.sum(Ken), 'Ken'
  
  
 



#  print sum(Ken), sum(pot)[0], toten[0]#, np.sum(mom)
# Plotting the positions



'''
  fig = pylab.figure()
  ax = Axes3D(fig) 
  ax.scatter(pos[:,0],pos[:,1],pos[:,2],c='b')
  ax.set_xlabel('X Label')
  ax.set_ylabel('Y Label')
  ax.set_zlabel('Z Label')
  plt.show()
 '''
