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
import correlationmodule

#flags
corrflag = 0

# independent parameters
dt = 0.004
N=256
lpnum = 3000
cutoff = lpnum/3.
density = 0.85
temp = 0.8
Ttarg = 1.8
Kb =  1
nbins = 100 #number of radial shells
alpha = 0.0001 #prerequisite for stability


# Loading initial conditions
mom = init_sys.init_mom(N, temp) 
pos, l, npdim = init_sys.init_pos(N, density) 
forces = init_sys.init_forc(N)
pot = init_sys.init_pot(N)
distances = init_sys.init_dist(N)
bin_vec_tot, finalbins = init_sys.init_bins(nbins,lpnum)
toten = init_sys.init_toten(lpnum)
presvirialtime = init_sys.init_presvirialtime(lpnum)
kenarray = init_sys.init_kenarray(lpnum)

# Iteration Verlet method
forces, pot, distances, presvirial = fm.calc_forces(pos,forces,pot,l,distances,[N])
formersummom = 0
presvirialtimelst = []
lamda=0
for lp in range(lpnum):
  mom = mom + forces*0.5*dt
  pos = (pos + mom*dt) % l           # % l means modulo of l, hence it adds/subtracts n*l untill 0<pos<l
  forces, pot, distances, presvirial = fm.calc_forces(pos,forces,pot,l,distances,[N])
  mom = mom + forces*0.5*dt
  Ken = np.sum(mom*mom*0.5, axis=1)
  kenarray[lp,:] = sum(Ken)
  toten[lp,:] = sum(Ken) + sum(pot)
  presvirialtime[lp,:] = presvirial
  if lp < (cutoff):
    if lp % 10 == 0:
      lamda = np.sqrt((Ttarg*3.*(N-1.)*Kb)/(np.sum(Ken)*2.))
      mom = mom*lamda
  elif corrflag == 1:
    plotflag = 0
    bin_vec_tot = bin_vec_tot + correlationmodule.cor(npdim,N,distances,nbins,finalbins,corrflag)
    



fig1 = plt.plot(kenarray)
fig2 = plt.plot(toten)
fig3 = plt.plot(toten-kenarray)
plt.show()

print np.mean(presvirialtime), 'pressure'

#print toten

# Calculate and plot correlationfunction
if corrflag == 1:
  finalbins = 2*bin_vec_tot/((lpnum - cutoff)*density*(N-1))
  bin_vec = correlationmodule.cor(npdim,N,distances,nbins,finalbins,corrflag)


























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
