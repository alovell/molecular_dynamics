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
corrflag = 2

# independent parameters
dt = 0.004
N=864
lpnum = 1000
cutoff = int(lpnum/5.)
density = 0.85
temp = 0.8
Ttarg = 1
Kb =  1
nbins = 100 #number of radial shells
n = 30 #number of timesteps per timeblock in pressure calculation


# Loading initial conditions
mom = init_sys.init_mom(N, temp) 
pos, l, npdim = init_sys.init_pos(N, density) 
forces = init_sys.init_forc(N)
distances = init_sys.init_dist(N)
bin_vec_tot, finalbins = init_sys.init_bins(nbins,lpnum)
toten = init_sys.init_toten(lpnum)
prestime = init_sys.init_presvirialtime(lpnum-cutoff)
kenarray = init_sys.init_kenarray(lpnum)

# Iteration Verlet method
  # Initializing the first values for the loop
pot = 0.0
forces, pot, distances, presvirial = fm.calc_forces(pos,forces,pot,l,distances,[N])
lamda=0
cnt = 0
for lp in range(lpnum):
  mom = mom + forces*0.5*dt
  pos = (pos + mom*dt) % l           # % l means modulo of l, hence it adds/subtracts n*l untill 0<pos<l
  forces, pot, distances, presvirial = fm.calc_forces(pos,forces,pot,l,distances,[N])
  mom = mom + forces*0.5*dt
  Ken = np.sum(mom*mom*0.5)
  kenarray[lp,:] = Ken
  toten[lp,:] = Ken + pot
    #print Ken, pot, Ken+pot
  if lp < (cutoff):
    if lp % 10 == 0:
      lamda = np.sqrt((Ttarg*3.*(N-1.)*Kb)/(np.sum(Ken)*2.))
      mom = mom*lamda
  elif corrflag == 1:
    corrflag = 0
    bin_vec_tot = bin_vec_tot + correlationmodule.cor(npdim,N,distances,nbins,finalbins,corrflag)
    corrflag = 1
  elif corrflag ==2:
    prestime[cnt,:] = density*((2*Ken/(3*N)) - presvirial[0]/(3*N)) #Here we assume that the integral term to calculate the pressure is 0.
    
    cnt = cnt + 1
#print prestime    
#print density*((2*Ken/(3*N)) - presvirial[0]/(3*N)), prestime[-1][0]
print density*(2*Ken/(3*N))
averagedp = np.zeros((lpnum-cutoff)/n, dtype=float)
for p in range((lpnum-cutoff)/n):
  #print presvirialtime[n*p:n*(p+1)], p
  averagedp[p] = np.mean(prestime[n*p:n*(p+1)])




#fig1 = plt.plot(kenarray)
#fig2 = plt.plot(toten)
#fig3 = plt.plot(toten-kenarray)
fig4 = plt.plot(averagedp)
plt.show()

mnp = np.mean(prestime)
sdomp = np.std(prestime)/np.sqrt(len(prestime))
mnavp = np.mean(averagedp)
sdomavp = np.std(averagedp)/np.sqrt(len(averagedp))

print mnp, sdomp, mnavp, sdomavp

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
