import numpy as np
import matplotlib.pyplot as plt
import pylab

def cor(npdim,N,distances,nbins,bin_vec):
#find the correlation function of the system
  dmax = np.sqrt(5.)*npdim

  for i in range(N):
    for j in range(N):
      bin_num = int(distances[i][j]*nbins/dmax)
      bin_vec[bin_num] = bin_vec[bin_num] + 1

#normalize based on number of particles
  bin_vec = bin_vec/N
#normalize based on the volume of the radial shell
#keeping order R^2*dR -> V=4*PI*R^2*dR
  dR = dmax/nbins
  Rin = np.zeros((nbins), dtype=float)
  for bin_num in range(nbins):
    Rout = (bin_num + 1)*dR
    Rin[bin_num] = (bin_num)*dR
    volume = 4.*np.pi*(Rout**3 - Rin[bin_num]**3)/3.
    bin_vec[bin_num] = bin_vec[bin_num]/(4*np.pi*Rin[bin_num]**2*dR)#volume#((4./3.)*np.pi*dR**3)
    print Rout,'Rout',Rin[bin_num],'Rin'
  print bin_vec

#print histogram
  plt.bar(Rin,bin_vec,width=dR)
  plt.show()
  return bin_vec
