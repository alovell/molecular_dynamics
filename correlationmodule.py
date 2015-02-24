import numpy as np
import matplotlib.pyplot as plt
import pylab

def cor(npdim,N,distances,nbins,finalbins,plotflag):
#find the correlation function of the system
  dmax = np.sqrt(5.)*npdim
  bin_vec = np.zeros(nbins, dtype = float)

  for i in range(N):
    for j in range(i+1,N):
      bin_num = int(distances[i][j]*nbins/dmax)
      if bin_num < nbins/2.:
        bin_vec[bin_num] = bin_vec[bin_num] + 1
      
      


#normalize based on the volume of the radial shell
#keeping order R^2*dR -> V=4*PI*R^2*dR
  dR = dmax/nbins
  Rin = np.zeros((nbins), dtype=float)
  for bin_num in range(nbins):
    Rout = (bin_num + 2)*dR
    Rin[bin_num] = (bin_num + 1)*dR
    #Rin[0] = dR
    #volume = 4.*np.pi*(Rout**3 - Rin[bin_num]**3)/3.
    bin_vec[bin_num] = bin_vec[bin_num]/(4*np.pi*Rin[bin_num]**2*dR)
    #print bin_vec
 
  
  

  if plotflag == 1:
    fig2 = plt.bar(Rin,finalbins,width=dR)
    plt.show()


  
  return bin_vec
  

