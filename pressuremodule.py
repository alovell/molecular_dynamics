import numpy as np
import matplotlib.pyplot as plt


def pres(lpnum,cutoff,n,prestime):
  averagedp = np.zeros((lpnum-cutoff)/n, dtype=float)
  for p in range((lpnum-cutoff)/n):
    averagedp[p] = np.mean(prestime[n*p:n*(p+1)])



  #fig1 = plt.plot(kenarray)
  #fig2 = plt.plot(toten)
  #fig3 = plt.plot(toten-kenarray)
  fig4 = plt.plot(averagedp)
  plt.show()

  print min(prestime), max(prestime)

  mnp = np.mean(prestime)
  sdomp = np.std(prestime)/np.sqrt(len(prestime))
  mnavp = sum(averagedp)/len(averagedp)
  sdomavp = np.std(averagedp)/np.sqrt(len(averagedp))

  print mnp, sdomp, mnavp, sdomavp
  return
