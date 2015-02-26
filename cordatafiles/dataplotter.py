import numpy as np
import matplotlib.pyplot as plt

Rin = np.loadtxt("cordataT0.5rho1.2N864Rin.txt")
finalbins = np.loadtxt("cordataT0.5rho1.2N864finbin.txt")
prestime = np.loadtxt("presdataT3rho0.3N864n100lpnum10000prestime.txt")
averagedp = np.loadtxt("presdataT3rho0.3N864n100lpnum10000averagedp.txt")

#values for the different parameters, need to be equal to those mentioned in the filenames.
Ttarg = 0.5
density = 1.2
n = 100 #width of pressure block

  
def corplotter(Rin,finalbins,Ttarg,density):
  plt.plot(Rin,finalbins)#,width=dR)
  plt.xlim([0,6])
  plt.ylabel('g(r)')
  plt.xlabel('r')
  plt.title('Correlationfunction')
  plt.text(4.2,max(finalbins)-0.1,r'$T$=%s, $\rho$=%s'%(Ttarg,density), fontsize=18) # x and y values for position of text are choosen for 864 particles
  plt.show()
  #plt.savefig('correlationfunctionT%sRho%s.jpg'%(Ttarg,density))
  return
  
  
def presplotter(prestime,averagedp,n):
  fig = plt.figure()
  ax1 = fig.add_subplot(2,1,1)
  ax2 = fig.add_subplot(1,1,1)
  fig.subplots_adjust(hspace=.35)
  ax1.plot(range(len(prestime)),prestime)
  ax2.plot(range(len(averagedp)),averagedp)
  ax1.set_ylabel('pressure')
  ax1.set_xlabel('time')
  ax2.set_ylabel('pressure')
  ax2.set_xlabel('time')
  ax1.set_title('pressure evolution')
  ax2.set_title('pressure evolution with pressure blocks')
  ax2.text(1,2.715,'blocksize = %s'%(n),fontsize=12)
  plt.show()
  return
  
#presplotter(prestime,averagedp,n)


# calculating the errors of the pressure and energy
mnp = np.mean(prestime)
sdp = np.std(prestime)
sdpb = np.std(averagedp)
sdompb =  np.std(averagedp)/np.sqrt(len(averagedp))

print 'mnp=',mnp,'sdp=',sdp,'sdpb=',sdpb,'sdompb=',sdompb 


