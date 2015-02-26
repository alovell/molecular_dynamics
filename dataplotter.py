import numpy as np
import matplotlib.pyplot as plt

Rin = np.loadtxt("cordataT3rho0.3N32Rin.txt")
finalbins = np.loadtxt("cordataT3rho0.3N32finbin.txt")
Ttarg = 3
density = 0.3


plt.plot(Rin,finalbins)#,width=dR)
plt.xlim([0,max(Rin)*(5/8.)])
plt.ylabel('g(r)')
plt.xlabel('r')
plt.title('Correlationfunction')
plt.text(max(Rin)*(1/4.),max(finalbins)-0.1,r'$T$=%s, $\rho$=%s'%(Ttarg,density)) # x and y values for position of text are choosen for 864 particles
plt.show()
print max(Rin)
plt.savefig('correlationfunctionT%sRho%s.jpg'%(Ttarg,density))
