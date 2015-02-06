#tes 1t

import numpy as np
import random as rn
import math
import forcemodule as fm

# Structure

pos = [[],[],[]]
mom = []
F = [[],[],[]]
dt = 1
l = 1.
N=21
gs = l*N**-(1./3.)
lpnum = 100

#print range(N)


'''
#Creating initial conditions
ICll = []
for z in range(3):
	ICl = []	
	for i in range(N):
		ICl.append(rn.normalvariate(0,1))
		#print ICl
		ICll.append(ICl)
'''
ICll = []
for i in range(N):
	ICL = []
	for dr in range(3):
		ICL.append(rn.normalvariate(0,1))
	ICll.append(ICL)

trx,trY,trz = 0,0,0
for i in range(len(ICll)):
	trx = trx + ICll[i][0]
	trY = trY + ICll[i][1]
	trz = trz + ICll[i][2]

mnx = trx/len(ICll)
mny = trY/len(ICll)
mnz = trz/len(ICll)

for i in range(len(ICll)):
	ICll[i][0] = ICll[i][0] - mnx
	ICll[i][1] = ICll[i][1] - mny
	ICll[i][2] = ICll[i][2] - mnz

trx2,trY2,trz2 = 0,0,0
for i in range(len(ICll)):
	trx2 = trx2 + ICll[i][0]
	trY2 = trY2 + ICll[i][1]
	trz2 = trz2 + ICll[i][2]

#print trx2,trY2,trz2, ICll



'''
#plugging the initial conditions in the momenta lists
for z in range(3):
	mom.append(ICll[z])

print mom


for i in range(3):
	tr = np.mean(mom[i])
	mom[i] = [x - tr for x in mom[i]]
	#print np.mean(mom[i])
'''	

	
#plugging the initial positions for x


xnum = int(round(N**(1./3.),0))

xplb = []
for i in range(xnum):
	xp = (gs/2)+i*gs
	xplb.append(xp)
	
if xplb[-1] > (l-(gs/2)):
	del xplb[-1]

xplm = []
xplm = [x+(gs/2) for x in xplb]
if xplm[-1] > (l-(gs/2)):
	del xplm[-1]
	
#plugging the intial positions for y

gsy = (np.sqrt(3)/2)*gs

ynum = int(round(l/((np.sqrt(3)/2)*gs)))

yplb = []
for i in range(ynum):
	yp = (gs/2)+i*gsy
	yplb.append(yp)

	
if yplb[-1] > (l-(gs/2)):
	del yplb[-1]
	
yplm = []
yplm = [y+(gs/2) for y in yplb]
if yplm[-1] > (l-(gs/2)):
	del yplm[-1]
	
#plugging the intitial positions for z

zpl = yplb

#print yplb
#print zpl

#making position vectors

partl=[]

for idx,z in enumerate(zpl):
	if idx%2==0:
		ypl = yplb
	else:
		ypl = yplm
	#print ypl, 'ypl'
	
	for idy,y in enumerate(ypl):
		if idy%2==0:
			xpl = xplb
		else:
			xpl = xplm
		for i in xpl:
			partx = i
			party = y
			partz = z
			part = [partx,party,partz]
			partl.append(part)
		
#print partl,len(partl)

#Creating the F0 list

rempart = len(partl)

forces = []
for i in range(rempart):
	forces.append([0,0,0])



'''
#test with the "known" formula

part_t = []
for k in range(ynum):
	for j in range(ynum):
		for i in range(xnum):
			partx = (gs/2)*(2*i + ((j+k)%2)) +0.15
			party = (gs/2)*(np.sqrt(3)*(j+(1./3.)*(k%2))) +0.15
			partz = (gs/2)*((2.*np.sqrt(6)/3.)*k) +0.15
			part=[partx,party,partz]
			part_t.append(part)
			
print part_t

 '''
# Moving particles

 
#print len(partl), len(partl[0]),range(3),len(range(len(partl))),len(mom)
for lp in range(lpnum):
	print partl[0]
	forces = fm.calc_forces(partl,forces)
	for i in range(len(partl)):
		for z in range(3):
			partl[i][z] = partl[i][z] + ICll[i][z]*dt + (forces[i][z]/2)*dt**2
			ICll[i][z] = ICll[i][z] + forces[i][z]*dt	


