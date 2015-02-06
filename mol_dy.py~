

import numpy as np
import random as rn
import math

# Structure

pos = [[],[],[]]
mom = []
F = [[],[],[]]
dt = 1
l = 1.
N=32.
gs = l*N**-(1./3.)

#print range(N)


'''
#Creating initial conditions
ICll = []
for z in range(3):
	ICl = []	
	for i in range(N):
		ICl.append(rn.normalvariate(0,1))
	ICll.append(ICl)
	


#plugging the initial conditions in the momenta lists
for z in range(3):
	mom.append(ICll[z])
	
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
		
print partl

#test with the "known" formula
'''
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

# Module to find 
