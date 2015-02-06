#Initial momenta module

import random as rn

def im(N):
	ICll = []
	for i in range(N):
		ICL = []
		for dr in range(3):
			ICL.append(rn.normalvariate(0,10))
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

	return ICll
