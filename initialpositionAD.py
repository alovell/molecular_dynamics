# Creating intial positions along the AMY-DAAN algorithm

def ip(N,gs):

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

	rempart = len(partl)
	forces = []
	for i in range(rempart):
		forces.append([0,0,0])
		

	return partl,forces,ynum,xnum

