#Finding intial positions according to the PERRY algorithm

import numpy as np

def ip(N,a):
	npdim = int((N/4)**(1./3.))
	partl = []
	for k in range(npdim):
		for j in range(npdim):
			for i in range(npdim):
				partl.append([0+2*a*i,0+2*a*j,0+2*a*k])
				partl.append([a+2*a*i,a+2*a*j,0+2*a*k])
				partl.append([a+2*a*i,0+2*a*j,a+2*a*k])
				partl.append([0+2*a*i,a+2*a*j,a+2*a*k])
	
	#print partl[0]
	for i in range(len(partl)):
		for z in range(3):
			#print partl[i][z], 'a'
			partl[i][z] = partl[i][z] - npdim*a/2. 
			#print partl[i][z], 'b'
	#print partl[0]


	
				

	rempart = len(partl)
	forces = []
	for i in range(rempart):
		forces.append([0,0,0])

				
			
#print part_t

	return partl,forces
