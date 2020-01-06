#------------------------------------------------
#pyORBIT error and correction example
#------------------------------------------------

import sys
import math
import matplotlib.pyplot as plt
import numpy as np

import orbit_mpi

from orbit.teapot import teapot
from orbit.teapot import TEAPOT_Lattice

from orbit.lattice import AccLattice, AccNode, AccActionsContainer
from bunch import Bunch

from orbit.matching import Optics, EnvelopeSolver

from orbit.teapot import TEAPOT_MATRIX_Lattice

from orbit.utils.orbit_mpi_utils import bunch_pyorbit_to_orbit

#---------------------------------------------Bunch init---------------------------------------------
print "Start."

b = Bunch()
energy = 1
syncPart=b.getSyncParticle()
syncPart.kinEnergy(energy)

#---------------------------------------------Make a Teapot Lattice----------------------------------
'''
print "Generate Lattice."
lattice = TEAPOT_Lattice("lattice")
#lattice.readMAD("/home/laptop/Documents/Scripts/FODO_thick.lat","fodo")
#lattice.readMAD("FODO_thick.lat","fodo")
lattice.readMAD("CRYRING_part.lat","cryring")
lattice.setUseRealCharge(useCharge = 1)
'''
#---------------------------------------------Make a Teapot Lattice----------------------------------

print "Generate Lattice."
latt = TEAPOT_Lattice("lattice")
#latt.readMADX("/home/laptop/Documents/Scripts/FODO_thick.seq","fodo")
#latt.readMADX("FODO_thick.seq","fodo")
#latt.readMADX("cryring.madx","cryring")
#latt.readMADX("sis100_full_thin_fix.seq","sis100ring")
latt.readMADX("sis100.madx","sis100ring")
latt.setUseRealCharge(useCharge = 1)

#---------------------------------------------SPLIT LONG ELEMENTS--------------------------------------

d = latt.getNodePositionsDict()

for node in latt.getNodes():
	L =node.getLength()
	val = d[node]
	print node.getName(),val,L, val[0]/2+val[-1]/2, val[0]+L/2
	node.setnParts(10)

#------------------------------------------------------------------------------------------------------
'''
dict1 =lattice.getNodePositionsDict()
dict2 =latt.getNodePositionsDict()
print(dict1)
print("SEQ vs LINE")
print(dict2)

#------------------------------------------------------------------------------------------------------
for node in lattice.getNodes():
	if node.getType()=="quad teapot":
		print(node.getParam("kq"))

print("SEQ vs LINE")

for node in latt.getNodes():
	if node.getType()=="quad teapot":
		print(node.getParam("kq"))

#------------------------------------------------------------------------------------------------------

matrix_lattice = TEAPOT_MATRIX_Lattice(lattice,b)

TwissDataX0,TwissDataY0 = matrix_lattice.getRingTwissDataX(),matrix_lattice.getRingTwissDataY()

beta_x0 = np.transpose(TwissDataX0[-1])
beta_y0 = np.transpose(TwissDataY0[-1])
'''
#------------------------------------------------------------------------------------------------------

matrix_latt = TEAPOT_MATRIX_Lattice(latt,b)

TwissDataX,TwissDataY = matrix_latt.getRingTwissDataX(),matrix_latt.getRingTwissDataY()

beta_x = np.transpose(TwissDataX[-1])
beta_y = np.transpose(TwissDataY[-1])

#------------------------------------------------------------------------------------------------------

def save_dict(data,filename):
	w= open(filename,"w")
	w.write(",".join(data.keys())+"\n")
	name=data.keys()[0]
	n = len(data[name])
	for i in range(n):
		out = []
		for key,arr in data.items():
			out.append("{}".format(arr[i]))
		w.write(",".join(out)+"\n")
	w.close()

d = {"s":beta_x[0],"betx":beta_x[1],"bety":beta_y[1]}
save_dict(d,"warmQuad.csv") 

#------------------------------------------------------------------------------------------------------
'''
plt.figure()
plt.plot(beta_x0[0],beta_x0[1])
plt.plot(beta_y0[0],beta_y0[1])

plt.plot(beta_x[0],beta_x[1], ls="--")
plt.plot(beta_y[0],beta_y[1], ls="--")
'''
#------------------------------------------------------------------------------------------------------

fig, ax1 = plt.subplots()

ax1.plot(beta_x[0],beta_x[1])
ax1.plot(beta_y[0],beta_y[1])

ax2 = ax1.twiny()

ax2.plot(beta_x[0],beta_x[1], ls="--")
ax2.plot(beta_y[0],beta_y[1], ls="--")

ax1.set_xlim([0,1083.6/6.0])
ax2.set_xlim([1083.6/6.0,1083.6/3.0])

#--------------------------------------------------------------------------------
plt.show()
#--------------------------------------------------------------------------------


