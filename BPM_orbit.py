import numpy as np
from matplotlib import pyplot as plt

from orbit.teapot import teapot
from orbit.teapot import TEAPOT_MATRIX_Lattice

from bunch import Bunch, BunchTwissAnalysis
from orbit.bunch_generators import TwissContainer, GaussDist2D, KVDist2D

from orbit.orbit_correction import orbit, correction

from orbit.errors import AddErrorNode, AddErrorSet

#----------------------------------------------------------------------------------------------------

# get lattice from madx
teapot_latt = teapot.TEAPOT_Lattice()
teapot_latt.readMAD("sis100.lat","SIS100")

#---------------------------------------------Bunch init---------------------------------------------
b = Bunch()

b.mass(0.93827231)
b.macroSize(1.0)

energy = 1 #Gev
b.getSyncParticle().kinEnergy(energy)
#----------------------------------------------------------------------------------------------------

#---------------------------------------------ORBIT ERRORS-------------------------------------------
# WE INTRODUCE MISALIGNEMENT IN THE QUADRUPOLES; dx, dy = HOR AND VER DISPLACEMENT OF QUADRUPOLES
setDict = {}
paramsDict = {}
positioni = 0.0
positionf = teapot_latt.getLength()
paramsDict["errtype"]  = "StraightError"
paramsDict["subtype"]  = "TransDisp"
paramsDict["sample"]      = "Uniform"
paramsDict["maximum"]        = 0.5
paramsDict["minimum"]       = 0.0
paramsDict["dx"]       = 0.001
paramsDict["dy"]       = 0.001

setDict["elementtype"] = "quad"
setDict["ringline"] = "ring"

ESet  = AddErrorSet(teapot_latt, positioni, positionf, setDict, paramsDict, seed_value=50)

paramsDict["errtype"]  = "FieldError"
paramsDict["subtype"]  = "QuadField"
paramsDict["sample"]      = "Uniform"
paramsDict["fracerr"]      = 0.05

setDict["elementtype"] = "quad"
setDict["ringline"] = "ring"

ESet  = AddErrorSet(teapot_latt, positioni, positionf, setDict, paramsDict, seed_value=50)
#---------------------------------------------ORBIT ERRORS-------------------------------------------

# get lattice function from transfer matrices
matrix_lattice = TEAPOT_MATRIX_Lattice(teapot_latt,b)
(muX, arrPosAlphaX, arrPosBetaX) = matrix_lattice.getRingTwissDataX()
(muY, arrPosAlphaY, arrPosBetaY) = matrix_lattice.getRingTwissDataY()
#----------------------------------Bunch-Distribusion------------------------------------------------
# machted beam
emittance_x, emittance_y = 35*10**(-6),15*10**(-6) # should be checked and corrected
twissX = TwissContainer(arrPosAlphaX[0][1], arrPosBetaX[0][1], emittance_x)
twissY = TwissContainer(arrPosAlphaY[0][1], arrPosBetaY[0][1], emittance_y)

#---------------------------------------------CALCULATE ORBIT------------------------------

OrbitX, OrbitY = orbit(teapot_latt,b).get_orbit() 


x0, px0 = OrbitX[0][1], OrbitX[0][2]
y0, py0 = OrbitY[0][1], OrbitY[0][2]
#----------------------------------------------------------------------------------------------------



n=10**(6) # num of particles
dist = KVDist2D(twissX,twissY)
for i in range(n):
    xi,pxi,yi,pyi = dist.getCoordinates()
    b.addParticle(x0+xi,px0+pxi,y0++yi,py0+pyi,0,0)

#----------------------------------------------------------------------------------------------------


nodes = teapot_latt.getNodes()
position = [teapot_latt.getNodePositionsDict()[node][0] for node in nodes if node.getType()=='monitor teapot'] # s coordinates of BPMs

bunchtwissanalysis = BunchTwissAnalysis()

#----------------------------------------------------------------------------------------------------
# the method calculates (<x>,<y>) of the beam distribution for all BPMs in the lattice
# the tracking takes place inside this method


def getAllBPMsInfo(beam, lattice_nodes):
    lst = list()
    for node in lattice_nodes:
        node.trackBunch(beam)
        if node.getType()=='monitor teapot':
            bunchtwissanalysis.analyzeBunch(beam)
            x = bunchtwissanalysis.getAverage(0)
            y = bunchtwissanalysis.getAverage(2)
            lst.append([x,y])
            
    return lst

#----------------------------------------------------------------------------------------------------

# in order to get beam oscillations, one can repeat getALLBPMsInfo method n_turn times

n_turn = 256
arr = [getAllBPMsInfo(b,nodes) for i in range(n_turn)]

        
lst = np.transpose(arr)
x_full = lst[0]
y_full = lst[1]



x_mean = [np.mean(elem) for elem in x_full]
y_mean = [np.mean(elem) for elem in y_full]

mx = [elem-np.mean(elem) for elem in x_full]
my = [elem-np.mean(elem) for elem in y_full]

U, s, V = np.linalg.svd(np.transpose(mx), full_matrices=True)
betaX = pow(s[0],2)*(V[0]**2)+pow(s[1],2)*(V[1]**2)
U, s, V = np.linalg.svd(np.transpose(my), full_matrices=True)
betaY = pow(s[0],2)*(V[0]**2)+pow(s[1],2)*(V[1]**2)

beta_orbit_x = np.transpose(arrPosBetaX)
beta_orbit_y = np.transpose(arrPosBetaY)

kappaX = beta_orbit_x[1][5]/betaX[0]
kappaY = beta_orbit_y[1][5]/betaY[0]

# ORBIT
x = []
y = []
s = []
for i in xrange(len(OrbitX)):
	s.append(OrbitX[i][0])
	x.append(OrbitX[i][1])
	y.append(OrbitY[i][1])

#----------------------------------------------------------------------------------------------------
#Comparison between pyORBIT and multi-turn tracking methods


plt.figure()
plt.plot(beta_orbit_x[0],beta_orbit_x[1],'-x', label="X-beta pyORBIT")
plt.scatter(position,kappaX*betaX, color = 'red', label="X-beta MIA")
plt.legend(loc="lower right")
plt.grid(True, ls = '--')
plt.xlim([0,100])
#plt.savefig('beta.pdf')

plt.figure()
plt.plot(beta_orbit_y[0],beta_orbit_y[1],'b-x', label="Y-beta pyORBIT")
plt.scatter(position,kappaY*betaY, color = 'blue', label="Y-beta MIA")
plt.legend(loc="lower right")
plt.grid(True, ls = '--')
plt.xlim([0,100])
#plt.savefig('beta.pdf')


plt.figure()
plt.plot(s,x,'r-x', label="x orbit, standart pyORBIT")
plt.scatter(position,x_mean, label="x at BPMs, multi-turn tracking method")
plt.legend(loc="lower right")
plt.grid(True, ls = '--')
#plt.ylim([min(x),max(x)])
plt.title('X-plane. {} particles, {} turns'.format(n,n_turn))
#plt.title('X-plane. Drift length = 4 instead of 3.062')
#plt.savefig('X-orbit_drift.pdf')

plt.figure()
plt.plot(s,y,'g-x', label="y, standart pyORBIT")
plt.scatter(position,y_mean,label="x at BPMs, multi-turn tracking method")
plt.legend(loc="lower right")
plt.grid(True, ls = '--')
#plt.ylim([min(y),max(y)])
plt.title('Y-plane. {} particles, {} turns'.format(n,n_turn))
#plt.title('Y-plane. Drift length = 4 instead of 3.062')
#plt.savefig('Y-orbit_drift.pdf')

plt.show()








