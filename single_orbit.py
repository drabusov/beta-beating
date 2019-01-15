import numpy as np
from matplotlib import pyplot as plt

from orbit.teapot import teapot
from orbit.teapot import TEAPOT_MATRIX_Lattice

from bunch import Bunch, BunchTwissAnalysis
from orbit.bunch_generators import TwissContainer, GaussDist2D

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

energy = 1.0 #Gev
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
paramsDict["dx"]       = 0.0095
paramsDict["dy"]       = 0.0095

setDict["elementtype"] = "quad"
setDict["ringline"] = "ring"

ESet  = AddErrorSet(teapot_latt, positioni, positionf, setDict, paramsDict, seed_value=50)
#ESet  = AddErrorSet(teapot_latt, positioni, positionf, setDict, paramsDict) # Random
#---------------------------------------------ORBIT ERRORS-------------------------------------------

# get lattice function from transfer matrices
matrix_lattice = TEAPOT_MATRIX_Lattice(teapot_latt,b)
(muX, arrPosAlphaX, arrPosBetaX) = matrix_lattice.getRingTwissDataX()
(muY, arrPosAlphaY, arrPosBetaY) = matrix_lattice.getRingTwissDataY()

#----------------------------------Bunch-Distribusion------------------------------------------------
# machted beam
emittance_x, emittance_y = 10**(-8),10**(-8) # should be checked and corrected
twissX = TwissContainer(arrPosAlphaX[0][1], arrPosBetaX[0][1], emittance_x)
twissY = TwissContainer(arrPosAlphaY[0][1], arrPosBetaY[0][1], emittance_y)



# num of particles = 1
dist = GaussDist2D(twissX,twissY)
#x,px,y,py = dist.getCoordinates()
#x,px,y,py = -0.00114903585122, 0.00124391282321, 0.00310204490361, 0.000305760909693
x,px,y,py = 0.,0.,0.,0.
b.addParticle(x,px,y,py,0,0)

#----------------------------------------------------------------------------------------------------


nodes = teapot_latt.getNodes()
position = [teapot_latt.getNodePositionsDict()[node][0] for node in nodes if node.getType()=='monitor teapot'] # s coordinates of BPMs

bunchtwissanalysis = BunchTwissAnalysis()


#----------------------------------------------------------------------------------------------------

def getAllBPMsInfo(beam, lattice_nodes):
    lst = list()
    for node in lattice_nodes:
        node.trackBunch(beam)
        if node.getType()=='monitor teapot':
            x = b.x(0)
            y = b.y(0)
            lst.append([x,y,px])
    return lst

#----------------------------------------------------------------------------------------------------

# in order to get beam oscillations, one can repeat getALLBPMsInfo method n_turn times

n_turn = 1024

arr = list() 
phase_space = list()
for i in range(n_turn):
    arr.append(getAllBPMsInfo(b,nodes))
    x = b.x(0)
    px = b.px(0)
    y = b.y(0)
    py = b.py(0)
    phase_space.append([x,px,y,py])
        
lst = np.transpose(arr)
x_full = lst[0]
y_full = lst[1]

ps = np.transpose(phase_space)



x_mean = [np.mean(x) for x in x_full]
y_mean = [np.mean(y) for y in y_full]

#---------------------------------------------CALCULATE ORBIT------------------------------

OrbitX, OrbitY = orbit(teapot_latt,b).get_orbit() 

x = []
y = []
s = []
for i in xrange(len(OrbitX)):
	s.append(OrbitX[i][0])
	x.append(OrbitX[i][1])
	y.append(OrbitY[i][1])



#----------------------------------------------------------------------------------------------------


plt.figure()
plt.plot(s,x,'-x', label="x")
plt.scatter(position,x_mean, color = 'red')
plt.legend(loc="lower right")
plt.grid(True, ls = '--')
#plt.ylim([min(x_mean),max(x_mean)])
#plt.savefig('X-orbit_quad.pdf')

plt.figure()
plt.plot(s,y,'b-', label="y")
plt.scatter(position,y_mean)
plt.legend(loc="lower right")
plt.grid(True, ls = '--')
#plt.ylim([min(y_mean),max(y_mean)])
#plt.savefig('Y-orbit_quad.pdf')



plt.figure()
plt.scatter(ps[0],ps[1],s=5)
plt.grid(True, ls = '--')
plt.xlabel('x, m')
plt.ylabel('$p_x$, rad')
plt.title('One particle in during {} consequent turns. X-plane'.format(n_turn))


plt.figure()
plt.scatter(ps[2],ps[3],s=5)
plt.grid(True, ls = '--')
plt.xlabel('y, m')
plt.ylabel('$p_y$, rad')
plt.title('One particle in during {} consequent turns. Y-plane'.format(n_turn))

plt.show()








