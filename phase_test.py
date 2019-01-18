import numpy as np
from matplotlib import pyplot as plt


from bunch import Bunch
from orbit.teapot import teapot
from orbit.teapot import TEAPOT_MATRIX_Lattice
from orbit.bunch_generators import TwissContainer, GaussDist2D

from orbit.orbit_correction import orbit

from orbit.errors import AddErrorNode
from orbit.errors import AddErrorSet

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

# get lattice function from transfer matrices
matrix_lattice = TEAPOT_MATRIX_Lattice(teapot_latt,b)
(muX, arrPosAlphaX, arrPosBetaX) = matrix_lattice.getRingTwissDataX()
(muY, arrPosAlphaY, arrPosBetaY) = matrix_lattice.getRingTwissDataY()


#----------------------------------Bunch-Distribution------------------------------------------------
# machted beam
emittance_x, emittance_y = 35*10**(-6),15*10**(-6) # this number I get just to check the idea
twissX = TwissContainer(arrPosAlphaX[0][1], arrPosBetaX[0][1], emittance_x)
twissY = TwissContainer(arrPosAlphaY[0][1], arrPosBetaY[0][1], emittance_y)

#----------------------------------------------------------------------------------------------------
dist = GaussDist2D(twissX,twissY)
x,px,y,py = dist.getCoordinates()
b.addParticle(x,px,y,py,0,0)
#----------------------------------------------------------------------------------------------------



turn = 512
x_turn, px_turn = list(), list()

for k in range(turn):
    teapot_latt.trackBunch(b)
    x_turn.append(b.x(0))
    px_turn.append(b.px(0))

plt.figure()
plt.scatter(x_turn,px_turn,s=5,color = 'blue', label = 'ideal orbit')
plt.grid(True, ls = '--')
plt.xlabel('x, m')
plt.ylabel('$p_x$, rad')
#plt.savefig('phase_space_x_ideal1.pdf')




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
paramsDict["dx"]       = 0.008
paramsDict["dy"]       = 0.008

setDict["elementtype"] = "quad"
setDict["ringline"] = "ring"

ESet  = AddErrorSet(teapot_latt, positioni, positionf, setDict, paramsDict, seed_value=50)
#ESet  = AddErrorSet(lattice, positioni, positionf, setDict, paramsDict) # Random

OrbitX, OrbitY = orbit(teapot_latt,b).get_orbit() 
x0, px0 = OrbitX[0][1], OrbitX[0][2]
y0, py0 = OrbitY[0][1], OrbitY[0][2]
#---------------------------------------------ORBIT ERRORS-------------------------------------------

b.deleteAllParticles()
b.addParticle(x0+x,px0+px,y0+y,py0+py,0,0)

x_turn1, px_turn1, y_turn1, py_turn1 = list(), list(), list(), list()

for k in range(turn):
    teapot_latt.trackBunch(b)
    x_turn1.append(b.x(0))
    px_turn1.append(b.px(0))
    y_turn1.append(b.y(0))
    py_turn1.append(b.py(0))

x00 = np.mean(x_turn1)
px00 = np.mean(px_turn1)
y00 = np.mean(y_turn1)
py00 = np.mean(py_turn1)

#----------------------------------------------------------------------------------------------------

b.deleteAllParticles()
b.addParticle(x0,px0,y0,py0,0,0)

x_turn2, px_turn2 = list(), list()

for k in range(turn):
    teapot_latt.trackBunch(b)
    x_turn2.append(b.x(0))
    px_turn2.append(b.px(0))


plt.scatter(x_turn2,px_turn2,s=5, color = 'black', label = 'start from 0')

b.deleteAllParticles()
b.addParticle(x00,px00,y00,py00,0,0)

x_turn3, px_turn3 = list(), list()

for k in range(turn):
    teapot_latt.trackBunch(b)
    x_turn3.append(b.x(0))
    px_turn3.append(b.px(0))


plt.scatter(x_turn3,px_turn3,s=5, color = 'grey', label = 'start from 00')
#----------------------------------------------------------------------------------------------------

plt.scatter(x_turn1,px_turn1,s=5, color = 'red', label = 'distorted orbit')
plt.scatter([x0],[px0], marker = 'x', color = 'blue', label = 'standart')
plt.scatter([x00],[px00], marker = 'x', color = 'red', label = 'averaging')


plt.grid(True, ls = '--')
plt.xlabel('x, m')
plt.ylabel('$p_x$, rad')
#plt.legend(loc="lower right")
plt.legend(loc="best")
plt.title('One particle in during 1024 consequent turns in ideal and distorted lattices')
plt.savefig('phase_space.pdf')
plt.show()













