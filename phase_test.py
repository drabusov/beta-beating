import numpy as np
from matplotlib import pyplot as plt


from bunch import Bunch
from orbit.teapot import teapot
from orbit.teapot import TEAPOT_MATRIX_Lattice
from orbit.bunch_generators import TwissContainer, GaussDist2D


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
emittance_x, emittance_y = 10**(-5),10**(-5) # this number I get just to check the idea
twissX = TwissContainer(arrPosAlphaX[0][1], arrPosBetaX[0][1], emittance_x)
twissY = TwissContainer(arrPosAlphaY[0][1], arrPosBetaY[0][1], emittance_y)



#----------------------------------------------------------------------------------------------------
dist = GaussDist2D(twissX,twissY)
x,px,y,py = dist.getCoordinates()
b.addParticle(x,px,y,py,0,0)
#----------------------------------------------------------------------------------------------------



turn = 1024
x_turn, px_turn = list(), list()

for k in range(turn):
    teapot_latt.trackBunch(b)
    x_turn.append(b.x(0))
    px_turn.append(b.px(0))

plt.figure()
plt.scatter(x_turn,px_turn,s=5)
plt.grid(True, ls = '--')
plt.xlabel('x, m')
plt.ylabel('$p_x$, rad')
plt.title('One particle in during 1024 consequent turns. Ideal orbit')
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
paramsDict["dx"]       = 0.006
paramsDict["dy"]       = 0.006

setDict["elementtype"] = "quad"
setDict["ringline"] = "ring"

ESet  = AddErrorSet(teapot_latt, positioni, positionf, setDict, paramsDict, seed_value=50)
#ESet  = AddErrorSet(lattice, positioni, positionf, setDict, paramsDict) # Random
#---------------------------------------------ORBIT ERRORS-------------------------------------------



x_turn, px_turn = list(), list()

for k in range(turn):
    teapot_latt.trackBunch(b)
    x_turn.append(b.x(0))
    px_turn.append(b.px(0))

plt.figure()
plt.scatter(x_turn,px_turn,s=5)
plt.grid(True, ls = '--')
plt.xlabel('x, m')
plt.ylabel('$p_x$, rad')
plt.title('One particle in during 1024 consequent turns. Distorted orbit')
#plt.savefig('phase_space_x_dist1.pdf')
plt.show()













