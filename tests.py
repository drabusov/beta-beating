import numpy as np
from matplotlib import pyplot as plt


from bunch import Bunch
from orbit.teapot import teapot
from orbit.teapot import TEAPOT_MATRIX_Lattice
from orbit.bunch_generators import TwissContainer, GaussDist2D


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
#----------------------------------------------------------------------------------------------------

# to make sure that it is really the beta-function
aa= np.transpose(arrPosBetaX)
bb= np.transpose(arrPosBetaY)
plt.plot(aa[0],aa[1], label = 'X-beta(s)')
plt.plot(bb[0],bb[1], label = 'Y-beta(s)')
plt.gca()
plt.xlim([0,200])
plt.xlabel("s [m]")
plt.ylabel("beta [m]")
plt.grid(True, ls = '--')
plt.legend()
#plt.savefig('beta.pdf')
# looks like everything is OK


#----------------------------------------------------------------------------------------------------

#----------------------------------Bunch-Distribusion------------------------------------------------
# machted beam
emittance_x, emittance_y = 10**(-5),10**(-5) # this number I get just to check the idea
twissX = TwissContainer(arrPosAlphaX[0][1], arrPosBetaX[0][1], emittance_x)
twissY = TwissContainer(arrPosAlphaY[0][1], arrPosBetaY[0][1], emittance_y)


n =10000 # num of particles

dist = GaussDist2D(twissX,twissY)
for i in range(n):
    x,px,y,py = dist.getCoordinates()
    b.addParticle(x,px,y,py,0,0)


#----------------------------------------------------------------------------------------------------

#-------------------------------Plot initial phase space----------------------------------------------
y_space =[ (b.y(i), b.py(i)) for i in range(n)]
plt.figure()
plt.scatter(*np.transpose(y_space))
plt.xlabel('y, m')
plt.ylabel('$p_y$, rad')
plt.grid(True, ls = '--')
plt.title('Phase Space in (y, $p_y$) plane, N = 10000 particles, turn = 0')
#plt.savefig('phase_space_0.pdf')

#-----------------------------One turn transport-----------------------------------------------------
teapot_latt.trackBunch(b)
y_space_new =[ (b.y(i), b.py(i)) for i in range(n)]
#----------------------------------------------------------------------------------------------------



#-------------------------------Plot final phase space-----------------------------------------------
plt.figure()
plt.scatter(*np.transpose(y_space_new))
plt.grid(True, ls = '--')
plt.xlabel('y, m')
plt.ylabel('$p_y$, rad')
plt.title('Phase Space in (y, $p_y$) plane after one turn')
#plt.savefig('phase_space_1.pdf')


#-------------------------------Checking--------------------------------------------------------------

# the method calculates rms emittance (<x^2>*<x'^2>-<x*x'>^2)^(-2)
def emitt_calc(space):
    arr =[ [a[0]*a[0], a[1]*a[1], a[0]*a[1]] for a in space]
    arr1 = np.transpose(arr)
    e =np.mean(arr1[0])*np.mean(arr1[1]) - np.mean(arr1[2])**2 
    return np.sqrt(e)

print(emitt_calc(y_space))
print(emitt_calc(y_space_new))
#----------------------------------------------------------------------------------------------------


#-------------------------------Tracking of one particle --------------------------------------------
b.deleteAllParticles()
x,px,y,py = dist.getCoordinates()
b.addParticle(x,px,y,py,0,0)

turn = 8024
x_turn, px_turn = list(), list()

for k in range(turn):
    teapot_latt.trackBunch(b)
    x_turn.append(b.x(0))
    px_turn.append(b.px(0))

plt.figure()
plt.scatter(x_turn,px_turn)
plt.grid(True, ls = '--')
plt.xlabel('x, m')
plt.ylabel('$p_x$, rad')
plt.title('One particle in (x, $p_x$) phase plane during 8024 consequent turns')
#plt.savefig('phase_space_x.pdf')
plt.show()













