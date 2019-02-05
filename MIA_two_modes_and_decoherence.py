#!/usr/bin/env python
# coding: utf-8


import numpy as np
from matplotlib import pyplot as plt
from matplotlib import rc
import random


#--------------------------------------------------------------------------------------
rc('text', usetex=False)
rc('font', serif ='Times')
rc('font', size=16)
rc('axes', linewidth=0.5)
rc('lines', linewidth=1.15	)
#rc('figure', figsize=(8.3,5.2))
rc('grid', c='0.5', ls='--', lw=0.5)
#--------------------------------------------------------------------------------------



# lets' define M BPMs readings
M = 84 # Number of BPMs
L = 1000 # the ring length
ds = L/M # The distance between BPMs

qb = 18.84 # Betatron tune (fraction part)
qs = 0.003 # Synchrotron tune. Replace by the formula sqrt(alpha_p*q*U/2/pi/E)

#lets assume the beta and the dispersion are equal at all BPMs
beta = 20
D = 1

N = 1000 # length of beamhistory
turn = np.array(range(N)) # turn indexes

# initial conditions (randomized)
J0 = random.random() # betatron motion
phi0 = 2*np.pi*random.random() # from zero to 2pi

print('Amplitude of betatrone osc = {} m; initial phase = {} rad'.format(J0, phi0))

delta0 = random.random() # synchrotron motion
psi0 = 2*np.pi*random.random() # from zero to 2pi

print('Amplitude of sync osc = {} dE/E; initial phase = {} rad'.format(delta0, psi0))

# Let's construct B matrix
b = list()
for m in range(M):
    Ab = np.sqrt(2*J0*beta)
    phi = 2*np.pi*qb*m/M + phi0 
    x_beta = Ab*np.cos(2*np.pi*qb*turn + phi) # betatron motion

    As = D*delta0
    psi = 2*np.pi*qs*m/M + psi0 
    x_sync = As*np.cos(2*np.pi*qs*turn + psi) # sync motion
    
    kappa = 0.002
    decoherence = np.exp(-kappa*turn)
    
    noise = Ab/100*np.random.rand(N)
    
    x_bpm = (x_beta + x_sync)*decoherence + noise
    b.append(x_bpm)




plt.figure()
plt.plot(b[0], label = 'full signal')
#plt.plot(As*np.cos(2*np.pi*qs*turn + psi), color = 'red', label = 'synchrotron mode')
plt.legend(loc='upper right',frameon=True)
plt.title('Beamhistory, 1st BPM')
plt.grid(True)




U, s, V = np.linalg.svd(np.transpose(b), full_matrices=True)


plt.figure()
plt.scatter(range(M),s, marker='x', color = 'red')
plt.title('Singular values')
plt.grid(True)





u1 = np.transpose(U)

mode1 = s[0]*u1[0]+s[1]*u1[1]
mode2 = s[2]*u1[2]+s[3]*u1[3]


spec_beta = np.abs(np.fft.rfft(mode1))
spec_sync = np.abs(np.fft.rfft(mode2))
q = np.linspace(0,0.5,len(spec_beta))




plt.figure()
plt.plot(q,spec_beta, color = 'red', label = '1st mode, betatron motion')
plt.plot(q,spec_sync, color = 'blue', label = '2nd mode, sync motion')
plt.legend(loc='best',frameon=True)

plt.grid(True)
plt.title('Fourier spectrum of the 1st and 2nd temporal modes')





spec_signal = np.abs(np.fft.rfft(b[1]))

plt.figure()
plt.plot(q,spec_signal, color = 'green')
plt.title('Fourier spectrum of the BPM signal')
plt.grid(True)





fig, axs = plt.subplots(2, 2)
fig.subplots_adjust(hspace=0.4,wspace = 0.3)
axs[0,0].plot(u1[0])
axs[0,1].plot(u1[1])
axs[1,0].plot(u1[2])
axs[1,1].plot(u1[3])
axs[0,0].set_title('first mode (cos)')
axs[0,1].set_title('first mode (sin)')
axs[1,0].set_title('second mode (cos)')
axs[1,1].set_title('second mode (sin)')


plt.show()





