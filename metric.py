# optimization of the sis100 structure with fielderr implemented

import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import least_squares as ls

from readMAD import *

#----------------------------------------------------------------------------
#filename = 'sis100.madx'
filename = 'cryring.madx'

def set_quads(theta, flag, r): 
    print('correctors are: {}'.format(theta))
    print('Field error is {}'.format(flag ==1))
    print('Isolated quadrupole {}'.format(r))
    g = open('optics_test.str', 'w')
    for i,v in enumerate(theta):
        g.write('{} := {};\n'.format(corr[i],v))

    g.write('{} := {};\n'.format('qv0',r)) # additional quads with error field 
    
    g.write('{} := {};\n'.format('flg',flag)) # in case of zero fiel errors has to be 0, change to 1 in another case    
    g.close()
    os.system("/home/dmitrii/MADX/madx < {} > ./out.dat".format(filename))
    print('twiss file with the {} lattice computed'.format('ideal' if flag ==0 else 'distorted'))

def cos_func(phi, tune):
    if phi <0:
        phi +=tune
    return np.cos(-2*np.pi*tune +2*np.pi*2*phi)

# old version
#def get_beat(set_corr,beta_b, beta_c, mu_beta, mu_corr, q, l):
#    print('beta at bpm is {} position of bpm {} tune {} len {}'.format(beta_b, mu_beta, q, l))
#    tmp = [k_i*beta_c[i]*cos_func(mu_beta-mu_corr[i],q) for i,k_i in enumerate(set_corr)]
#    return -np.sum(tmp)/2./np.sin(2*np.pi*q)


def get_beat(set_corr,beta_b, beta_c, mu_beta, mu_corr, q, l):
    
    tmp=0
    for i,k_i in enumerate(set_corr):

        phi0 = 2*np.pi*q

        phi = 2*np.pi*(mu_beta-mu_corr[i])

        if phi<0:
            phi+=phi0



        tmp += k_i*beta_c[i]*np.cos(2*phi-phi0)

    return -tmp/2./np.sin(phi0)



def metric_counter(theta, betabpm, betacorr,mubpm, mucorr, Q1, L, beat_initial):

    beat = [get_beat(theta, b, betacorr, mubpm[j], mucorr, Q1, L) for j,b in enumerate(betabpm)] # count beating initiated by quads
    beat += beat_initial # take initial beating into the account

    beta_period = np.transpose(np.reshape(beat, (6,14)))
    metric = np.mean([np.std(b) for b in beta_period])

#    print('metric is {}    metric0 is {}'.format(metric, np.std(beat0)))


    return metric


def first_assumption(theta, betabpm, betacorr,mubpm, mucorr, Q1, L, beat_initial):

    beat = [get_beat(theta, b, betacorr, mubpm[j], mucorr, Q1, L) for j,b in enumerate(betabpm)] 
    beat += beat_initial

    metric = np.std(beat)    
#    print('metric is {}    metric0 is {}'.format(metric, np.std(beat0)))

    return metric

# this method probably will be tested with new ideas
def combined_metrics(theta, betabpm, betacorr,mubpm, mucorr, Q1, L, beat_initial):

    beat = [get_beat(theta, b, betacorr, mubpm[j], mucorr, Q1, L) for j,b in enumerate(betabpm)] # count beating initiated by quads
    beat += beat_initial # take initial beating into the account
    metric_d = np.std(beat)    # deviation

    beta_period = np.transpose(np.reshape(beat, (6,14)))
    metric_p = np.mean([np.std(b) for b in beta_period]) # periodicity

#    print('metric is {}    metric0 is {}'.format(metric, np.std(beat0)))


    return metric_d**2 + metric_p**2*0


def beat_counter(theta):
    print(theta)
    betabpm, betacorr,mubpm, mucorr, Q1, L = read_twiss('twiss.txt',SP2)
    beat = [get_beat(theta, b, betacorr, mubpm[j], mucorr, Q1, L) for j,b in enumerate(betabpm)]
    return beat


def optimize(vector, beta_ideal, n_steps, precision, opt_type):
    
    print('Start optimization')

    arg = read_twiss('twiss.txt',SP2)
    beat_initial = (np.array(arg[0]) - np.array(beta_ideal))/np.array(beta_ideal)
    arg += (beat_initial,)

    if opt_type == 0:
        metric0 = np.std(beat_initial) 
        print('Start from metric {}'.format(metric0))
        print('Try to make zero deviation of beta from ideal')

        vec = ls(first_assumption, vector, args = arg, max_nfev = n_steps, ftol = precision,xtol = precision)
    else:
        beat_periodic = np.transpose(np.reshape(beat_initial, (6,14)))
        metric0 = np.mean([np.std(b) for b in beat_periodic])
        print('Start from metric {}'.format(metric0))
        print('Make beta-function periodic again!')

        vec = ls(metric_counter, vector, args = arg, max_nfev = n_steps, ftol = precision)

    metric = np.mean(vec.fun)
    print('Optimum {}found'.format('' if vec.success== True else 'not'))
    print('Final metric is {}'.format(metric))
    print('Optimization rate is {}'.format(metric0/metric))
    set_quads(vec.x, 1, R)

    return vec.x

#----------------------------------------------------------------------------

SP2 = 48 # header of twiss file


#corr = ['q1','q2','q3','q4','q5','q6']
corr = ['q1','q2','q3','q4','q5','q6','q7','q8','q9','q10','q11','q12']

R = np.random.normal(0,0.05) # quad field violation
 
zeros = [0.0 for elem in corr] 
set_quads(zeros, 0, 0) # in order to get twiss file with the ideal lattice
betabpm0, betacorr0,mubpm0, mucorr0, Q0, L0 = read_twiss('twiss.txt',SP2) # read beta of the ideal lattice
s, beta0 = read_beta('twiss.txt',SP2)

# theta is the vector of parameters to optimize, theta0 is the first point to start
theta0 = [np.random.normal(0,0.001) for elem in corr] 
#theta0 = [np.random.normal(0,0.005) if elem =='q6' else 0.0 for elem in corr] # use only one corrector
#theta0 = zeros # simple test




'''
beat =beat_counter(theta0) # beating computed by analytic formula

set_quads(theta0, 1, R) # in order to get twiss file with the lattice with corr set errors
betabpm, betacorr, mubpm, mucorr, Q, L = read_twiss('twiss.txt',SP2) # read beta of the distorted lattice
beat0 = (np.array(betabpm) - np.array(betabpm0))/np.array(betabpm0) # beating computed by MADX


beatc = (np.array(betacorr) - np.array(betacorr0))/np.array(betacorr0) # beating corr computed by MADX
print(beatc)
print(1/np.tan(2*np.pi*Q0))
[print(2*beatc[i]/betacorr0[i]/v) for i,v in enumerate(theta0) if v]

plt.figure()
plt.plot(beat0)
plt.plot(beat)
#[plt.axvline(mu) for mu in mucorr]
'''

set_quads(zeros, 1, R) # in order to get twiss file with distorted lattice without correction
beta_dist = read_twiss('twiss.txt',SP2)[0] # read beta the distorted lattice




#theta1 = optimize(theta0, betabpm0, 100, 10**(-3), 0)
#theta2 = optimize(theta1, betabpm0, 100, 10**(-4), 1)
#theta3 = optimize(theta2, betabpm0, 100, 10**(-5), 0)
#theta4 = optimize(theta3, betabpm0, 2000, 10**(-6), 1)
theta4 = optimize(theta0, betabpm0, 2000, 10**(-8), 0)


betabpm, betacorr,mubpm, mucorr, Q1, L = read_twiss('twiss.txt',SP2)
s, beta = read_beta('twiss.txt',SP2)

plt.figure()
plt.scatter(mubpm, betabpm, marker = 'x', color = 'red')
plt.scatter(mubpm, betabpm0)
plt.scatter(mubpm, beta_dist, marker = '^', color = 'green')




beat = (np.array(beta) -np.array(beta0))/np.array(beta0)
idx = int(len(s)/6.)

fig, axs = plt.subplots(3, 2)
fig.subplots_adjust(hspace=0.4,wspace = 0.3)
axs[0,0].plot(s[0:idx],beat[0:idx])

#axs[0,0].scatter(np.array(mubpm[0:14])*L/Q1, (np.array(betabpm[0:14]) -np.array(betabpm0[0:14]))/np.array(betabpm0[0:14]))

axs[0,1].plot(s[idx:2*idx],beat[idx:2*idx])
axs[1,0].plot(s[2*idx:3*idx],beat[2*idx:3*idx])
axs[1,1].plot(s[3*idx:4*idx],beat[3*idx:4*idx])
axs[2,0].plot(s[4*idx:5*idx],beat[4*idx:5*idx])
axs[2,1].plot(s[5*idx:6*idx],beat[5*idx:6*idx])




plt.show()








