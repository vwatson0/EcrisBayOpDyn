

import numpy as np
import matplotlib.pyplot as plt
import EcrisBayOpDyn
import Fecris as venus
import GenKFlib as kfilt

import warnings
warnings.filterwarnings("ignore")

''' tests of the functionalities of the EcrisBayOpDyn and GenKFlib on synthetic 2D + time Fake ecris

 Create ecris and randomly search 3 points
 Start the optimizer and create time vector as column 0 of Settings
 While loop with KF and initial Buffer zone to monitor transition
 When transition over grab statistics and re-optimize
 
 The first value of Klen is the length of the kernel following time. In essence the value of previous queries on the current state
 of the system. Since the system moves slowly, this value should be high compared tto the others .1 -> 1. or 2. maybe.
 If the previous queries happened a long time ago, the user should reduce the first value of Klen to relate the lack
 of confidence in the information held by such measures.
 
 The computation of the Gaussian process can be expensive, when the number measures accumulate. If the system runs for 
 an extensive amount of time, the user should reduce the size of the Settings and associates Y and S sent to optimize
 By removing older queries especially in areas that have been revisited since.
 
'''
# parameters of the Kalman filter
cMeas = 1e4 # filter parameter the increase to filter more decrerasse to track faster
slopeMax  = .001 # slope used to determine if the system is settled

# parameters for the bayesian optimizer
Klen = [.1, .02, .02] # kernel length (1st value is time kernel depends on time span bigger time span -> smaller kernel)
alpha = [1e-6, 1e-5, 1e-4, 1e-3, 1e-2] # possible values for alpha
exBias = 2.6 # exploration exploitation bias
risk = .1 # acceptable rsik of evaluating an unstable area
thresh = .05 # stability threshold
expNoise = .01 # expected noise level
limits = np.asarray([[0,0,0],[1,1,1]]) # limits of the exploration space

Disp = True # Turns on display at every step
Dynamic = True
# lists for data storage and display
storeMeas = []
slope = []
FilteredI = []


Settings = np.random.rand(3, 2) # Settings where the Fecris will ber evaluated

myEcris = venus.FECRIS2D(Settings[0]) # creating the ion source object with the initial state
Y = [] # list for output (beam current)
S = [] # list for stability
tmp = []
# getting the fisrt statistics the old way
for k in range(100):
    tmp.append(myEcris.read())
Y.append(np.mean(tmp))
S.append(np.std(tmp))

angle = 0
imgNum = 0

# evaluating the first 3 random Settings
for k in range(len(Settings)-1):
    unsettled = True
    KF = kfilt.KFobject(np.atleast_1d(myEcris.read()), cMeas) # reinitialize Kalman filter for every new settings
    inftyLoop = 0
    while(unsettled): # going through the transition
        inftyLoop += 1
        measure = myEcris.Transition(Settings[k+1])# getting a transition measurement from the Fecris
        deltaT = np.random.randn() * .3 + .1 # forcing irregular sampling
        KF.EstimateState(np.atleast_1d(measure), deltaT) # Filtering the new measure
        if inftyLoop > 250: # after 20 call we start to evaluate the slope to see if the system has settled
            if np.prod(np.abs(slope[len(slope)-101::]) < slopeMax): # if the last 10 values of the slope meet the settling criterion
                unsettled = False
                Y.append(KF.X[0])# grabing the filtered output value
                S.append(KF.Sig[0]/KF.X[0]) # calculating the stability value
                myEcris.SetState(Settings[k+1]) # !!!! Important !!!! before moving to a new transition we force Fecris to new state
                print('settled at t=', inftyLoop) # how many call were needed for the system to settle

        if inftyLoop > 5000: # if the system never settles
            print('timeout')
            # grabing values but will pollute data set
            Y.append(KF.X[0])
            S.append(KF.Sig[0]/KF.X[0])
            unsettled = False
            myEcris.SetState(Settings[k + 1])
        # storing values for display
        storeMeas.append(measure)
        FilteredI.append(KF.X[0])
        slope.append(KF.X[1])

Settings = np.concatenate((np.atleast_2d(np.arange(len(Settings))), Settings.T), axis = 0).T # adding time stamps column 0
#print(Settings)
#Opt = EcrisBayOps.Optimizer(Klen, alpha, exBias, expNoise, risk, thresh) # building the optimizer
Opt = EcrisBayOpDyn.Optimizer(Klen, alpha, exBias, expNoise, risk, thresh, Dynamic = Dynamic) # building the optimizer

timeout = False
for k in range(35):
    if k == 25:
        Opt.expBias = 2.3
    if timeout:
        nextSet = np.random.rand(2)
        timeout = False
    else:
        #nextSet = EcrisBayOps.NextPointQuery(Settings, Y, S, Opt, limits) # getting next settings
        nextSet = EcrisBayOpDyn.NextPointQuery(Settings, Y, S, Opt, limits) # getting next settings
    # going through the same steps as followed during the random search
    unsettled = True
    KF = kfilt.KFobject(np.atleast_1d(myEcris.read()), cMeas)
    inftyLoop = 0
    while unsettled:
        inftyLoop += 1
        if Dynamic:
            measure = myEcris.Transition(nextSet[1::])
        else:
            measure = myEcris.Transition(nextSet)
        deltaT = np.random.randn() * .3 + .1
        KF.EstimateState(np.atleast_1d(measure), deltaT)
        if inftyLoop > 250 :
            if np.prod(np.abs(slope[len(slope)-101::]) < slopeMax):
                unsettled = False
                Y.append(KF.X[0])
                S.append(KF.Sig[0]/KF.X[0])
                Settings = np.concatenate((Settings, [nextSet]))
                if Dynamic:
                    myEcris.SetState(nextSet[1::])
                else:
                    myEcris.SetState(nextSet)
                print('settled at t=', inftyLoop)

        if inftyLoop > 50000 :
            print('timeout')
            timeout = True
            Settings = np.concatenate((Settings, [nextSet]))
            Y.append(KF.X[0])
            S.append(KF.Sig[0] / KF.X[0])
            unsettled = False
            if Dynamic :
                myEcris.SetState(nextSet[1 : :])
            else :
                myEcris.SetState(nextSet)
        Settings[len(Settings)-1, 0] = k + 3
        storeMeas.append(measure)
        FilteredI.append(KF.X[0])
        slope.append(KF.X[1])
    print('Fecris evaluated at', Settings[len(Settings) - 1])
if Disp:
    if 1:
    #for l in range(30):
        fig = plt.figure()
        ax = fig.add_subplot(122, projection='3d')

        inds = np.where(np.asarray(S[0:k+4]) > .05)[0]
        inds2 = np.where(np.asarray(S[0:k+4]) <= .05)[0]
        ax.scatter(.5, .5, 1., marker='*', color='black')
        p = ax.scatter(Settings[inds2, 1], Settings[inds2, 2], np.asarray(Y)[inds2])
        ax.scatter(Settings[inds2[len(inds2)-1], 1], Settings[inds2[len(inds2)-1], 2], np.asarray(Y)[inds2[len(inds2)-1]], s = 50, marker = '^', color = 'green')
        if len(inds):
            ax.scatter(Settings[inds, 1], Settings[inds, 2], np.asarray(Y)[inds], s = 40, marker = 'x', color = 'red')


        ax.set_xlabel(r'$\theta_1$')
        ax.set_ylabel(r'$\theta_2$')
        ax.set_zlabel(r'$\mathcal{I}$')
        ax.legend([ 'Target', 'Stable', 'last Eval','Unstable'])
        ax.set_xlim([0., 1.])
        ax.set_ylim([0., 1.])
        plt.title('Bay Opt on Ecris simulation')
        ax.view_init(elev=30, azim=angle)
        #cbar  = fig.colorbar(p, orientation = "horizontal")
        #cbar.set_clim(.0, .05)
        ax1 = fig.add_subplot(221)
        ax1.plot(storeMeas)
        ax1.plot(FilteredI)
        ax1.plot(1. * np.ones(len(storeMeas)), 'k--')
        if k >1:
            ax1.plot( np.max(np.asarray(Y)[inds2]) * np.ones(len(storeMeas)), 'g--')
        ax1.legend(['Output', 'Filtered', 'target', 'maxStable'])
        ax1.set_ylabel('Beam')
        ax2 = fig.add_subplot(223)
        ax2.plot(slope)
        ax2.plot(np.ones(len(slope)) * slopeMax, '--')
        ax2.plot(np.ones(len(slope)) * slopeMax, 'r--')
        ax2.plot(-np.ones(len(slope)) * slopeMax, 'r--')
        ax2.legend(['slope', 'limit'])
        ax2.set_ylabel('system settle')
        ax2.set_xlabel('time')
        fig.set_size_inches(12, 7)
        plt.show()
        # this part is just to make a video following the test
        #angle += 5
        #if angle >= 360:
        #    angle -= (360)
        #if imgNum < 10:
        #    plt.savefig('pics4vid/000'+str(imgNum)+'.png')
        #elif imgNum < 100:
        #    plt.savefig('pics4vid/00'+str(imgNum)+'.png')
        #elif imgNum < 1000:
        #    plt.savefig('pics4vid/0'+str(imgNum)+'.png')
        #else:
        #    plt.savefig('pics4vid/' + str(imgNum) + '.png')


        #imgNum+=1
        #plt.close()