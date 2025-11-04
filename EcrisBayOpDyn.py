import numpy as np
from scipy.integrate import quad
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import Matern
import scipy as sp

'''Optimizer for ecr ion source control (static version)
This library contains the necessary function to optimize  and ECR ion source given a number of input parameters
under stability constraints.

This adds new features to a standard Bayesian optimizer:
- Auto-tuning of alpha according to the expected noise stdev
- Considering stability constraint and avoid unstable areas for the search using a second gaussian process

The object (Optimizer) uses the bayesian optimization framework of scikit-learn and needs:
Klen: list containing the kernel length along all the dimensions of parameters
AlphaVect: list of test values for the auto-tuning of the smoothness of the Gaussian process 
exBias: the exploitation vs exploration bias for the acquisition function typically 2.5
expNoise: the expected stdev of the noise
risk: the admissible risk to set the ecris in a potential unstable area
threshold: the value of the constraint over which the system is considered unstable

The function (NextPointQuery) returns the parameters for the next settings according to the bayesian optimizer
X: ndarray with the series of previous settings evaluated Dim0 is time and Dim1 is parameters
Y: 1darray with the series of output associated for the objective function (beam current mean) len = Dim0
S: 1darray with the series of output associated for the control function (Stability measure) len = Dim0
OptPar: Optimizer object
limits: 2darray with the boundaries of the search along each parameters np.asarray([[x0min,x1min,...],[x0max,x1max,...]])
'''

def Gaussian(x, m, s) :
    return 1 / (np.sqrt(2 * np.pi) * s) * np.exp(-0.5 * ((x - m) / (s)) ** 2)

def InvLCB(Param, gpr, kappa, gprStab, risk, StabLimit = 0.05, Dynamic = False) :
    if Dynamic:
        meanS, stdS = gprStab.predict([0, Param], return_std=True)
    else:
        meanS, stdS = gprStab.predict([Param], return_std=True)
    Sum = quad(Gaussian, StabLimit, np.inf, args=(meanS, stdS))
    if Sum[0] <= risk :
        if Dynamic:
            mean, std = gpr.predict([0, Param], return_std=True)
        else:
            mean, std = gpr.predict([Param], return_std=True)
        return -(mean + kappa * std)
    else :
        return np.inf

class Optimizer:
    "Making an optimizer object to update and give next coordinates"
    def __init__(self, Klen, AlphaVect, expBias, expNoise =.01, risk =.2, threshold = 0.05, Dynamic = False):
        self.Klen = Klen
        self.AlphaVect = AlphaVect
        self.Alpha = 0
        self.expBias = expBias
        self.risk = risk
        self.threshold = threshold
        self.expNoise = expNoise
        self.Dynamic = Dynamic



def NextPointQuery(X, Y, S, OptPar, limits):

    # normalization
    if OptPar.Dynamic:
        X[:,0] = np.abs(X[:,0] - X[len(X)-1,0])# setting différential time with 0 for current time and >0 for past
        limits[1, 0] = max(X[:, 0])
    X = (X - limits[0]) / (limits[1] - limits[0]) # to test


    YdataW = Y / np.max(Y)
    StdDataW = S / np.max(Y)

    # optimizing alpha
    # Setting up Anisotropic Kernel function
    kernel = Matern(length_scale=OptPar.Klen, nu=2.5)
    # Initializing regressor
    sigmas = []
    for alphaParam in OptPar.AlphaVect :  # Testing different values for alpha VW Mod Sept 22
        gpr = GaussianProcessRegressor(kernel=kernel, random_state=0, alpha=alphaParam)
        # Fitting the data
        gpr.fit(X, YdataW)
        sigmas.append(np.std(gpr.predict(X) - YdataW))  # Extrtacting the residue stdev
    sigmas_array = np.array(sigmas)
    ind = np.argmin(np.abs(sigmas_array - OptPar.expNoise))  # determining the best alpha with the expected noise level (might do better later with inflextion point)
    OptPar.alpha = OptPar.AlphaVect[ind]


    gpr = GaussianProcessRegressor(kernel=kernel, random_state=0, alpha=OptPar.alpha)
    gprStab = GaussianProcessRegressor(kernel=kernel, random_state=0, alpha=OptPar.alpha)

    gpr.fit(X, YdataW)
    gprStab.fit(X, StdDataW / YdataW)
    XdataSpaceDimension = len(X[0])

    mins = np.zeros(10)
    sol = np.zeros([10, XdataSpaceDimension])
    bounds = []
    for k in range(XdataSpaceDimension):
        bounds.append([0, 1])
    if OptPar.Dynamic:
        bounds[0] = [0., 0.]


    for k in range(10) :
        InitParam = np.random.rand(XdataSpaceDimension)
        if OptPar.Dynamic:
            InitParam[0] = 0.

        res = sp.optimize.minimize(InvLCB, InitParam, args=(gpr, OptPar.expBias, gprStab, OptPar.risk), bounds=bounds,
                                   method='Nelder-Mead') # to test
        mins[k] = res.fun
        sol[k] = res.x

    GlobMin = np.argmin(mins)
    NPQ = sol[GlobMin]


    return (NPQ * (limits[1] - limits[0])) + limits[0]