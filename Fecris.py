import numpy as np

def gaussian_2d(x, y, amplitude, x_mean, y_mean, x_stddev, y_stddev, theta=0):
    """
    Calculates a 2D Gaussian function.

    Args:
        x (numpy.ndarray): X-coordinates of the grid.
        y (numpy.ndarray): Y-coordinates of the grid.
        amplitude (float): The peak value of the Gaussian.
        x_mean (float): The mean of the Gaussian in the x-direction.
        y_mean (float): The mean of the Gaussian in the y-direction.
        x_stddev (float): The standard deviation in the x-direction (before rotation).
        y_stddev (float): The standard deviation in the y-direction (before rotation).
        theta (float, optional): The rotation angle in radians (counterclockwise). Defaults to 0.

    Returns:
        numpy.ndarray: A 2D array representing the Gaussian function.
    """
    # Rotate the coordinates if theta is not zero
    x_prime = (x - x_mean) * np.cos(theta) - (y - y_mean) * np.sin(theta)
    y_prime = (x - x_mean) * np.sin(theta) + (y - y_mean) * np.cos(theta)

    # Calculate the exponent term
    exponent = -0.5 * ((x_prime**2 / x_stddev**2) + (y_prime**2 / y_stddev**2))

    # Calculate the 2D Gaussian
    return amplitude * np.exp(exponent)



class FECRIS2D:
    '''Objects simulating the response of an ECR ion source with 2 input parameters. 2 2D gaussian are definde in the init
    for the mean and stdev of the beram current statistics as function of the input.
    Methods:
    SetState(input) will force the ion source in the state of the given input
    no return
    Transition(newState) can be called recursively to simulate a transition from the current state to the one of
    the input given in argument and return the measure along the transition
    read() will simulate a measure at the current ion source state
    '''
    def __init__(self, State):
        self.State = State
        self.BeamCurrP = [1., .5, .5, .3, .3] # Parameters of the 2D Gaussian simulating the beam current mean
        self.StabP = [.5, .6, .6, .06, .06] # Parameters of the 2D Gaussian simulating the beam current stdev
        # param order is [Amplitude, x_mean, y_mean, x_span, y_span]
        self.time = 0
        self.Bcurr = gaussian_2d(State[0], State[1], self.BeamCurrP[0], self.BeamCurrP[1], self.BeamCurrP[2], self.BeamCurrP[3], self.BeamCurrP[4]) + np.random.randn() * gaussian_2d(State[0], State[1], self.StabP[0], self.StabP[1], self.StabP[2], self.StabP[3], self.StabP[4])

    def Transition(self, NewState):
        self.time += 1. #time stem between two calls of the transition
        stretch = .02 # violence of the transition
        P = self.BeamCurrP
        S = self.StabP
        arrival = (gaussian_2d(NewState[0], NewState[1], P[0], P[1], P[2], P[3], P[4]))
        arrStab = (gaussian_2d(NewState[0], NewState[1], S[0], S[1], S[2], S[3], S[4]))
        departure = (gaussian_2d(self.State[0], self.State[1], P[0], P[1], P[2], P[3], P[4]))
        depStab = (gaussian_2d(self.State[0], self.State[1], S[0], S[1], S[2], S[3], S[4]))
        ret = departure + (arrival - departure)/2 * (np.arctan(stretch * (self.time - 150))/ np.pi *2 + 1. ) + np.random.randn() * ((np.arctan(stretch*(self.time - 150)) *2/np.pi +1 ) * (arrStab - depStab)/2 + depStab)
        return  ret

    def SetState(self, State):
        self.State = State
        self.time = 0.
    def read(self):
        self.Bcurr = gaussian_2d(self.State[0], self.State[1], self.BeamCurrP[0], self.BeamCurrP[1], self.BeamCurrP[2], self.BeamCurrP[3], self.BeamCurrP[4]) + np.random.randn() * gaussian_2d(self.State[0], self.State[1], self.StabP[0], self.StabP[1], self.StabP[2], self.StabP[3], self.StabP[4])
        return self.Bcurr


