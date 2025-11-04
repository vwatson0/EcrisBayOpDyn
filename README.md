Bayesian dynamic optimization for non-invariant system under control constraint with Gaussian process smoothness auto-tunning 

This pakage includes:
- 2 libraries necessary to run a bayasian optimizer on an ECR ion source like VENUS (EcrisByOpDyn.py and GenKFlib.py)
- 1 library to simulate a 2dimensional input fake ion source (Fecris.py) // note: In this simulation the source is static even is considered dynamic during the test
- 1 test file runing the simulation and example of the use of the functionalities of EcrisByOpDyn.py and GenKFlib.py

  To use the kalman filter for time series tracking:
  
    import GenKFlib

  Create an object:

    myObj = GenKFlib.KFobject(FirstMeasure(vector), FilterParameter(scalar or vector))

        The filter parameter is the covariance of the noise in each dimension. If different depending on the dimension, send a vector of the same size as the measure

  To estimate recusively the stats of the time series:

    myObj.EstimateState(CurrentMeasure(vector), ElaspsedTimeSinceLastCall(scalar))

  To have the current statistics:
    MeanValue = myObj.X[0:len(CurrentMeasure)]
    StDev = myObj.Sig[0:len(CurrentMeasure)]
    slope = myObj.X[len(CurrentMeasure)::]


  To use the optimizer:

    import EcrisBayOpDyn

  Create an object:

    Opt = EcrisBayOpDyn.Optimizer(Klen, alpha, exBias, expNoise, risk, thresh, Dynamic)

  The parameters are:

  - Klen: kernel length (vect)
  - alpha: possible values for alpha (vect)
  - exBias: exploration exploitation bias (scalar)
  - risk: acceptable rsik of evaluating an unstable area (scalar)
  - thresh: stability threshold (scalar)
  - expNoise: expected noise level (scalar)
  - limits: limits of the exploration space (2xn vect)
  - Dynamic is a bool depending if the system is considered dynamic (time variant model) or static 
 
  Get the next settings for the system:

    nextSettings = EcrisBayOpDyn.NextPointQuery(PreviousSettings, Y, S, Opt, limits)

  The parameters are:

    - nextSettings: Settings the Bayesian optimizer is suggesting (vect - len:d) - d: number of parameters
    - PreviousSettings: Series of settings previously measured used to build the gaussian processes (vect - nxd)
      In the dynamic case, the first column of PreviousSettings should be the time (absolute or relative) of the associated query
    - Y: series of output (beam current) associated with the previous settings
    - S: series of control (stability) associated with the prervious settings
    - Opt: optimizer object
    - limits: search boundaries for the parameters (2xd)
  
