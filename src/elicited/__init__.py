##########################################################################
#
# Functions in this library:
#
# elicitLogNormal
# elicitPERT
# elicitPareto
# elicitZipf
#
##########################################################################



##########################################################################
#
# PYTHON 3 FUNCTION
#
# logN_mu,logN_sig = elicitLogNormal(modeX,quantX,quantP=0.95)
#
# This function solves for the parameters describing the mean and standard
# deviation of a log-normal distribution of X from eliciting two values of 
# the distribution: the most common value of X (modeX) and the 
# quantile value (quantX) for the (quantP) quantile. We can understand
# quantX as the value for which P(X<=quantX) = quantP, or the value whose
# exceedance probability is 1 - quantP.
#
# The values of logN_mu and logN_sig correspond to the distribution of
# Y = exp(X), or the values of the distribution projected into log-space.
# Once these values have been computed, a log-normal distribution can be
# described using scipy.stats.lognorm(s=logN_sig,scale=exp(logN_mu)).
#
# INPUTS:
#    modeX ..................................... most common value of X
#    quantX .................................... quantile value of X
#    quantP .................................... quantile probability of quantX (i.e. exceedance probability of quantX is 1 - quantP): Default value of quantP=0.95
#
# OUTPUTS:
#    logN_mu ................................... mean of Y = exp(X)
#    logN_sig .................................. standard deviation of Y = exp(X)
#
##########################################################################
#
# DESCRIPTION OF ELICITATION STRATEGY
#
# The log-normal distribution has two parameters, logN_mu and logN_sig,
# which correspond directly to the mean and standard deviation of 
# Y = exp(X). These are not eaily elicited values from an expert, and we
# need to rely on a few properties of the log-normal distribution to allow
# for a translation between the questions we *can* ask an expert and the
# parameter values we need.
#
# The mode of X (modeX), which is the peak of the log-normal PDF and
# represents the most common value of X. It is related to logN_mu and
# logN_sigma by:
#
# (1) modeX = exp(logN_mu - logN_sig**2)
#
# We will define some quantile value of X (quantX; i.e., values of 
# X<=quantX have a probability of quantP, and quantX has an exceedance
# probability of 1 - quantP). This quantile of the log-normal (logN_q) is
# related to logN_mu and logN_sigma by:
#
# (2) quantX = logN_q(quantP) = exp(logN_mu + logN_sig * N_q(quantP))
#
# Where N_q(quantP) is the quantP quantile value of a standard
# normal distribution, which is easily computed.  
#
# Equations (1) and (2) have two unknowns (logN_mu and logN_sig), making
# this system solvable through combining equations. We can combine by
# subtracting Equation (2) from Equation (1), which produces a single
# equation containing only logN_sig as an unknown, and yield the
# quadratic expression:
#
# (3) logN_sig**2 + N_q(quantP) * logN_sig + log(modeX) - log(quantX) = 0
#
# Thus, Equation (3) can be solved via the quadratic equation, yielding 
# two roots as solutions for logN_sigma as:
#
# (4) logN_sigma = (-N_q(quantP) +/- sqrt(N_q(quantP)**2 - 4*(log(modeX)-log(quantX))))/2
#
# This yields two potential values for logN_sigma that are positive, 
# negative, or complex. We can disregard negative and complex solutions as 
# nonphysical, and retain any positive+real solutions (ideally there will be
# one solution of this kind expected, although I have no proof of that).
#
# With a valid logN_sig, we can then combine:
#
# N_q(quantP) * (Equation-1) + logN_sig * (Equation-2) and solve for
# logN_mu:
#
# (5) logN_mu = (N_q(quantP) * log(modeX) + logN_sig * log(quantX))/(N_q(quantP)+logN_sig)
#
# Yielding the required parameters for the system. The resulting log-normal
# distribution will peak at X=modeX and its quantP quantile will be quantX.
#
##########################################################################
def elicitLogNormal(modeX,quantX,quantP=0.95):
    ######################################################################
    #
    # Load required modules
    #
    import numpy as np #.................................................. array module
    from scipy.stats import norm #........................................ normal distribution module
    #
    ######################################################################
    #
    # Initialize output as None
    #
    logN_mu = None #...................................................... mean of log-normal (initialized to None)
    logN_sig = None #..................................................... st. dev. of log-normal (initialized to None)
    #
    ######################################################################
    #
    # Solve for logN_sig as roots of polynomial (Equation-3, above)
    #
    # Define N_q(quantP) and generate polynomial
    N_q = norm.ppf(quantP) #.............................................. quantile of standard normal distribution at P=quantP
    p=np.polynomial.polynomial.Polynomial( #.............................. polnomial expression
                                           [
                                             np.log(modeX)-np.log(quantX) ,
                                             N_q                        ,
                                             1.
                                           ]
                                         )
    # The roots may be unsolvable, so we are using a try/except block and
    # defaulting to generating an empty list for roots if the polynomial
    # roots are not solved.
    try:
        p_roots = p.roots() #............................................. roots of p
        # Remove any complex roots
        p_roots = np.delete(p_roots,np.iscomplex(p_roots))
        # Remove any negative or zero-roots
        p_roots = np.delete(p_roots,np.where(p_roots<=0.))
    except:
        p_roots = []
    #
    ######################################################################
    #
    # If any positive+real roots were discovered, solve for logN_mu for
    # each of them (will yield multiple solutions if len(p_roots)>1)
    #
    if (len(p_roots)==0):
        print('no real roots >0 found - UNSOLVABLE DISTRIBUTION')
        # Return
        return logN_mu, logN_sig
    else:
        # Compute mean values for each valid root
        logN_mu = np.nan * np.ones((len(p_roots))) #...................... vector of logN_mu values (initialized to NaN)
        logN_sig = np.nan * np.ones((len(p_roots))) #..................... vector of logN_sig values (initialized to NaN)
        for i in range(len(p_roots)):
            logN_sig[i] = p_roots[i]
            # Solve for mean
            m = (N_q*np.log(modeX)+p_roots[i]*np.log(quantX))/(N_q+p_roots[i])
            logN_mu[i] = m
    # Return
    return logN_mu, logN_sig










##########################################################################
#
# PYTHON 3 FUNCTION
#
# PERT_alpha,PERT_beta = elicitPERT(minX,modeX,maxX)
#
# This function solves for the shape parameters (alpha, beta) describing 
# the beta distribution of X from eliciting three values of the 
# distribution: the minimum value of X (minX), the most common value of X 
# (modeX) and the maximum value of X (maxX). The values of alpha and beta
# are further constrained by the PERT distribution's requirement that the
# mean of the distribution is defined by:
#
# meanX = (minX+4*modeX+maxX)/6
#
# The PERT distribution can then be described using 
# scipy.stats.beta(a=PERT_alpha,b=PERT_beta,loc=minX,scale=maxX-minX).
#
# INPUTS:
#    minX ...................................... smallest possible value of X
#    modeX ..................................... most common value of X
#    maxX ...................................... largest possible value of X
#
# OUTPUTS:
#    PERT_alpha ................................ alpha shape-parameter for (PERT-)beta distribution of X
#    PERT_beta ................................. beta shape-parameter for (PERT-)beta distribution of X
#
##########################################################################
#
# DESCRIPTION OF ELICITATION STRATEGY
#
# The elicitation strategy of a PERT distribution is straightforward,
# since the distribution is directly controlled by the elicitable terms
# of (minX,modeX,maxX). The values of PERT_alpha and PERT_beta are then
# computed by utilizing the following equations:
#
# (1) meanX = (minX+4*modeX+maxX)/6
# (2) PERT_alpha = (4*modeX+maxX-5*minX)/(maxX-minX)
# (3) PERT_beta = (5*maxX-minX-4*modeX)/(maxX-minX)
#
# For both PERT_alpha and PERT_beta, we will define using minX, maxX, and
# meanX, via use of Equation-1. For PERT_alpha, via Equation-1 let:
#
# (4) 4*modeX = 6*meanX-minX-maxX
#
# Applying to Equation-2:
#
# (5) PERT_alpha = ((6*meanX-minX-maxX)+maxX-5*minX)/(maxX-minX)
#
# Which simplifies to:
#
# (5b) PERT_alpha = 6*(meanX-minX)/(maxX-minX)
#
# Applying to Equation-3:
#
# (6) PERT_beta = (5*maxX-minX-(6*meanX-minX-maxX))/(maxX-minX)
#
# Which simplifies to
#
# (6b) PERT_beta = 6*(maxX-meanX)/(maxX-minX)
#
# Yielding the required parameters for the system.
#
def elicitPERT(minX,modeX,maxX):
    ######################################################################
    #
    # Load required modules
    #
    import numpy as np #.................................................. array module
    #
    ######################################################################
    #
    # Initialize output as None
    #
    PERT_alpha = None #................................................... alpha shape-parameter for (PERT-)beta distribution (initialized to None)
    PERT_beta = None #.................................................... alpha shape-parameter for (PERT-)beta distribution (initialized to None)
    #
    ######################################################################
    #
    # Compute meanX via Equation-1:
    #
    meanX = (minX + 4*modeX + maxX)/6.
    #
    #
    # Compute PERT_alpha and PERT_beta in try/except blocks to handle
    # errors
    #
    # PERT_alpha:
    try:
        PERT_alpha = 6*((meanX - minX)/(maxX - minX))
    except:
        print('Error computing PERT_alpha')
        PERT_alpha = None
    # PERT_beta:
    try:
        PERT_beta = 6*((maxX - meanX)/(maxX - minX))
    except:
        print('Error computing PERT_beta')
        PERT_beta = None
    #
    # Return
    #
    return PERT_alpha, PERT_beta










##########################################################################
#
# PYTHON 3 FUNCTION
#
# Pareto_b = elicitPareto(minX,quantX,quantP=0.95)
#
# This function solves for the shape parameter (b) describing 
# the Pareto distribution of X from eliciting two values of the 
# distribution: the minimum value of X (minX) and the quantile value of X 
# (quantX) for the (quantP) quantile, or the value for which X>quantX
# would be an event with an exceedance probability of 1 - quantP. 
#
# The Pareto distribution can then be described using 
# scipy.stats.pareto(b=Pareto_b,loc=minX-1.,scale=1.).
#
# INPUTS:
#    minX ...................................... smallest possible value of X
#    quantX .................................... quantile value of X at P=quantP
#    quantP .................................... quantile probability of quantX: Default quantP=0.95
#
# OUTPUTS:
#    Pareto_b .................................. b shape-parameter for Pareto distribution of X
#
##########################################################################
#
# DESCRIPTION OF ELICITATION STRATEGY
#
# The elicitation strategy of a Pareto distribution makes use of the
# definition of the quantile of the distribution for some probability P,
# assuming that the scale parameter is 1. (standard) and the distribution
# minimum value is a combination of some standard minimum of X=1 and a
# shift of the distribution by a value of loc (e.g. minX=loc+1):
#
# q(P) = (1.-P)**(-1./b) + loc
#
# We are eliciting a value of quantX, representing the quantile value for
# the quantP quantile, or the value of X for which a higher value of X 
# would be an event with an exceedance probability of 1 - quantP:
#
# q(quantP) = (1.-quantP)**(-1./b) + loc
#
# We can solve for b, the Pareto shape parameter, as:
#
# b = -1./log_[1.-quantP](quantX-loc)
#
# Where log_[1.-quantP]() is the [1.-quantP]-base logarithm. By the logarithm 
# base-change rule, this becomes:
#
# b = -1./(log(quantX-loc)/log(1.-quantP))
#
# Where log() is the natural logarithm, yielding the required parameter 
# for the system.
#
def elicitPareto(minX,quantX,quantP=0.95):
    ######################################################################
    #
    # Load required modules
    #
    import numpy as np #.................................................. array module
    #
    ######################################################################
    #
    # Initialize output as None
    #
    Pareto_b = None #..................................................... b shape-parameter for Pareto distribution (initialized to None)
    #
    ######################################################################
    #
    # Define necessary constants
    #
    loc_val=minX - 1. #................................................... loc parameter for (sciPy) Pareto distribution
    a=(quantX-loc_val) #.................................................. logarithm quantity (see Elicitation Strategy, above)
    exceedP=1. - quantP #................................................. exceedance probability of quantX
    #
    # Define Pareto_b, use try/except block to handle errors
    #
    try:
        Pareto_b=-1./(np.log(a)/np.log(exceedP))
    except:
        print('Error defining Pareto_b, exiting')
    #
    # Return
    #
    return Pareto_b
   









##########################################################################
#
# PYTHON 3 FUNCTION
#
# zipf_s = elicitZipf(minX,quantX,quantP=0.95,report=False)
#
# This function solves for the single shape parameter (s) for a Zipf (or
# zeta) distribution of X from eliciting a minimum value of X, a quantile
# value of X, and the quantile probability of quantX. In this way, quantX
# can be thought of as a very large value of X for which P(x>maxX) is a very
# small probability 1.-quantP.
#
# A Zipf distribution can then be defined via:
# scipy.stats.zipf(zipf_s,quantX-1).
#
# INPUTS:
#    minX ...................................... smallest possible value of X
#    quantX .................................... quantile value of X
#    quantP ...................................... quantile probability for X=quantX: Default quantX=0.95
#    report .................................... boolean for reporting error of solution: Default report=False
#
# OUTPUTS:
#    zipf_s .................................... shape parameter of zipf distribution
#
##########################################################################
#
# DESCRIPTION OF ELICITATION STRATEGY
#
# The cumulative distribution function that defines quantiles for the Zipf
# distribution relies on the Riemann zeta function, an infinite series. I
# am not qualified to attempt an algebraic solution to a problem involving
# such a series, so unlike other distributions the solution here is solved
# empirically rather than analytically. This technique has a few
# characteristics that differ slightly from the analytic solutions for
# other distributions:
#
# 1. The function can take much longer to solve - Typically only a few
#    seconds in testing, but the amount of time it takes to solve the
#    system is related to the amount of tolerable error in the solution
#    (see below)
#
# 2. The function solves for the shape-parameter (s) with some amount of
#    error. The error comes from selecting a value of s which minimizes
#    error, but only searching for a value of s to within some chosen
#    number of significant digits.
#
# Error can be reported (as a percentage of residual P(x<maxX)) through
# use of report=True. Error can be reduced through a smaller value of the
# internal variable s_step, which defines the resolution of the search for
# a best-fit value of s. A smaller value of s_step requires a search
# across a larger number of possible values of s, which reduces the time 
# efficiency of the function.
#
# The elicitation strategy involves computing many Zipf distributions
# across a range of values of s, and then selecting the value that best 
# approximates P(x>quantX)=1.-quantP. This is accomplished through two 
# separate searches:
#
# The first search is a coarse-resolution search to identify a range of 
# values for s down to a value of +/- 0.1. This does not solve the value
# of s, but provides a range in which to perform a higher-resolution
# search. The strategy involves starting with an extremely small value
# for s (s must be >1, starting value is 1 plus some very small value),
# and then looping in increments of 0.1 until a value for s that is too
# large is found. A value too large will have P(x>quantX) < 1.-quantP.
# The first value meeting this criteria is defined as the ceiling for s,
# and the prior value (smaller by 0.1) is defined as the floor.
#
# The second search is performed only within the 0.1 range defined by
# the first search, at a resolution defined by the value of s_step.
# All values between the floor and ceiling at a resolution of s_step
# are used to define Zipf distributions, and P(x>quantX) is computed for
# each of them. The value of s for which abs(P(x>quantX)-[1.-quantP]) is
# minimized is the solved value for s.
#
# Error is defined as 100.*abs(P(x>quantX)-[1.-quantP])/[1.-quantP], or 
# the percent difference between the solved probability of exceedance and
# the elicited value for the probability of exceedance. This can be 
# reported to the user through report=True.
#
##########################################################################
def elicitZipf(minX,quantX,quantP=0.95,report=False):
    ######################################################################
    #
    # Load required modules
    #
    import numpy as np #.................................................. array module
    from scipy.stats import zipf #........................................ Zipf distribution module
    #
    ######################################################################
    #
    # Define k, the localization parameter for the Zipf distribution, as
    # minX - 1. The default is k=0, indicating that the minimum value of
    # the distribution is 1. For a distribution whose minimum (and mode)
    # value is zero, minX=0 will properly adjust k.
    #
    k=minX-1 #............................................................ localization parameter of Zipf distribution 
    #
    # Define exceedP, the exceedance probability for X=quantX
    exceedP=1.-quantP #................................................... exceedance probability P(X>quantX)
    #
    ######################################################################
    #
    # First search: Find search-volume for fitting parameter s
    #
    # Starting with a value of s > 1.0, compute P(x>quantX) in 0.1 
    # increments until a value P<exceedP is discovered. The value of s 
    # where P(x<quantX)<exceedP is found represents a ceiling on s, with
    # the prior searched value representing a floor on s.
    #
    s_floor = 1.0 + 1.0e-05 #............................................. initial floor for search-volume of s (tiny bit larger than 1)
    # Define P(x>quantX) for s = s_floor
    s = s_floor #......................................................... test-value for s
    d = zipf(s,k) #....................................................... test Zipf distribution
    p = 1. - d.cdf(quantX) #.............................................. test P(x>quantX)
    # While p >= exceedP, increment s by 0.1 and update s_floor to track
    # a 0.1-range search-volume for s until a ceiling is found with
    # P(x>quantX) < exceedP
    while p >= exceedP:
        s_floor = s
        s = s + 0.1
        d = zipf(s,k)
        p = 1. - d.cdf(quantX)
    s_ceiling = s #....................................................... ceiling of search-volume for s
    #
    ######################################################################
    #
    # Second search: Perform a detailed search for s between s_floor and 
    #                s_ceiling for best-fit value of s
    #
    # Compute P(x>quantX) for s in s_step increments between s_floor and 
    # s_ceiling and assign s to the value for which 
    # abs(P(x>quantX)-exceedP) is minimized.
    #
    s_step = 0.0001 #..................................................... resolution of search
    s_range = np.arange(s_floor,s_ceiling+s_step/10.,s_step) #............ search range
    p = np.nan * np.ones(np.shape(s_range)) #............................. P(x>maxX) for all search-values (initialized to NaN)
    # Loop through all values of s in s_range, compute P(x>quantX)
    for i in range(np.size(s_range)):
        si = s_range[i]
        d = zipf(si,k)
        p[i] = 1. - d.cdf(quantX)
    # Find value of s in search for minimum error in P(x>quantX)
    idx = np.argmin(np.abs(p-exceedP)) #.................................. index of best-fit s
    s = s_range[idx] #.................................................... value of best-fit s
    #
    ######################################################################
    #
    # For report=True, report value of s and error in best-fit
    #
    if report:
        print('s={:12.6f}'.format(s),'error={:10.4f}%'.format(100.*np.abs(p[idx]-exceedP)/exceedP))
    return s

