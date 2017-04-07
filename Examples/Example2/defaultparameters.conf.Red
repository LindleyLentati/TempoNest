
# General parameters:


# useGPUS = 1 for no, 1 for yes
useGPUS = 0;
sample = 1;
sampler = 0;


# Root of the results files,relative to the directory in which TempoNest is run. This will be followed by the pulsar name, and then the individual output file extensions.
root = results/Example-Red-4-

# numTempo2its - sets the number of iterations Tempo2 should do before setting the priors.  Should only be set to 0 if all the priors are set in setTNPriors
numTempo2its = 2;


# ModelChoice

#White noise parameters

incEFAC = 0; # include EFAC: 0 = none, 1 = one for all residuals, 2 = one for each 'whiteflag' (see below)
incEQUAD = 0; # include EQUAD: 0 = no, 1 = yes
incNGJitter = 0;
whiteflag = -group  #whiteflag defines the flag in the tim file that will be used to determine how many EFACS/EQUADs are used when incEFAC or incEQUAD are set to 2
whitemodel = 0; #whitemodel = 0: uses the temponest defintion of EFAC/EQUAD - err^2 = (efac*toaerr)^2 + equad^2, whitemodel=1: uses tempo2 definition of EFAC/EQUAD - err^2 = efac^2*(toaerr^2 + equad^2)
useOriginalErrors = 0;


#Red noise and DM 

incRED = 3; # include Red Noise model: 0 = no, 1 = power law model (vHL2013), 2 = model independent (L2013), 3 = power law model (L2013)
incDM = 0; # include Red Noise model: 0 = no, 1 = power law model (vHL2013), 2 = model independent (L2013), 3 = power law model (L2013)
incBandNoise = 0;
incGroupNoise = 0;

incGlitch = 0;
incGlitchTerms = 0;
GlitchFitSig = 100;

doTimeMargin = 2; # 0=No Analytical Marginalisation over Timing Model. 1=Marginalise over QSD. 2=Marginalise over all Model params excluding jumps.
doJumpMargin = 1; # 0=No Analytical Marginalisation over Jumps. 1=Marginalise over Jumps.



# Priors

RedPriorType = 0; # 0 = Log, 1 = Uniform
DMPriorType = 0;   # 0 = Log, 1 = Uniform
EQUADPriorType = 0;   # 0 = Log, 1 = Uniform
EFACPriorType = 1;   # 0 = Log, 1 = Uniform

# Which priors to use: customPriors=0 uses the Priors from tempo2 fit, along with values set in this function, =1:set priors for specific parameters in setTNPriors
customPriors = 0; 


# FitSig sets the priors for all timing model and jump parameters for both non linear and linear timing models.
# For the non linear fit, Fitsig multiples the error returned by Tempo2, and sets the prior to be the best fit value returned by tempo2 +/- the scaled error.
#  For the linear fit, multiplies the ratio of the rms of the designmatrix vector for each timing model parameter, and the rms of the residuals returned by Tempo2.
FitSig = 10;

# Remaining priors for the stochastic parameters.  
EFACPrior[0] = -1;
EFACPrior[1] = 1.2;


EQUADPrior[0] = -10;
EQUADPrior[1] = -5;

numRedPL = 1;
numDMPL = 1;

varyRedCoeff = 0;
varyDMCoeff = 0;

AlphaPrior[0] = 0;
AlphaPrior[1] = 7;


AmpPrior[0] = -20;
AmpPrior[1] = -8;

numRedCoeff = 120;
numDMCoeff = 120;


CoeffPrior[0] = -10;
CoeffPrior[1] = 0;

DMAlphaPrior[0] = 0;
DMAlphaPrior[1] = 7;


DMAmpPrior[0] = -20;
DMAmpPrior[1] = -8;

FloatingDM = 0;
DMFreqPrior[0] = 0;
DMFreqPrior[1] = 100;

FloatingRed = 0;
RedFreqPrior[0] = 0;
RedFreqPrior[1] = 100;


numBandNoiseCoeff = 60;
BandNoiseAmpPrior[0] = -20;
BandNoiseAmpPrior[1] = -10;
BandNoiseAlphaPrior[0] = 0;
BandNoiseAlphaPrior[1] = 7;

FitForBand[0][0] = 0;
FitForBand[0][1] = 730;
FitForBand[0][2] = 0;
FitForBand[0][3] = 0;

FitForBand[1][0] = 750;
FitForBand[1][1] = 890;
FitForBand[1][2] = 0;
FitForBand[1][3] = 0;





numGroupCoeff = 120;
GroupNoiseFlag = -rnoise
GroupNoiseAmpPrior[0] = -20;
GroupNoiseAmpPrior[1] = -10;
GroupNoiseAlphaPrior[0] = 0;
GroupNoiseAlphaPrior[1] = 7;
FitForGroup[0][0] = 0;
FitForGroup[0][1] = 0;

FitForGroup[1][0] = 1;
FitForGroup[1][1] = 0;


yearlyDM = 0;
DMCoeffPrior[0] = -10;
DMCoeffPrior[1] = 0;

incDMShapeEvent = 0;
numDMShapeCoeff = 0;
DMShapeCoeffPrior[0] = -0.01;
DMShapeCoeffPrior[1] = 0.01;
DMEventStartPrior[0] = 52905;
DMEventStartPrior[1] = 55600;
DMEventLengthPrior[0] = 14;
DMEventLengthPrior[1] = 365.25;


#Set Custom priors to  overwrite the default values for the priors sent to multinest, and the long double priors used by tempo2, you need to be aware of what dimension is what if you use this function.

#THe order of the parameters is always the same:
#Timing Model parameters (linear or non linear)
#Jumps
#EFAC(s)
#EQUAD
#Red Noise Parameters (Amplitude then Alpha for incRed=1, coefficients 1..n for incRed=2)


IS = 1;
ceff = 0;
efr = 0.3;
nlive = 500;
