# adiabatic molecular dynamics of NO scattering off of an Au(111) surface

# initialize all functions
include("functions/InitMDUnits.jl")
include("functions/InitAllParams.jl")
include("functions/InitAllAuNOVFunc.jl")
include("functions/InitAllAuVFunc.jl")

#initialize variables from input files
au=initAuParams()
param=initSimParams()
PES_GS=initPESParamGS()
PES_ionic=initPESParamIonic()
PES_coup=initPESParamCoup()
no=initNOParams()

#other variables
mAu=196.966569
mN=14.00674
mO=15.9994

Av_Et = 0
Av_theta = 0
Av_weighted_theta = 0
Norm_weighted_theta = 0
Av_Evib = 0
Av_Er = 0
Av_EAu = 0
NTRAP = 0

rAu=initAuCoords() # 3 x N matrix. xyz coords of each Au atom stored in columns
nn=getnn() # nearest neighbor array for each Au atom. ith row corresponds to Au atom i's nearest neighbors (in terms of atom number)
r0ij=getdrij(rAu) # array of vectors to each Au atom's nearest neighbors. ith row cooresponds to all vectors (in [x,y,z] format) between Au atom i and its nearest neighbors
nncase=getcase() # case # array for each Au atom. ith row cooresponds to the case #'s between Au atom i and its nearest neighbors 
Aij=initAij() # force matrix as a function of case #
Aijarray=initAijarray()

run(`pwd`)