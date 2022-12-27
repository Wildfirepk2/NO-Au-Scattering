# adiabatic molecular dynamics (MD) of NO scattering off of an Au(111) surface

# checked: 10/27/22
############################################################################################################

# initialize all functions
include("functions/Dependencies.jl")
include("functions/InitMDUnits.jl")
include("functions/InitAllParams.jl")
include("functions/InitAllAuNOVFunc.jl")
include("functions/InitAllAuVFunc.jl")
include("functions/InitCustomMollyFunc.jl")
include("functions/InitSimAnalysisFunc.jl")

############################################################################################################

# parameter variables (type: DataFrame) from input files

# Au parameters: coordinates (x, y, z), number of atoms (N), distance to nearest neighbor (dnn), unit cell length (a), simulation box dimensions (aPBCx, aPBCy, aPBCz), force matrix constants (α, β, γ), mass (m)
au=initAuParams()

# simulation parameters: number of trajectories (Ntraj), temperature (T), time step (dt), number of steps for Au slab equilibration (Nsteps_eq), number of steps for MD (Nsteps_dyn)
param=initSimParams()

# ground state PES parameters
PES_GS=initPESParamGS()

# excited state PES parameters
PES_ionic=initPESParamIonic()

# coupled state PES parameters
PES_coup=initPESParamCoup()

# NO parameters: velocity of vibration (v) as a function of bond length (r), initial vibrational state (Nv), number of v,r pairs (Nvib), initial translational energy (Et_i), initial incidence angle (θi), mass of nitrogen and oxygen (mN, mO)
no=initNOParams()

############################################################################################################

# MD related variables

# first Au atom of last layer. last layer (atoms 397-528) is frozen. may put in au var
const auatomcutoff=397

# log parameters after n steps. may put in param var
const stepslogging=10

### scale down factor on steps for debugging. f=1 no scaling. f=5000: 1 step
scalefactor=1000

# actual steps for au equilibration. maybe edit later to always be divisable by 10
const steps_eq::Int64=param.Nsteps_eq[1]/scalefactor

# actual steps for logging. small number of actual steps: log every step. otherwise use default step log value
const actsteplog = steps_eq<=100 ? 1 : stepslogging

### description of run
aurundesc="Au slab"

############################################################################################################

# other variables

# copied from fortran
Av_Et = 0
Av_theta = 0
Av_weighted_theta = 0
Norm_weighted_theta = 0
Av_Evib = 0
Av_Er = 0
Av_EAu = 0
NTRAP = 0

############################################################################################################

# check if Au slab is equilibrated. if not, equilibrate Au slab

# make helper function for outputting index of au folder. adapt for outputting most recent run?
resultsdir=readdir("results";join=true)
i_au=findfirst(contains.(resultsdir,"Au slab"))

# leave in main
if i_au isa Nothing
    runAuSlabEquilibration()
end

# make file for holding no/au related functions
# make function for no/au system init. function composed initNO and initeqAu
resultsdir=readdir("results";join=true)
i_au=findfirst(contains.(resultsdir,"Au slab"))
audir=resultsdir[i_au]
coordsfile="$audir/syscoords.xlsx"
xfcoord=XLSX.readxlsx(coordsfile)
sheets=XLSX.sheetnames(xfcoord)
sheetlastcoord=sheets[end]
dfcoord=DataFrame(XLSX.readtable(coordsfile,sheetlastcoord))
aueqcoords=[SA[dfcoord[i,1],dfcoord[i,2],dfcoord[i,3]] for i in 1:nrow(dfcoord)]

# make simulator/system variables for no/au in global. later on put in separate file like au slab equilibration