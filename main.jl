# adiabatic molecular dynamics (MD) of NO scattering off of an Au(111) surface

############################################################################################################

# global settings

### if doing a NO/Au run
const runningnoau=false

### if doing a O/Au run
const runningoau=false

# if debugging, multirun: 1 set varied params (T,E,orient,xy). if isaac, false
debug=false

# if wanting 1 step only on trajs. if isaac, false
shortrun=false

# if wanting simple results. no animation, no fv excels, etc. if isaac, true
simplerun=false

############################################################################################################

# if running on isaac or not
const isaac = Sys.total_memory()/10^9 > 100 # GB

# ensure right isaac settings
if isaac;debug=false;shortrun=false;simplerun=true;end

# description of NO/Au run. \fix
const noaurundesc = isaac ? "NO-Au_sc-ISAAC" : "NO-Au_sc"

# description of O/Au run. \fix
const oaurundesc = isaac ? "O-Au_sc-ISAAC" : "O-Au_sc"

############################################################################################################

# initialize all functions
include("functions/Dependencies.jl")
include("functions/MDUnits.jl")
include("functions/InitSimParams.jl")
include("functions/NOAuVFunc.jl")
include("functions/NOAuFFunc.jl")
include("functions/InitAuSys.jl")
include("functions/InitNOAuSys.jl")
include("functions/InitOAuSys.jl")
include("functions/SupportMollyFunc.jl")
include("functions/CustomAuMollyFunc.jl")
include("functions/CustomNOAuMollyFunc.jl")
include("functions/CustomOAuMollyFunc.jl")
include("functions/SimAnalysisFunc.jl")
include("functions/AuSimAnalysisFunc.jl")
include("functions/NOAuSimAnalysisFunc.jl")
include("functions/OAuSimAnalysisFunc.jl")
include("functions/MultirunSupport.jl")

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

# Molly object. simulation box dimensions. periodic in x,y directions. needed in other areas outside Molly
# virtual box dims. needed for certain molly functions
simboxdims=CubicBoundary(au.aPBCx[1], au.aPBCy[1], Inf*u"Å")
virtboxdims=CubicBoundary(au.aPBCx[1], au.aPBCy[1], au.aPBCz[1])

# log parameters after n steps. may put in param var
const stepslogging=10

# actual steps for au equilibration
const steps_eq::Int64 = shortrun ? 1 : param.Nsteps_eq[1]

# actual steps for no/au scattering
const steps_dyn::Int64= shortrun ? 1 : param.Nsteps_dyn[1]

# actual steps for o/au scattering
const steps_dyn_OAu::Int64 = shortrun ? 1 : param.Nsteps_dyn_OAu[1]

# actual steps for logging. small number of actual steps: log every step. otherwise use default step log value
const actsteplog = steps_eq<=100 ? 1 : stepslogging

### description of Au run \fix const?
aurundesc="Au_slab"

# choosing PESs for NO/Au scattering. all true: diabatic PES
const neutral_PES_active=true
const ionic_PES_active=true
const coupled_PES_active=true

# temp to test NO orientations
const Torient=300u"K"

# actual no of trajectories
acttraj=getacttraj()

############################################################################################################

# main variables

# store eigenvalues for f func
storeEs=DataFrame()

# store NO forces from AuN
FNO_AuN=SVector[]

# first Au atom of last layer. last layer (atoms 399-530) is frozen. redefined as 397 for Au eq (since theres no N,O) and turned back to 399 after
auatomcutoff=399

# 3 x N matrix. xyz coords of each Au atom stored in columns
rAu=initAuCoords()

# nn: array of arrays. nearest neighbors for each Au atom. ith row corresponds to Au atom i's nearest neighbors (in terms of atom number)
# nn_molly: molly neighbor object. same as nn. for molly compatibility
nn,nn_molly=getnn()
nntuples_au=[(nn_molly.list[i][1]+2,nn_molly.list[i][2]+2,nn_molly.list[i][3]) for i in eachindex(nn_molly.list)]
nntuples_oau=[(nn_molly.list[i][1]+1,nn_molly.list[i][2]+1,nn_molly.list[i][3]) for i in eachindex(nn_molly.list)]

# array of matrices. initial distances to nearest neighbors for each atom. nn xyz coords stored in columns
r0ij=getdrij(rAu)

# array of arrays. case number for nearest neighbors for each atom. 
nncase=getcase()

# Dict mapping nn case number to force matrix.
Aij=initAij()

# array of arrays of matrices. force matrix for nearest neighbors for each atom. 
Aijarray=initAijarray()

############################################################################################################
#\fix
# collectallresults=makeresultsfolder("results/run") 
# path=collectallresults

# NO/Au scattering (MD with velocity verlet)
if runningnoau
    runMultiNOAuTrajectory()
end

# NO/Au scattering. fixed orientation runs
if !all(ismissing,no.θorient)
    runMultiNOAuTrajectory(fixorient=true,)
end

# # NO/Au scattering. fixed vib phase runs
# if !all(ismissing,no.vib_phase)
#     runMultiNOAuTrajectory(fixorient=true,fixvib=true,)
# end

# O/Au scattering
if runningoau
    runMultiOAuTrajectory()
end

# #\debug
# sys=runNOAuTrajectory()
# sys=runOAuTrajectory()
# runAuSlabEquilibration(force=true)
