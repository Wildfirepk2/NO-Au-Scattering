# adiabatic molecular dynamics (MD) of NO scattering off of an Au(111) surface

############################################################################################################

# global settings

# if running on isaac or not
const isaac=false

# if generating multiple no-au trajectories
const multirun=false

# if running random trajectories (changes initial NO pos)
const randomtraj=true

# energy in kJ/mol of NO runs. \fix
const erun=25

### description of NO/Au run. \fix
const noaurundesc="NO-Au_sc_E$erun"

### description of O/Au run. \fix
const oaurundesc="O-Au_sc_E$erun"

# if debugging, cut steps to 1, multirun: 2 params T,E,xy
const debug=false

# if wanting simple results. no animation, no fv excels, etc
const simplerun=false

# if doing a O/Au run
const runningoau=false

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

### scale down factor on steps for debugging. f=1 -> no scaling.
# au: f=500 -> 1 step
# no/au: f=1e4 -> 1 step
scalefactor = debug ? 1e4 : 1

# actual steps for au equilibration. maybe edit later to always be divisable by 10
const steps_eq::Int64 = scalefactor>5000 ? 1 : param.Nsteps_eq[1]/scalefactor

# actual steps for no/au scattering. maybe edit later to always be divisable by 10
const steps_dyn::Int64=param.Nsteps_dyn[1]/scalefactor

# actual steps for o/au scattering. maybe edit later to always be divisable by 10
const steps_dyn_OAu::Int64 = scalefactor>5000 ? 1 : param.Nsteps_dyn_OAu[1]/scalefactor

# actual steps for logging. small number of actual steps: log every step. otherwise use default step log value
const actsteplog = steps_eq<=100 ? 1 : stepslogging

### description of Au run \fix const
aurundesc="Au_slab"

# choosing PESs for NO/Au scattering. all true: diabatic PES
const neutral_PES_active=true
const ionic_PES_active=true
const coupled_PES_active=true

############################################################################################################

# main variables

# store eigenvalues for f func
headers = ["Eg", "λ1", "λ2"]
# fix later for specific types
storeEs=DataFrame([name => [] for name in headers])

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

if runningoau
    # O/Au scattering (MD with velocity verlet)
    if multirun
        ts, eis, xs, ys = initTExy()
        vary="T_Ei_xy"
        trajtrap, trajscatter, counttraj = inittrajcontainers()

        desc = isaac ? "$(oaurundesc)_multi_runs-ISAAC" : "$(oaurundesc)_multi_runs"
        outpath=makeresultsfolder(desc,vary)

        runMultiOAuTrajectory(ts, eis, xs, ys)
    else
        # check if Au slab is equilibrated. if not, equilibrate Au slab (MD with velocity verlet)
        runAuSlabEquilibration()

        if randomtraj
            xNOi=au.aPBCx[1]*rand()
            yNOi=au.aPBCy[1]*rand()
        else
            xNOi=25u"Å"
            yNOi=xNOi
        end
        sys=runOAuTrajectory(xNOi,yNOi)
    end
else
    # NO/Au scattering (MD with velocity verlet)
    if multirun
        ts, eis, xs, ys = initTExy()
        vary="T_Ei_xy"
        trajtrap, trajscatter, counttraj = inittrajcontainers()

        desc = isaac ? "$(noaurundesc)_multi_runs-ISAAC" : "$(noaurundesc)_multi_runs"
        outpath=makeresultsfolder(desc,vary)

        runMultiNOAuTrajectory(ts, eis, xs, ys)
    else
        # check if Au slab is equilibrated. if not, equilibrate Au slab (MD with velocity verlet)
        runAuSlabEquilibration()

        if randomtraj
            xNOi=au.aPBCx[1]*rand()
            yNOi=au.aPBCy[1]*rand()
        else
            xNOi=25u"Å"
            yNOi=xNOi
        end
        sys=runNOAuTrajectory(xNOi,yNOi)
    end
end
