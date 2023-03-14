# adiabatic molecular dynamics (MD) of NO scattering off of an Au(111) surface

############################################################################################################

# if running on isaac or not
isaac=true

############################################################################################################

# initialize all functions
include("functions/Dependencies.jl")
include("functions/MDUnits.jl")
include("functions/InitSimParams.jl")
include("functions/NOAuFVFunc.jl")
include("functions/InitAuSys.jl")
include("functions/InitNOAuSys.jl")
include("functions/SupportMollyFunc.jl")
include("functions/CustomMollyFunc.jl")
include("functions/SimAnalysisFunc.jl")

############################################################################################################

# if generating multiple no-au trajectories
multirun=true
xypos=25u"Å"
globalpath="results"
vary="T Ei xy"
ntrap=0
nscatter=0
trajtrap=Tuple[]
trajscatter=Tuple[]

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
# au: f=5000 -> 1 step
# no/au: f=1e5 -> 1 step
scalefactor=20

# actual steps for au equilibration. maybe edit later to always be divisable by 10
const steps_eq::Int64 = scalefactor>5000 ? 1 : param.Nsteps_eq[1]/scalefactor*2

# actual steps for no/au scattering. maybe edit later to always be divisable by 10
const steps_dyn::Int64=param.Nsteps_dyn[1]/scalefactor

# actual steps for logging. small number of actual steps: log every step. otherwise use default step log value
const actsteplog = steps_eq<=100 ? 1 : stepslogging

### description of Au run
aurundesc="Au slab-recent-last layer froz"

### description of NO/Au run
noaurundesc="NO-Au sc-d-DIA-Au move-100ei-xy25-react F Au-VAu"

multirundesc="NO-Au sc multi runs-$vary"

# choosing PESs for NO/Au scattering. all true: diabatic PES
const neutral_PES_active=true
const ionic_PES_active=true
const coupled_PES_active=true

############################################################################################################

# running MD with velocity verlet

# check if Au slab is equilibrated. if not, equilibrate Au slab
runAuSlabEquilibration()

# store eigenvalues for f func
headers = ["Eg", "λ1", "λ2"]
# fix later for specific types
storeEs=DataFrame([name => [] for name in headers])

# store NO forces from AuN
FNO_AuN=SVector[]

# NO scattering off of eq Au surface. \debug
if multirun
    ts=(300:50:300)u"K"
    xys=(3:2:21)u"Å"
    eis=(25:25:100)u"kJ/mol"

    # # debug vals
    # ts=(300:50:350)u"K"
    # xys=(1:2:1)u"Å"
    # eis=(20:20:20)u"kJ/mol"

    date=Dates.format(now(), "yyyy-mm-ddTHHMMSS")
    outpath=mkpath("$globalpath/$multirundesc--$date")

    for i in eachindex(ts)
        param.T[1]=ts[i]
        T=Int64(ustrip(u"K",param.T[1]))
        global aurundesc="Au slab-T $T"
        runAuSlabEquilibration()
        tpath=mkpath("$outpath/T $T")

        for j in eachindex(xys)
            global xypos=xys[j]
            xy=ustrip(u"Å",xypos)
            xypath=mkpath("$tpath/xy $xy")

            for k in eachindex(eis)
                no.Et_i[1]=eis[k]
                ei=Int64(ustrip(u"kJ/mol",no.Et_i[1]))
                eipath=mkpath("$xypath/Ei $ei")
                global globalpath="$eipath"
                global noaurundesc="NO-Au sc-T $T-$(ei)ei-xy$xy"

                runNOAuTrajectory()

                trajparam=(ts[i], xys[j], eis[k])
                if checkscattering(sys_NOAu)
                    global nscatter+=1
                    push!(trajscatter, trajparam)
                else
                    global ntrap+=1
                    push!(trajtrap, trajparam)
                end
            end
        end
        global globalpath="results"
    end

    outputmultirunsummary(outpath)
else
    runNOAuTrajectory()
end

############################################################################################################

# make file for holding no/au related functions

# make simulator/system variables for no/au in global. later on put in separate file like au slab equilibration

# other variables
# copied from fortran
# Av_Et = 0
# Av_theta = 0
# Av_weighted_theta = 0
# Norm_weighted_theta = 0
# Av_Evib = 0
# Av_Er = 0
# Av_EAu = 0
# NTRAP = 0
