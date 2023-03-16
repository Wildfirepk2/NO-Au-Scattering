# adiabatic molecular dynamics (MD) of NO scattering off of an Au(111) surface

############################################################################################################

# if running on isaac or not
isaac=true

# if generating multiple no-au trajectories
multirun=true

### scale down factor on steps for debugging. f=1 -> no scaling.
# au: f=500 -> 1 step
# no/au: f=1e4 -> 1 step
scalefactor=1 #\fix to 1

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

# actual steps for au equilibration. maybe edit later to always be divisable by 10
const steps_eq::Int64 = scalefactor>5000 ? 1 : param.Nsteps_eq[1]/scalefactor

# actual steps for no/au scattering. maybe edit later to always be divisable by 10
const steps_dyn::Int64=param.Nsteps_dyn[1]/scalefactor

# actual steps for logging. small number of actual steps: log every step. otherwise use default step log value
const actsteplog = steps_eq<=100 ? 1 : stepslogging

### description of Au run
aurundesc="Au_slab"

### description of NO/Au run
noaurundesc="NO-Au_sc"

# choosing PESs for NO/Au scattering. all true: diabatic PES
const neutral_PES_active=true
const ionic_PES_active=true
const coupled_PES_active=true

############################################################################################################

# check if Au slab is equilibrated. if not, equilibrate Au slab (MD with velocity verlet)
runAuSlabEquilibration()

############################################################################################################

# NO/Au scattering (MD with velocity verlet)

# store eigenvalues for f func
headers = ["Eg", "λ1", "λ2"]
# fix later for specific types
storeEs=DataFrame([name => [] for name in headers])

# store NO forces from AuN
FNO_AuN=SVector[]

# first Au atom of last layer. last layer (atoms 397-528) is frozen. may put in au var
auatomcutoff=399

# 3 x N matrix. xyz coords of each Au atom stored in columns
rAu=initAuCoords()

# nn: array of arrays. nearest neighbors for each Au atom. ith row corresponds to Au atom i's nearest neighbors (in terms of atom number)
# nn_molly: molly neighbor object. same as nn. for molly compatibility
nn,nn_molly=getnn()
nntuples_au=[(nn_molly.list[i][1]+2,nn_molly.list[i][2]+2,nn_molly.list[i][3]) for i in eachindex(nn_molly.list)]

# array of matrices. initial distances to nearest neighbors for each atom. nn xyz coords stored in columns
r0ij=getdrij(rAu)

# array of arrays. case number for nearest neighbors for each atom. 
nncase=getcase()

# Dict mapping nn case number to force matrix.
Aij=initAij()

# array of arrays of matrices. force matrix for nearest neighbors for each atom. 
Aijarray=initAijarray()

# NO scattering off of eq Au surface. \debug
if multirun
    ts=(300:50:300)u"K"
    xys=(3:2:21)u"Å"
    eis=(25:25:100)u"kJ/mol"

    # # \debug vals
    # ts=(300:50:350)u"K"
    # eis=(20:20:40)u"kJ/mol"
    # xys=(1:2:3)u"Å"

    vary="T_Ei_xy"
    h=["T", "Ei", "xypos", "xNf", "yNf", "zNf", "xOf", "yOf", "zOf", "vxNf", "vyNf", "vzNf", "vxOf", "vyOf", "vzOf", "KEtot", "Etrans", "Erot", "KEvib"]
    trajtrap=DataFrame([name => [] for name in h])
    trajscatter=DataFrame([name => [] for name in h])
    h2=["T", "Ei", "n_scatter", "n_trap", "frac_scatter", "frac_trap"]
    counttraj=DataFrame([name => [] for name in h2])

    date=Dates.format(now(), "yyyy-mm-ddTHHMMSS")
    desc = isaac ? "NO-Au_sc_multi_runs-ISAAC" : "NO-Au_sc_multi_runs"
    outpath=mkpath("results/$desc-$vary--$date")

    for i in eachindex(ts)
        param.T[1]=ts[i]
        T=Int64(ustrip(u"K",param.T[1]))
        global aurundesc="Au slab-T $T"
        runAuSlabEquilibration()
        tpath=mkpath("$outpath/T $T")

        for j in eachindex(eis)
            no.Et_i[1]=eis[j]
            ei=round(ustrip(u"e_MD",no.Et_i[1]);digits=2)
            eipath=mkpath("$tpath/Ei $ei")

            nscatter=0
            ntrap=0
            for k in eachindex(xys)
                xypos=xys[k]
                xy=ustrip(u"Å",xypos)
                xypath=mkpath("$eipath/xy $xy")

                sys=runNOAuTrajectory(xypath,xypos)
                if !checkEconserved(sys)
                    error("Energy not conserved")
                end

                finalNOinfo=finalE_NO(sys)
                allinfo=vcat([T, ei, xy], finalNOinfo)
                if checkscattering(sys)
                    nscatter+=1
                    push!(trajscatter, allinfo)
                else
                    ntrap+=1
                    push!(trajtrap, allinfo)
                end
            end

            totaltraj=ntrap+nscatter
            fracscatter=nscatter/totaltraj
            fractrap=ntrap/totaltraj
            push!(counttraj, [T, ei, nscatter, ntrap, fracscatter, fractrap])
        end
    end

    outputmultirunsummary(vary,counttraj,outpath)
    outputtrajinfo(trajscatter,trajtrap,outpath)
else
    xypos=25u"Å"
    sys=runNOAuTrajectory(xypos)
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
