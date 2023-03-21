# adiabatic molecular dynamics (MD) of NO scattering off of an Au(111) surface

############################################################################################################

# if running on isaac or not
const isaac=true

# if generating multiple no-au trajectories
const multirun=true

# if running random trajectories (changes initial NO pos)
const randomtraj=true

### description of NO/Au run
const noaurundesc="NO-Au_sc_E25"

# if debugging
const debug=false

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

### scale down factor on steps for debugging. f=1 -> no scaling.
# au: f=500 -> 1 step
# no/au: f=1e4 -> 1 step
scalefactor = debug ? 1e4 : 1

# actual steps for au equilibration. maybe edit later to always be divisable by 10
const steps_eq::Int64 = scalefactor>5000 ? 1 : param.Nsteps_eq[1]/scalefactor

# actual steps for no/au scattering. maybe edit later to always be divisable by 10
const steps_dyn::Int64=param.Nsteps_dyn[1]/scalefactor

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

# array of matrices. initial distances to nearest neighbors for each atom. nn xyz coords stored in columns
r0ij=getdrij(rAu)

# array of arrays. case number for nearest neighbors for each atom. 
nncase=getcase()

# Dict mapping nn case number to force matrix.
Aij=initAij()

# array of arrays of matrices. force matrix for nearest neighbors for each atom. 
Aijarray=initAijarray()

############################################################################################################

# NO/Au scattering (MD with velocity verlet)

if multirun
    if debug
        ts=(300:50:350)u"K"
        eis=(20:20:40)u"kJ/mol"
        if randomtraj
            ntraj=2
            xs=[au.aPBCx[1]*rand() for _ in Base.OneTo(ntraj)]
            ys=[au.aPBCy[1]*rand() for _ in Base.OneTo(ntraj)]
        else
            xs=(1:2:3)u"Å"
            ys=xs
        end
    else
        ts=(300:50:300)u"K"
        eis=(25:25:25)u"kJ/mol"
        if randomtraj
            ntraj=200
            xs=[au.aPBCx[1]*rand() for _ in Base.OneTo(ntraj)]
            ys=[au.aPBCy[1]*rand() for _ in Base.OneTo(ntraj)]
        else
            xs=(3:2:21)u"Å"
            ys=xs
        end
    end

    vary="T_Ei_xy"
    h=["T", "Ei", "xNOi", "yNOi", "xNf", "yNf", "zNf", "xOf", "yOf", "zOf", "vxNf", "vyNf", "vzNf", "vxOf", "vyOf", "vzOf", "KEtot", "Etrans", "Erot", "KEvib"]
    trajtrap=DataFrame([name => [] for name in h])
    trajscatter=DataFrame([name => [] for name in h])
    h2=["T", "Ei", "n_scatter", "n_trap", "frac_scatter", "frac_trap"]
    counttraj=DataFrame([name => [] for name in h2])

    desc = isaac ? "$(noaurundesc)_multi_runs-ISAAC" : "$(noaurundesc)_multi_runs"
    outpath=makeresultsfolder(desc,vary)

    for i in eachindex(ts)
        param.T[1]=ts[i]
        T=Int64(ustrip(u"K",param.T[1]))
        global aurundesc="Au_slab-T $T"
        runAuSlabEquilibration()
        tpath=mkpath("$outpath/T $T")

        for j in eachindex(eis)
            no.Et_i[1]=eis[j]
            ei=round(ustrip(u"e_MD",no.Et_i[1]);digits=2)
            eipath=mkpath("$tpath/Ei $ei")

            nscatter=0
            ntrap=0
            for k in eachindex(xs)
                xNOi=xs[k]
                x=round(ustrip(u"Å",xNOi);digits=3)
                yNOi=ys[k]
                y=round(ustrip(u"Å",yNOi);digits=3)
                xypath=makeresultsfolder("$eipath/x $x y $y")

                sys=runNOAuTrajectory(xNOi,yNOi,xypath)
                if !checkEconserved(sys)
                    error("Energy not conserved")
                end

                finalNOinfo=finalE_NO(sys)
                allinfo=vcat([T, ei, x, y], finalNOinfo)
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

    # zip all folders
    zipfolders(outpath)
    outputmultirunsummary(vary,counttraj,outpath)
    outputtrajinfo(trajscatter,trajtrap,outpath)
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
