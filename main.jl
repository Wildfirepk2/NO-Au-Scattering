# adiabatic molecular dynamics (MD) of NO scattering off of an Au(111) surface

# checked: 10/27/22
############################################################################################################

# initialize all functions
include("functions/InitMDUnits.jl")
include("functions/InitAllParams.jl")
include("functions/InitAllAuNOVFunc.jl")
include("functions/InitAllAuVFunc.jl")
include("functions/InitCustomMollyFunc.jl")

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

# main variables

# 3 x N matrix. xyz coords of each Au atom stored in columns
rAu=initAuCoords()

# array of arrays. nearest neighbors for each Au atom. ith row corresponds to Au atom i's nearest neighbors (in terms of atom number)
nn=getnn()

# array of matrices. initial distances to nearest neighbors for each atom. nn xyz coords stored in columns
r0ij=getdrij(rAu)

# array of arrays. case number for nearest neighbors for each atom. 
nncase=getcase()

# Dict mapping nn case number to force matrix.
Aij=initAij()

# array of arrays of matrices. force matrix for nearest neighbors for each atom. 
Aijarray=initAijarray()

############################################################################################################

# other variables
Av_Et = 0
Av_theta = 0
Av_weighted_theta = 0
Norm_weighted_theta = 0
Av_Evib = 0
Av_Er = 0
Av_EAu = 0
NTRAP = 0

# first Au atom of last layer. last layer (atoms 397-528) is frozen. may put in au var
auatomcutoff=397

# log parameters after n steps. may put in param var
stepslogging=10

############################################################################################################

# Au slab equilibration using Molly

# defining MD propagation method (velocity verlet)
simulator = VelocityVerlet(
    # Time step
    dt=param.dt[1],

    # # random scaling of atom velocities for thermal equilibration? setting time constant to 500*dt (same as Molly example). may change later
    # coupling=AndersenThermostat(param.T[1], 500*param.dt[1]),
)

# defining system
sys_Au = System(
    # initializing atoms in system
    atoms=[Atom(index=i, mass=au.m[1]) for i in 1:au.N[1]],

    # system bound by custom Au slab interactions 
    pairwise_inters=(AuSlabInteraction(false),),

    # initial atom coordinates. using static arrays (SA) for Molly compatibility
    coords=[SA[au.x[i],au.y[i],au.z[i]] for i in 1:au.N[1]],

    # initial atom velocities based on maxwell-Boltzmann distribution at system temp. freezing back layer (velocity at 0K is 0). using velocity function for back layer for consistent units
    velocities=[i<auatomcutoff ? velocity(au.m[1], param.T[1]) : velocity(1u"u", 0u"K") for i in 1:au.N[1]],

    # system boundary. is periodic
    boundary=CubicBoundary(au.aPBCx[1], au.aPBCy[1], au.aPBCz[1]),

    # tracking parameters wrt time. value in parentheses is number of time steps
    loggers=(
        # capture velocities and forces at last time step
        velocity=VelocityLogger(param.Nsteps_eq[1]),
        forces=ForceLogger(param.Nsteps_eq[1]),

        # # checking energy conservation
        # te=TotalEnergyLogger(stepslogging),
        # pe=PotentialEnergyLogger(stepslogging),
        # ke=KineticEnergyLogger(stepslogging),

        # for animation
        coords=CoordinateLogger(stepslogging),
    ),
)

# # may use instead of r0ij, rAu in future
# initcoords=copy(sys_Au.coords)

# tmp step counter
step_no=1

# run MD. cutting down steps for debugging
simulate!(sys_Au, simulator, param.Nsteps_eq[1]/1000)

# output all system data: animation, coords, last velocities/forces
outputsysinfo(sys_Au,"au slab equilibration")
