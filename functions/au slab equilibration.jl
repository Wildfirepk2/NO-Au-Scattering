# setup and running of Au slab equilibration

############################################################################################################

# main variables

# 3 x N matrix. xyz coords of each Au atom stored in columns
rAu=initAuCoords()

# Molly object. simulation box dimensions. periodic in x,y directions. needed in other areas outside Molly
simboxdims=CubicBoundary(au.aPBCx[1], au.aPBCy[1], au.aPBCz[1])

# nn: array of arrays. nearest neighbors for each Au atom. ith row corresponds to Au atom i's nearest neighbors (in terms of atom number)
# nn_molly: molly neighbor object. same as nn. for molly compatibility
nn,nn_molly=getnn()

# array of matrices. initial distances to nearest neighbors for each atom. nn xyz coords stored in columns
r0ij=getdrij(rAu)

# array of arrays. case number for nearest neighbors for each atom. 
nncase=getcase()

# Dict mapping nn case number to force matrix.
Aij=initAij()

# array of arrays of matrices. force matrix for nearest neighbors for each atom. 
Aijarray=initAijarray()

############################################################################################################

# Au slab equilibration (MD) using Molly

# defining MD propagation method (velocity verlet)
simulator = VelocityVerlet(
    # Time step
    dt=param.dt[1],

    # dont remove center of mass motion to keep layer fixed. may revert.
    remove_CM_motion=false,

    # # random scaling of atom velocities for thermal equilibration? setting time constant to 500*dt (same as Molly example). may change later
    # coupling=AndersenThermostat(param.T[1], 500*param.dt[1]),
)

# defining system
sys_Au = System(
    # initializing atoms in system
    atoms=[Atom(index=i, mass=au.m[1]) for i in 1:au.N[1]],

    # system bound by custom Au slab interactions. using neighbor list=true
    pairwise_inters=(AuSlabInteraction(true),),

    # initial atom coordinates. using static arrays (SA) for Molly compatibility
    coords=[SA[au.x[i],au.y[i],au.z[i]] for i in 1:au.N[1]],

    # initial atom velocities based on maxwell-Boltzmann distribution at system temp. freezing back layer (velocity at 0K is 0). using velocity function for back layer for consistent units
    velocities=[i<auatomcutoff ? velocity(au.m[1], param.T[1]) : velocity(1u"u", 0u"K") for i in 1:au.N[1]],

    # system boundary. is periodic in x,y
    boundary=simboxdims,

    # using custom neighbor finder
    neighbor_finder=AuNeighborFinder(),

    # tracking parameters wrt time. value in parentheses is number of time steps. log at last step: set to steps_eq, default steps: set to actsteplog
    loggers=(
        # capture velocities and forces at last time step
        velocities=VelocityLogger(steps_eq),
        forces=ForceLogger(steps_eq),

        # checking energy conservation
        et=TotalEnergyLogger(actsteplog),
        pe=PotentialEnergyLogger(actsteplog),
        ke=KineticEnergyLogger(actsteplog),

        # for animation
        coords=CoordinateLogger(actsteplog),
    ),
)

# # may use instead of r0ij, rAu in future
# initcoords=copy(sys_Au.coords)

# tmp step counter
step_no=1

# running MD + output results
AuSlabEquilibration(sys_Au, aurundesc, simulator, steps_eq)
