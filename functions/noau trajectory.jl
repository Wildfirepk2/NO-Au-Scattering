# setup and running of NO/Au scattering

############################################################################################################

# main variables

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

# slab potential

############################################################################################################

# defining MD propagation method (velocity verlet)
simulator_NOAu = VelocityVerlet(
    # Time step
    dt=param.dt[1],

    # dont remove center of mass motion to keep layer fixed. may revert.
    remove_CM_motion=false,
)

# defining system
sys_NOAu = System(
    # initializing atoms in system
    atoms=initNOAuAtoms(),

    # system bound by custom NO/Au interactions. using neighbor list=true
    pairwise_inters=(NOAuInteraction(true),),

    # initial atom coordinates. using static arrays (SA) for Molly compatibility
    coords=initNOAuCoords(),

    # initial atom velocities. NO vel based on √(2E/m). Au slab set to last vels
    velocities=initNOAuVelocities(),

    # system boundary. is periodic in x,y
    boundary=simboxdims,

    # using custom neighbor finder
    neighbor_finder=NONeighborFinder(
        N_cutoff=PES_GS.AuNcutoff[1],
        O_cutoff=PES_GS.AuOcutoff[1],
    ),

    # tracking parameters wrt time. value in parentheses is number of time steps. log at last step: set to steps_dyn, default steps: set to actsteplog
    loggers=(
        # checking energy conservation. E needs to be calculated before F
        et=TotalEnergyLogger(actsteplog),
        pe=PotentialEnergyLogger(actsteplog),
        ke=KineticEnergyLogger(actsteplog),

        # capture velocities and forces at last time step
        velocities=VelocityLogger(steps_dyn),
        forces=ForceLogger(steps_dyn),

        # for animation
        coords=CoordinateLogger(actsteplog),
    ),
)

# running MD + output results
runMDprintresults(sys_NOAu, noaurundesc, simulator_NOAu, steps_dyn)
