# setup and running of NO/Au scattering

############################################################################################################

# setup

# defining MD propagation method (velocity verlet)
simulator_NOAu = VelocityVerlet(
    # Time step
    dt=param.dt[1],

    # dont remove center of mass motion to keep layer fixed. may revert.
    remove_CM_motion=false,

    # # random scaling of atom velocities for thermal equilibration? setting time constant to 500*dt (same as Molly example). may change later
    # coupling=AndersenThermostat(param.T[1], 500*param.dt[1]),
)

# defining system
sys_NOAu = System(
    # initializing atoms in system
    atoms=initNOAuAtoms(),

    # system bound by custom NO/Au interactions. using neighbor list=true
    pairwise_inters=(NOAuInteraction(true),),

    # initial atom coordinates. using static arrays (SA) for Molly compatibility
    coords=initNOAuCoords(),

    # initial atom velocities. NO vel based on √(2E/m). Au slab set to 0
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

# # may use instead of r0ij, rAu in future
# initcoords=copy(sys_Au.coords)

# tmp step counter
step_no=1

# running MD + output results
runMDprintresults(sys_NOAu, noaurundesc, simulator_NOAu, steps_dyn)
