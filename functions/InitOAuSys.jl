# functions for setting up O/Au system in Molly.

############################################################################################################

function initOAuCoords(xOi=au.aPBCx[1]*rand(),yOi=au.aPBCy[1]*rand())
    zO=12u"Å"+maximum(au.z)
    ocoords=[SA[xOi,yOi,zO]]
    aucoords=getEquilAuCoords()
    vcat(ocoords,aucoords)
end

############################################################################################################

function initOAuAtoms()
    oatom=[Atom(index=1, mass=no.mO[1])]
    auatoms=[Atom(index=i, mass=au.m[1]) for i in 2:au.N[1]+1]
    vcat(oatom,auatoms)
end

############################################################################################################

function initOVelocities(Ei=no.Et_i[1])
    # convert incident molecule energy to velocity
    θ=no.θi[1]
    mass=no.mO[1]*N_A # now in kg/mol
    vmag=sqrt(2*Ei/mass)
    [SA[vmag*sin(θ),0u"m/s",-vmag*cos(θ)]] # negative: pointing down
end

############################################################################################################

function initOAuVelocities(Ei=no.Et_i[1])
    vo=initOVelocities(Ei)
    auv=getLastAuVelocities()
    vcat(vo,auv)
end

############################################################################################################

"""
initiallize O/Au system. O placed at specified x,y position
"""
function initOAuSys(xcom::Unitful.Length=au.aPBCx[1]*rand(),
                            ycom::Unitful.Length=au.aPBCy[1]*rand(),
                            Ei::EnergyPerMole=no.Et_i[1])
    # Molly sys params
    atoms=initOAuAtoms()
    pairwise_inters=(OAuInteraction(true),) # using neighbor list=true
    coords=initOAuCoords(xcom,ycom)
    velocities=initOAuVelocities(Ei) # O vel based on √(2Ei/m). Au slab set to last vels
    boundary=simboxdims # periodic in x,y
    neighbor_finder=ONeighborFinder(O_cutoff=PES_GS.AuOcutoff[1],)
    if simplerun                                          
        loggers=(
                    et=TotalEnergyLogger(steps_dyn), # checking energy conservation AT LAST TIME STEP ONLY
                    pe=PotentialEnergyLogger(steps_dyn),
                    ke=KineticEnergyLogger(steps_dyn),
                    velocities=VelocityLogger(steps_dyn),
                    charge=ChargeLogger(actsteplog),
                    coords=CoordinateLogger(actsteplog),
                    )
    else
        loggers=(
            et=TotalEnergyLogger(actsteplog), # checking energy conservation 
            pe=PotentialEnergyLogger(actsteplog),
            ke=KineticEnergyLogger(actsteplog),
            velocities=VelocityLogger(steps_dyn),
            forces=ForceLogger(steps_dyn), # E needs to be calculated before F
            charge=ChargeLogger(actsteplog),
            coords=CoordinateLogger(actsteplog),
        )
    end
    dt=param.dt[1]
    remove_CM_motion=false # dont remove center of mass motion to keep layer fixed. may revert.

    # defining system/simulator
    s = System(
                atoms=atoms,
                pairwise_inters=pairwise_inters,
                coords=coords,
                velocities=velocities,
                boundary=boundary,
                neighbor_finder=neighbor_finder,
                loggers=loggers,
                )
    simul = VelocityVerlet(
                                    dt=dt,
                                    remove_CM_motion=remove_CM_motion,
                                )

    return s, simul
end

############################################################################################################

"""
run o/au trajectory and output run info to results folder
"""
function runOAuTrajectory(x::Unitful.Length=au.aPBCx[1]*rand(),
                                    y::Unitful.Length=au.aPBCy[1]*rand(),
                                    T::Unitful.Temperature=param.T[1],
                                    Ei::EnergyPerMole=no.Et_i[1],
                                    path::String=makeresultsfolder(oaurundesc,steps_dyn_OAu)
                                    )
    # redefine for O/Au
    global auatomcutoff=398

    # initialize system
    sys_OAu, simulator_OAu = initOAuSys(x,y,Ei)

    # running MD + output results
    t=@elapsed runMDprintresults(sys_OAu, oaurundesc, simulator_OAu, steps_dyn_OAu, path, T, Ei)
    checkEconserved(sys_OAu)

    println("O/Au trajectory is complete")
    println("Conditions: T=$T, Ei=$Ei, x=$x, y=$y")
    println("Time to run: $t seconds")
    println()

    # back to original
    global auatomcutoff=399

    return sys_OAu
end
