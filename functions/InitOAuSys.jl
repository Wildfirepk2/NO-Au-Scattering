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

function initOAuVelocities()
	uv=normalize(SVector{3}(rand(-10:10,3))) # unit vector
	uvneg=uv[3]<0 ? uv : -uv # downward velocity. moving toward surface

	# convert incident molecule energy to velocity
	e=no.Et_i[1]
	mass=no.mO[1]*N_A # now in kg/mol
	v=uvneg*sqrt(2*e/mass)
	
	nov=[v,v]
	auv=[velocity(1u"u", 0u"K") for _ in 1:au.N[1]]
	vcat(nov,auv)

	# \debugging
	t=[SA[0u"Å/ps",0u"Å/ps",-sqrt(2*e/mass)]]
	auvt=getLastAuVelocities()
	vcat(t,auvt)
end

############################################################################################################

"""
initiallize O/Au system. O placed at specified x,y position
"""
function initOAuSys(xOi=au.aPBCx[1]*rand(),yOi=au.aPBCy[1]*rand())
	# defining MD propagation method (velocity verlet)
	simul = VelocityVerlet(
		# Time step
		dt=param.dt[1],

		# dont remove center of mass motion to keep layer fixed. may revert.
		remove_CM_motion=false,
	)

	# defining system
   if simplerun
      s = System(
         # initializing atoms in system
         atoms=initOAuAtoms(),
   
         # system bound by custom O/Au interactions. using neighbor list=true
         pairwise_inters=(OAuInteraction(true),),
   
         # initial atom coordinates. using static arrays (SA) for Molly compatibility
         coords=initOAuCoords(xOi,yOi),
   
         # initial atom velocities. NO vel based on √(2E/m). Au slab set to last vels
         velocities=initOAuVelocities(),
   
         # system boundary. is periodic in x,y
         boundary=simboxdims,
   
         # using custom neighbor finder
         neighbor_finder=ONeighborFinder(
            O_cutoff=PES_GS.AuOcutoff[1],
         ),
   
         # tracking parameters wrt time. value in parentheses is number of time steps. log at last step: set to steps_dyn_OAu, default steps: set to actsteplog
         loggers=(
            # checking energy conservation AT LAST TIME STEP ONLY. E needs to be calculated before F
            et=TotalEnergyLogger(steps_dyn_OAu),
            pe=PotentialEnergyLogger(steps_dyn_OAu),
            ke=KineticEnergyLogger(steps_dyn_OAu),
   
            # capture ONLY velocities at last time step
            velocities=VelocityLogger(steps_dyn_OAu),
   
            # for animation
            coords=CoordinateLogger(actsteplog),
         ),
      )
   else
      s = System(
         # initializing atoms in system
         atoms=initOAuAtoms(),
   
         # system bound by custom O/Au interactions. using neighbor list=true
         pairwise_inters=(OAuInteraction(true),),
   
         # initial atom coordinates. using static arrays (SA) for Molly compatibility
         coords=initOAuCoords(xOi,yOi),
   
         # initial atom velocities. NO vel based on √(2E/m). Au slab set to last vels
         velocities=initOAuVelocities(),
   
         # system boundary. is periodic in x,y
         boundary=simboxdims,
   
         # using custom neighbor finder
         neighbor_finder=ONeighborFinder(
            O_cutoff=PES_GS.AuOcutoff[1],
         ),
   
         # tracking parameters wrt time. value in parentheses is number of time steps. log at last step: set to steps_dyn_OAu, default steps: set to actsteplog
         loggers=(
            # checking energy conservation. E needs to be calculated before F
            et=TotalEnergyLogger(actsteplog),
            pe=PotentialEnergyLogger(actsteplog),
            ke=KineticEnergyLogger(actsteplog),
   
            # capture velocities and forces at last time step
            velocities=VelocityLogger(steps_dyn_OAu),
            forces=ForceLogger(steps_dyn_OAu),
   
            # for animation
            coords=CoordinateLogger(actsteplog),
         ),
      )
   end

	return s, simul
end

############################################################################################################

"""
run o/au trajectory and output run info to results folder
"""
function runOAuTrajectory(xOi=au.aPBCx[1]*rand(),yOi=au.aPBCy[1]*rand(),path::String=makeresultsfolder(oaurundesc,steps_dyn_OAu))
   # redefine for O/Au
	global auatomcutoff=398

	# initialize system
	sys_OAu, simulator_OAu = initOAuSys(xOi,yOi)

	# running MD + output results
	t=@elapsed runMDprintresults(sys_OAu, oaurundesc, simulator_OAu, steps_dyn_OAu, path)

	println("O/Au trajectory is complete")
	println("Time to run: $t seconds")

   # back to original
   global auatomcutoff=399

	return sys_OAu
end