# functions for setting up NO/Au system in Molly.

############################################################################################################

"""
get equilibrated au coords from previous run
"""
function getEquilAuCoords()
   audir=getAuDirPath("results")
	readaudir=readdir(audir)
	i_coords=findfirst(contains.(readaudir,"syscoords.xlsx"))
	if i_coords isa Int
		coordsfile="$audir/syscoords.xlsx"
		xfcoord=XLSX.readxlsx(coordsfile)
		sheets=XLSX.sheetnames(xfcoord)
		sheetlastcoord=sheets[end]
		dfcoord=DataFrame(XLSX.readtable(coordsfile,sheetlastcoord))
		[SA[dfcoord[i,1],dfcoord[i,2],dfcoord[i,3]]u"d_MD" for i in 1:nrow(dfcoord)]
	else
		readcoorddir=readdir("$audir/syscoords") # incorrectly sorts by number. eg 1,10,11,12,...,2,etc
      sortedcoorddir = sort(readcoorddir, by = x -> parse(Int, match(r"\d+", x).match)) # sorted by name. eg syscoords 1,2,3,4,etc
		lastcoords=sortedcoorddir[end]
		dfcoord=CSV.read("$audir/syscoords/$lastcoords",DataFrame)
		[SA[dfcoord[i,1],dfcoord[i,2],dfcoord[i,3]]u"d_MD" for i in 1:nrow(dfcoord)]
	end
end

############################################################################################################

"""
get last au velocities from previous run
"""
function getLastAuVelocities()
    audir=getAuDirPath("results")
    velfile="$audir/sysvelocities.xlsx"
    xfvel=XLSX.readxlsx(velfile)
    sheets=XLSX.sheetnames(xfvel)
    sheetlastcoord=sheets[end]
    dfvel=DataFrame(XLSX.readtable(velfile,sheetlastcoord))
    [SA[dfvel[i,1],dfvel[i,2],dfvel[i,3]]u"v_MD" for i in 1:nrow(dfvel)]
end

############################################################################################################

"""initial NO bond length/vib velocity"""
function getrvNO()
   idx=rand(1:length(no.r))
   r=no.r[idx]
   v=no.v[idx]
   return r,v
end

############################################################################################################

"""random initial NO orientation"""
function getNOorient()
   vec=2rand(3).-1
   normalize(vec)
end

"""
specified initial NO orientation

θ: angle between NO bond and surface normal (0°: N pointing down)
"""
function getNOorient(θ::Unitful.DimensionlessQuantity)
   [sin(θ),0,cos(θ)]
end

############################################################################################################

function initNOCoords(rNOi,u::Vector,xcom=au.aPBCx[1]*rand(),ycom=au.aPBCy[1]*rand(),zcom=12u"Å")
   # masses
   mN=no.mN[1]
   mO=no.mO[1]
   μ=mN*mO/(mN+mO)

   zcomact=zcom+maximum(au.z)
   com=SA[xcom,ycom,zcomact]

   rN=com-rNOi*μ/mN*u
   rO=com+rNOi*μ/mO*u
   [rN,rO]
end

############################################################################################################

function initNOAuCoords(rNOi,u::Vector,xcom=au.aPBCx[1]*rand(),ycom=au.aPBCy[1]*rand())
	nocoords=initNOCoords(rNOi,u,xcom,ycom)
	aucoords=getEquilAuCoords()
	vcat(nocoords,aucoords)
end

############################################################################################################

function initNOAuAtoms()
	noatoms=[Atom(index=1, mass=no.mN[1]), Atom(index=2, mass=no.mO[1])]
	auatoms=[Atom(index=i, mass=au.m[1]) for i in 3:au.N[1]+2]
	vcat(noatoms,auatoms)
end

############################################################################################################

function initNOVelocities(vrel,u::Vector,Ei=no.Et_i[1])
   # masses
   mN=no.mN[1]
   mO=no.mO[1]
   mt=mN+mO
   μ=mN*mO/mt

   # convert incident molecule energy to velocity
   θ=no.θi[1]
	mass=mt*N_A # now in kg/mol
	vmag=-sqrt(2*Ei/mass)
   vcom=SA[vmag*sin(θ),0u"m/s",vmag*cos(θ)]

   vN=vcom-μ/mN*vrel*u
   vO=vcom+μ/mO*vrel*u
   [vN,vO]
end

############################################################################################################

function initNOAuVelocities(vrel,u::Vector,Ei=no.Et_i[1])
	vno=initNOVelocities(vrel,u,Ei)
   auv=getLastAuVelocities()
	vcat(vno,auv)
end

############################################################################################################

"""
helper function: get path of au folder in directory. if not found, return nothing
"""
function getAuDirPath(path::String=".")
    resultsinpath=readdir(path;join=true)
    i_au=findfirst(contains.(resultsinpath,aurundesc))
    i_au isa Nothing ? i_au : resultsinpath[i_au]
end

############################################################################################################

"""
initiallize NO/Au system. NO placed at specified x,y position
"""
function initNOAuSys(xcom::Unitful.Length=au.aPBCx[1]*rand(),
                     ycom::Unitful.Length=au.aPBCy[1]*rand(),
                     Ei::EnergyPerMole=no.Et_i[1]
                     )
   # initial NO BL/vib vel
   rNOi,vrel=getrvNO()

   # initial NO orientation
   u=getNOorient()

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
         atoms=initNOAuAtoms(),
   
         # system bound by custom NO/Au interactions. using neighbor list=true
         pairwise_inters=(NOAuInteraction(true),),
   
         # initial atom coordinates. using static arrays (SA) for Molly compatibility
         coords=initNOAuCoords(rNOi,u,xcom,ycom),
   
         # initial atom velocities. NO vel based on √(2E/m). Au slab set to last vels
         velocities=initNOAuVelocities(vrel,u,Ei),
   
         # system boundary. is periodic in x,y
         boundary=simboxdims,
   
         # using custom neighbor finder
         neighbor_finder=NONeighborFinder(
            N_cutoff=PES_GS.AuNcutoff[1],
            O_cutoff=PES_GS.AuOcutoff[1],
         ),
   
         # tracking parameters wrt time. value in parentheses is number of time steps. log at last step: set to steps_dyn, default steps: set to actsteplog
         loggers=(
            # checking energy conservation AT LAST TIME STEP ONLY. E needs to be calculated before F
            et=TotalEnergyLogger(steps_dyn),
            pe=PotentialEnergyLogger(steps_dyn),
            ke=KineticEnergyLogger(steps_dyn),
   
            # capture ONLY velocities at last time step
            velocities=VelocityLogger(steps_dyn),
   
            # for animation
            charge=ChargeLogger(actsteplog),
            coords=CoordinateLogger(actsteplog),
         ),
      )
   else
      s = System(
         # initializing atoms in system
         atoms=initNOAuAtoms(),
   
         # system bound by custom NO/Au interactions. using neighbor list=true
         pairwise_inters=(NOAuInteraction(true),),
   
         # initial atom coordinates. using static arrays (SA) for Molly compatibility
         coords=initNOAuCoords(rNOi,u,xcom,ycom),
   
         # initial atom velocities. NO vel based on √(2E/m). Au slab set to last vels
         velocities=initNOAuVelocities(vrel,u),
   
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
            charge=ChargeLogger(actsteplog),
            coords=CoordinateLogger(actsteplog),
         ),
      )
   end

	return s, simul
end

"""
initiallize NO/Au system. NO placed at specified x,y position AND specified orientation
"""
function initNOAuSys(θorient::Unitful.DimensionlessQuantity,
                     xcom::Unitful.Length=au.aPBCx[1]*rand(),
                     ycom::Unitful.Length=au.aPBCy[1]*rand(),
                     Ei::EnergyPerMole=no.Et_i[1]
                     )
   # initial NO BL/vib vel
   rNOi,vrel=getrvNO()

   # initial NO orientation
   u=getNOorient(θorient)

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
         atoms=initNOAuAtoms(),
   
         # system bound by custom NO/Au interactions. using neighbor list=true
         pairwise_inters=(NOAuInteraction(true),),
   
         # initial atom coordinates. using static arrays (SA) for Molly compatibility
         coords=initNOAuCoords(rNOi,u,xcom,ycom),
   
         # initial atom velocities. NO vel based on √(2E/m). Au slab set to last vels
         velocities=initNOAuVelocities(vrel,u,Ei),
   
         # system boundary. is periodic in x,y
         boundary=simboxdims,
   
         # using custom neighbor finder
         neighbor_finder=NONeighborFinder(
            N_cutoff=PES_GS.AuNcutoff[1],
            O_cutoff=PES_GS.AuOcutoff[1],
         ),
   
         # tracking parameters wrt time. value in parentheses is number of time steps. log at last step: set to steps_dyn, default steps: set to actsteplog
         loggers=(
            # checking energy conservation AT LAST TIME STEP ONLY. E needs to be calculated before F
            et=TotalEnergyLogger(steps_dyn),
            pe=PotentialEnergyLogger(steps_dyn),
            ke=KineticEnergyLogger(steps_dyn),
   
            # capture ONLY velocities at last time step
            velocities=VelocityLogger(steps_dyn),
   
            # for animation
            charge=ChargeLogger(actsteplog),
            coords=CoordinateLogger(actsteplog),
         ),
      )
   else
      s = System(
         # initializing atoms in system
         atoms=initNOAuAtoms(),
   
         # system bound by custom NO/Au interactions. using neighbor list=true
         pairwise_inters=(NOAuInteraction(true),),
   
         # initial atom coordinates. using static arrays (SA) for Molly compatibility
         coords=initNOAuCoords(rNOi,u,xcom,ycom),
   
         # initial atom velocities. NO vel based on √(2E/m). Au slab set to last vels
         velocities=initNOAuVelocities(vrel,u),
   
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
            charge=ChargeLogger(actsteplog),
            coords=CoordinateLogger(actsteplog),
         ),
      )
   end

	return s, simul
end

############################################################################################################

"""
run no/au trajectory and output run info to results folder
"""
function runNOAuTrajectory(xcom::Unitful.Length=au.aPBCx[1]*rand(),
                           ycom::Unitful.Length=au.aPBCy[1]*rand(),
                           T::Unitful.Temperature=param.T[1],
                           Ei::EnergyPerMole=no.Et_i[1],
                           path::String=makeresultsfolder(noaurundesc,steps_dyn)
                           )
	# initialize system
	sys_NOAu, simulator_NOAu = initNOAuSys(xcom,ycom,Ei)

	# running MD + output results
	t=@elapsed runMDprintresults(sys_NOAu, noaurundesc, simulator_NOAu, steps_dyn, T, Ei, path)
   checkEconserved(sys_NOAu)

	println("NO/Au trajectory is complete")
   println("Conditions: T=$T, Ei=$Ei, xcom=$xcom, ycom=$ycom")
	println("Time to run: $t seconds")
   println()

	return sys_NOAu
end

function runNOAuTrajectory(θorient::Unitful.DimensionlessQuantity,
                           xcom::Unitful.Length=au.aPBCx[1]*rand(),
                           ycom::Unitful.Length=au.aPBCy[1]*rand(),
                           T::Unitful.Temperature=param.T[1],
                           Ei::EnergyPerMole=no.Et_i[1],
                           path::String=makeresultsfolder(noaurundesc,steps_dyn)
                           )
	# initialize system
	sys_NOAu, simulator_NOAu = initNOAuSys(θorient,xcom,ycom,Ei)

	# running MD + output results
	t=@elapsed runMDprintresults(sys_NOAu, noaurundesc, simulator_NOAu, steps_dyn, T, Ei, path)
   checkEconserved(sys_NOAu)

	println("NO/Au trajectory is complete")
   println("Conditions: T=$T, Ei=$Ei, θorient=$θorient, xcom=$xcom, ycom=$ycom")
	println("Time to run: $t seconds")
   println()

	return sys_NOAu
end