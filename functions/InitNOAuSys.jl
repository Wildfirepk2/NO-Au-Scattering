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
		readcoorddir=readdir("$audir/syscoords")
		lastcoords=readcoorddir[end]
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
# # \old
# function initNOCoords()
#     r=no.r[1]
#     n_molec=1
#     nocoords=place_diatomics(n_molec,virtboxdims,r)

# 	# place N 12 above Au surface
# 	placeat=12u"Å"+maximum(au.z)
# 	nz=nocoords[1][3]
# 	delta=placeat-nz
# 	[nocoords[i]+[0u"Å",0u"Å",delta] for i in eachindex(nocoords)]
# end

############################################################################################################

# \debugging
function initNOCoords(xNi=au.aPBCx[1]*rand(),yNi=au.aPBCy[1]*rand())
    r=no.r[1]
	zN=12u"Å"+maximum(au.z)
	zO=zN+r

	# place N and O at same spot
	[SA[xNi,yNi,zN],SA[xNi,yNi,zO]]
end

############################################################################################################

function initNOAuCoords(xNi=au.aPBCx[1]*rand(),yNi=au.aPBCy[1]*rand())
	nocoords=initNOCoords(xNi,yNi)
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

function initNOAuVelocities()
	uv=normalize(SVector{3}(rand(-10:10,3))) # unit vector
	uvneg=uv[3]<0 ? uv : -uv # downward velocity. moving toward surface

	# convert incident molecule energy to velocity
	e=no.Et_i[1]
	mass=(no.mN[1]+no.mO[1])*N_A # now in kg/mol
	v=uvneg*sqrt(2*e/mass)
	
	nov=[v,v]
	auv=[velocity(1u"u", 0u"K") for _ in 1:au.N[1]]
	vcat(nov,auv)

	# \debugging
	t=[SA[0u"Å/ps",0u"Å/ps",-sqrt(2*e/mass)],SA[0u"Å/ps",0u"Å/ps",-sqrt(2*e/mass)]]
	auvt=getLastAuVelocities()
	vcat(t,auvt)
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
initiallize NO/Au system. N placed at specified x,y position
"""
function initNOAuSys(xNi=au.aPBCx[1]*rand(),yNi=au.aPBCy[1]*rand())
	# defining MD propagation method (velocity verlet)
	simul = VelocityVerlet(
		# Time step
		dt=param.dt[1],

		# dont remove center of mass motion to keep layer fixed. may revert.
		remove_CM_motion=false,
	)

	# defining system
   if multirun
      s = System(
         # initializing atoms in system
         atoms=initNOAuAtoms(),
   
         # system bound by custom NO/Au interactions. using neighbor list=true
         pairwise_inters=(NOAuInteraction(true),),
   
         # initial atom coordinates. using static arrays (SA) for Molly compatibility
         coords=initNOAuCoords(xNi,yNi),
   
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
            # checking energy conservation AT LAST TIME STEP ONLY. E needs to be calculated before F
            et=TotalEnergyLogger(steps_dyn),
            pe=PotentialEnergyLogger(steps_dyn),
            ke=KineticEnergyLogger(steps_dyn),
   
            # capture ONLY velocities at last time step
            velocities=VelocityLogger(steps_dyn),
   
            # for animation
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
         coords=initNOAuCoords(xNi,yNi),
   
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
   end

	return s, simul
end

############################################################################################################

"""
run no/au trajectory and output run info to results folder
"""
function runNOAuTrajectory(xNi=au.aPBCx[1]*rand(),yNi=au.aPBCy[1]*rand(),path::String=makeresultsfolder(noaurundesc,steps_dyn))
	# initialize system
	sys_NOAu, simulator_NOAu = initNOAuSys(xNi,yNi)

	# running MD + output results
	t=@elapsed runMDprintresults(sys_NOAu, noaurundesc, simulator_NOAu, steps_dyn, path)

	println("NO/Au trajectory is complete")
	println("Time to run: $t seconds")

	return sys_NOAu
end
