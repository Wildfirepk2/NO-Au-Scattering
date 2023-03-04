# functions for setting up NO/Au system in Molly.

############################################################################################################

"""
get equilibrated au coords from previous run
"""
function getEquilAuCoords()
    audir=getAuDirPath("results")
    coordsfile="$audir/syscoords.xlsx"
    xfcoord=XLSX.readxlsx(coordsfile)
    sheets=XLSX.sheetnames(xfcoord)
    sheetlastcoord=sheets[end]
    dfcoord=DataFrame(XLSX.readtable(coordsfile,sheetlastcoord))
    [SA[dfcoord[i,1],dfcoord[i,2],dfcoord[i,3]]u"d_MD" for i in 1:nrow(dfcoord)]
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

function initNOCoords()
    r=no.r[1]
    n_molec=1
    nocoords=place_diatomics(n_molec,virtboxdims,r)

	# place N 12 above Au surface
	placeat=12u"Å"+maximum(au.z)
	nz=nocoords[1][3]
	delta=placeat-nz
	[nocoords[i]+[0u"Å",0u"Å",delta] for i in eachindex(nocoords)]

	# \debugging
	zN=12u"Å"+maximum(au.z)
	zO=zN+r
	[SA[3u"Å",3u"Å",zN],SA[3u"Å",3u"Å",zO]]
end

############################################################################################################

function initNOAuCoords()
	nocoords=initNOCoords()
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
run no/au trajectory and output run info to results folder
"""
function runNOAuTrajectory()
    t=@elapsed include("functions/noau trajectory.jl")
    println("NO/Au trajectory is complete")
    println("Time to run: $t seconds")
end
