# functions for initializing variables from input files. 
# variables are DataFrames. 
# some DataFrames' columns have only 1 entry so must be referenced by VAR.COL[1]

# checked: 11/1/22 (including all input files)

############################################################################################################

function initAuParams()
	# initialize Au parameters from input file
	file="input files/au params.csv"
	au=CSV.read(file,DataFrame)

	# add units. assuming A from unit system. see notebook for slab info
	au.x*=u"Å"
	au.y*=u"Å"
	au.z*=u"Å"
	au.dnn*=u"Å"
	au.a*=u"Å"
	au.aPBCx*=u"Å" # 1 layer of space beyond slab boundary for periodic boundary conditions
	au.aPBCy*=u"Å" # 1 layer of space beyond slab boundary for periodic boundary conditions
	au.aPBCz*=u"Å"
	au.m*=u"u"

	# shift up atoms such that back layer at y=0. needed for Molly compatibility/giving 1 layer of space
	au.y.-=minimum(au.y)

	# shift up atoms such that back layer at z=0. needed for Molly compatibility
	au.z.-=minimum(au.z)

	# # set infinite height for system. may revert later. best for Molly
	# au.aPBCz[1]=Inf*u"Å"

	# see roy art table 3 (p8). see fortran getVNN2 for adding NA
	au.α*=u"N/m"*N_A
	au.β*=u"N/m"*N_A
	au.γ*=u"N/m"*N_A

	return au
end

############################################################################################################

function initSimParams()
	# initialize simulation parameters from input file
	file="input files/input.csv"
	input=CSV.read(file,DataFrame)

	# add units. T: assuming K from unit system. see fortran dt for unit
	input.T*=u"K"
	input.dt*=u"fs"

	return input
end

############################################################################################################

function initPESParamGS()
	# initialize ground state PES parameters from input file
	file="input files/pes params-ground.csv"
	gs=CSV.read(file,DataFrame)

	# add units. see roy art table 2 (p6). 1/A for V function (eg eq 14) compatibility
	gs.A0*=u"kJ/mol"
	gs.α0*=u"1/Å"
	gs.AuOcutoff*=u"Å"
	gs.B0*=u"kJ/mol"
	gs.β0*=u"1/Å"
	gs.AuNcutoff*=u"Å"
	gs.F0*=u"kJ/mol"
	gs.γ0*=u"1/Å"
	gs.r0_NO*=u"Å"

	return gs
end

############################################################################################################

function initPESParamIonic()
	# initialize excited state PES parameters from input file
	file="input files/pes params-ionic.csv"
	ion=CSV.read(file,DataFrame)

	# add units. see roy art table 2 (p6). p7 for ϕ and Ea eV values. 1/A for V function compatibility
	ion.A1*=u"kJ/mol"
	ion.α1*=u"1/Å"
	ion.AuOcutoff*=u"Å"
	ion.B1*=u"kJ/mol"
	ion.β1*=u"1/Å"
	ion.AuNcutoff*=u"Å"
	ion.r1_AuN*=u"Å"
	ion.D*=u"kJ/mol"
	ion.C*=u"Å"
	ion.zimage*=u"Å"
	ion.F1*=u"kJ/mol"
	ion.γ1*=u"1/Å"
	ion.r1_NO*=u"Å"

	# adding NA for compatibility with kJ/mol
	ion.ϕ*=u"eV"*N_A
	ion.Ea*=u"eV"*N_A # is negative

	return ion
end

############################################################################################################

function initPESParamCoup()
	# initialize coupling PES parameters from input file
	file="input files/pes params-coupled.csv"
	coup=CSV.read(file,DataFrame)

	# add units. see roy art table 2 (p6). 1/A for V function compatibility. A3/B3 unitless based on V function form
	coup.A2*=u"kJ/mol"
	coup.γ2*=u"1/Å"
	coup.AuOcutoff*=u"Å"
	coup.B2*=u"kJ/mol"
	coup.γ3*=u"1/Å"
	coup.AuNcutoff*=u"Å"

	return coup
end

############################################################################################################

function initNOParams()
	# initialize NO parameters from input file
	file="input files/no params.csv"
	no=CSV.read(file,DataFrame)

	# add units. r/v assumed to be same as unit system. see fortran Et_i for unit
	no.r*=u"Å"
	no.v*=u"v_MD"#check up
	no.Et_i*=u"kJ/mol"
	no.θi*=u"°"
	no.mN*=u"u"
	no.mO*=u"u"

	return no
end

############################################################################################################

"""initialize Au coordinates into 3 x 528 matrix. xyz coordinates stored in columns"""
function initAuCoords()
	transpose([au.x au.y au.z])
end

############################################################################################################

"""
get equilibrated au coords
"""
function getEquilAuCoords()
    audir=getAuDirPath("results")
    coordsfile="$audir/syscoords.xlsx"
    xfcoord=XLSX.readxlsx(coordsfile)
    sheets=XLSX.sheetnames(xfcoord)
    sheetlastcoord=sheets[end]
    dfcoord=DataFrame(XLSX.readtable(coordsfile,sheetlastcoord))
    [SA[dfcoord[i,1],dfcoord[i,2],dfcoord[i,3]] for i in 1:nrow(dfcoord)]
end

############################################################################################################

function initNOCoords()
    r=no.r[1]
    molec=1
    nocoords=place_diatomics(molec,virtboxdims,r)

	# place N 12 above Au surface
	placeat=12u"Å"+maximum(au.z)
	nz=nocoords[1][3]
	delta=placeat-nz
	[nocoords[i]+[0u"Å",0u"Å",delta] for i in eachindex(nocoords)]
	# nocoords.+[0u"Å",0u"Å",delta]
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
	uv=normalize(SVector{3}(rand(-10:10,3)))
	uvneg=uv[3]<0 ? uv : -uv
	e=no.Et_i[1]
	mass=(no.mN[1]+no.mO[1])*N_A
	v=uvneg*sqrt(2*e/mass)
	
	nov=[v,v]
	auv=[velocity(1u"u", 0u"K") for _ in 1:au.N[1]]
	vcat(nov,auv)
end