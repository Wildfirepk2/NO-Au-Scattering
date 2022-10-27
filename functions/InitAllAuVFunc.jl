# functions leading up to and including Au-Au interactions.

# checked: 10/27/22

############################################################################################################

using NearestNeighbors
using LinearAlgebra

############################################################################################################

"""
removes same atom pairs from nearest neighbor list. 

called in getnn() only.
"""	
function removeiipairs!(nn)
	for i in eachindex(nn)
		for j in eachindex(nn[i])
			if nn[i][j]==i
				splice!(nn[i],j)
				break
			end
		end
	end
end

############################################################################################################

"""
outputs 1st nearest neighbors (nn) for each Au atom as array of arrays. 

employs NearestNeighbors package. Au's crystal structure is fcc. thus, each atom (except those at the edges) has 12 nearest neighbors.
"""
function getnn()
	# stripping units for NearestNeighbors compatibility
	rAu_nounit=ustrip(rAu)
	dnn_nounit=ustrip(au.dnn[1])

	# initializing coordinates for nearest neighbor search
	tree=KDTree(rAu_nounit)

	# search for all points within distance r. +1 for rounding effects
	r=dnn_nounit+1
	nn=inrange(tree,rAu_nounit,r)

	# remove same atom pairs from nn list
	removeiipairs!(nn)

	# output as 528 length array of arrays
	return nn
end

############################################################################################################

"""
transforms Au 111 coordinates (in 3 x N matrix form) to coordinates w.r.t. 100. 

for use in figuring out which case of Aij (force matrix) to use.
"""
function transformAuCoord(rAu_111)
	# transformation matrix. xyz unit vectors stored in rows. negative wrt to U in fortram program (eg see GetNN2 function) but is ok. see notebook for derivation.
	u=[1/√2 -1/√6 1/√3; 0 2/√6 1/√3; -1/√2 -1/√6 1/√3]

	# transform Au coordinates to 100
	u*rAu_111
end

############################################################################################################

"""
outputs nearest neighbor distances from atom i as matrix. nn xyz coordinates stored in columns. 

called in getdrij() only

nni: nn[i] (ith entry in nearest neighbor array)
"""
function getnnimatrix(rAu,nni,i)
	# atom i coordinates
    icoords=rAu[:,i]

	# atom i nn coordinates in array
    nnicoords=[rAu[:,j] for j in nni]

	# format nn coordinates as matrix. xyz coords stored in columns
    nnimatrix=reduce(hcat,nnicoords)

	# output nn distances from atom i as matrix. xyz coords stored in columns
    nnimatrix=nnimatrix.-icoords
end

############################################################################################################

"""
outputs drij (ri-rj) for all nearest neighbors of each atom as array of matrices. 

input is matrix containing coords of Au atoms (xyz stored in columns). 

see roy art, p7, eq 20
"""
function getdrij(rAu)
    [getnnimatrix(rAu,nn[i],i) for i in eachindex(nn)]
end

############################################################################################################

"""
initializes variable mapping nearest neighbor pair distances to case number as Dict. 

see begbie p191
"""
function initmapcases()
	# possible cases of nn distances
	possiblers=[[0,1,1],[0,1,-1],[1,0,1],[-1,0,1],[1,1,0],[1,-1,0],[0,-1,-1],[0,-1,1],[-1,0,-1],[1,0,-1],[-1,-1,0],[-1,1,0]]
	
	# store cases/case numbers (1-12) in dictionary
	Dict(zip(possiblers,eachindex(possiblers)))
end

#initialize cases
map_r_to_cases=initmapcases()

############################################################################################################

"""
outputs case number for each nn of atom i as array. 
	
input is matrix with nn coords (xyz stored in columns).

called only in getcase(). 

see begbie p191 for cases
"""
function getnnicasearray(r_nni)
	# rounding distances to Int
	r_nni_rounded=round.(Int,r_nni)

	# match nn pairs to case numbers
	[map_r_to_cases[i] for i in eachcol(r_nni_rounded)]
end

############################################################################################################

"""
outputs case number for each atom's nn's as array of arrays.
	
see begbie p191 for cases
"""
function getcase()
	# transform Au coordinates to 100
	rAut=transformAuCoord(rAu)

	# find nn distances for each atom
	rijt=getdrij(rAut)

	# scale distances by a/2
	a=au.a[1]
	rijt/=a/2

	# match nn pairs to case numbers for all atoms
    [getnnicasearray(rijt[i]) for i in eachindex(rijt)]
end

############################################################################################################

"""
initializes variable mapping case number to Aij (force matrix) as Dict. 

Aij depending on case. see begbie p193 for cases
"""
function initAij()
	# assign variables for convienence
	α=au.α[1]
	β=au.β[1]
	γ=au.γ[1]

	# grabbing units for α. needed for giving 0's the right units
	aunit=unit(α)

	# possible force matrices wrt 100
	m1=[α 0*aunit 0*aunit; 0*aunit β γ; 0*aunit γ β]
	m7=[α 0*aunit 0*aunit; 0*aunit β γ; 0*aunit γ β]

	m2=[α 0*aunit 0*aunit; 0*aunit β -γ; 0*aunit -γ β]
	m8=[α 0*aunit 0*aunit; 0*aunit β -γ; 0*aunit -γ β]

	m3=[β 0*aunit γ; 0*aunit α 0*aunit; γ 0*aunit β]
	m9=[β 0*aunit γ; 0*aunit α 0*aunit; γ 0*aunit β]

	m4=[β 0*aunit -γ; 0*aunit α 0*aunit; -γ 0*aunit β]
	m10=[β 0*aunit -γ; 0*aunit α 0*aunit; -γ 0*aunit β]

	m5=[β γ 0*aunit; γ β 0*aunit; 0*aunit 0*aunit α]
	m11=[β γ 0*aunit; γ β 0*aunit; 0*aunit 0*aunit α]

	m6=[β -γ 0*aunit; -γ β 0*aunit; 0*aunit 0*aunit α]
	m12=[β -γ 0*aunit; -γ β 0*aunit; 0*aunit 0*aunit α]

	# all force matrices in array
	A_100=[m1,m2,m3,m4,m5,m6,m7,m8,m9,m10,m11,m12]

	# convert force matrices to 111 and store in array. U'*A_100*U=A_111. U is transformation matrix (see notebook for derivation). ' is transpose
	u=[1/√2 -1/√6 1/√3; 0 2/√6 1/√3; -1/√2 -1/√6 1/√3]
	A_100U=[A_100[i]*u for i in eachindex(A_100)]
	A_111=[u'*A_100U[i] for i in eachindex(A_100U)]
	
	# store force matrices/case numbers (1-12) in dictionary
	Dict(zip(1:12,A_111))
end

############################################################################################################

"""
outputs force matrix (Aij) for each atom's nn's (depending on case number).

see begbie p193 for matrices
"""
function initAijarray()
	[[Aij[j] for j in nncase[i]] for i in eachindex(nncase)]
end

############################################################################################################

"""
outputs forces on atom i due to distance away from each nn as array of arrays

called only in F_Au()
"""
function getnnifarray(rij_i,i)
	# F=−Aij(nn case no)*rij
	[-Aijarray[i][j]*rij_i[:,j] for j in axes(rij_i,2)]
end

############################################################################################################

"""
outputs forces on each Au Atom (due to nn's) as array of arrays

input is matrix with current Au coords (xyz stored in columns)

see roy art, p7, eq 20. see fortran GetFNN2
"""
function F_Au(rAu;inVfunc::Bool=false)
	# get current nn distances for each atom
	drij=getdrij(rAu)

	# find distances away from equilibrium positions. array of matrices
	rij=drij-r0ij

	# temp array with forces on each atom due to nn's. array of arrays of arrays
	farray=[getnnifarray(rij[i],i) for i in eachindex(rij)]

	if inVfunc
		# only if called in V_AuAu function. option added for speed
		return rij,farray
	else
		# normal behavior. return only forces

		# sum forces on each atom
		f=sum.(farray)
		return f
	end
end

############################################################################################################

"""
gets v for atom i. 
	
called only in V_AuAu()
"""
function getvi(rij_i,farray_i)
	# array of v's from each nn. V=r*-F. F=-A*r
	varray=[dot(rij_i[:,j],-farray_i[j]) for j in axes(rij_i,2)]

	# total v for atom i is sum of v's due to each nn
	sum(varray)
end

############################################################################################################

"""
outputs potential energy of Au slab due to nn's

input is matrix with current Au coords (xyz stored in columns)

see roy art, p7, eq 20. see fortran GetVNN2
"""
function V_AuAu(rAu;fv::Bool=false)
	# V_AuAu employs same temp variables as F_Au. rij: array of matrices. farray: array of arrays of arrays
	rij,farray=F_Au(rAu;inVfunc=true)

	# v for each atom. array
	varray=[getvi(rij[i],farray[i]) for i in eachindex(rij)]

	# total v is sum of v's for each atom. divide by 2 because of v definition
	v=sum(varray)/2

	if fv
		# return force and potential energy simulataneously. option added for speed
		f=sum.(farray)
		return f,v
	else
		# normal behavior. return only potential energy
		return v
	end
end