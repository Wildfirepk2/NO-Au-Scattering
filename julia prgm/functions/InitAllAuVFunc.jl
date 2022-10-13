using NearestNeighbors
using LinearAlgebra

"""
input is vector of vectors. removes same atom pairs from nearest neighbor list. called in getnn() only.
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

"""
identifies 1st nearest neighbors for set of Au atoms. Au's crystal structure is fcc. thus, each atom (except those at the edges) has 12 nearest neighbors.
"""
function getnn()
	tree=KDTree(rAu)
	# println(tree)
	r=au.dnn[1]+1 # nearest neighbor dist for fcc. +1 for rounding effects
	# println(r)
	# i,d=knn(tree,pts,13)
	# println(i[123],"\n",d[123])
	nn=inrange(tree,rAu,r)
	# println(nn[123])
	removeiipairs!(nn)
	# println(nn[123])
	return nn
end

"""
transforms Au 111 to 100. for use in figuring out which case of Aij to use. see notebook for u derivation.
"""
function transformAuCoord(m)
	u=[1/√2 -1/√6 1/√3; 0 2/√6 1/√3; -1/√2 -1/√6 1/√3] # transformation matrix. xyz unit vectors stored in rows
	rt=u*m
end

"""
creates a matrix of atom i's nearest neighbor coordinates. xyz stored in rows. called only in getdrij()

nni: nn[i]
"""
function getnnimatrix(rAu,nni,i)
    icoords=rAu[:,i]
    nnicoords=[rAu[:,j] for j in nni]
    nnimatrix=reduce(hcat,nnicoords)
    nnimatrix=nnimatrix.-icoords
end

"""
finds drij (ri-rj) for all nearest neighbor pairs. input is matrix containing coords of Au atoms

article, p7, eq 20
"""
function getdrij(rAu)
    [getnnimatrix(rAu,nn[i],i) for i in eachindex(nn)]
end

"""
possible cases for nearest neighbor pairs

see begbie pt 2 p 191
"""
function initmapcases()
	possiblers=[[0,1,1],[0,1,-1],[1,0,1],[-1,0,1],[1,1,0],[1,-1,0],[0,-1,-1],[0,-1,1],[-1,0,-1],[1,0,-1],[-1,-1,0],[-1,1,0]]
	mapcases=Dict(zip(possiblers,eachindex(possiblers)))
end
mapcases=initmapcases()

"""
gets cases for each col in ith matrix of rij. output is array. called only in getcase()

see begbie pt 2 p 191 for cases
"""
function getnnicasearray(mi)
	mirounded=round.(Int,mi)
	nnicasearray=[mapcases[i] for i in eachcol(mirounded)]
end

"""
finds case number for each nn pair. see begbie pt 2 p 191 for cases
"""
function getcase()
	rAut=transformAuCoord(rAu)
	rijt=getdrij(rAut)
	a=au.a[1]
	rijt/=a/2

    cases=[getnnicasearray(rijt[i]) for i in eachindex(rijt)]
end

"""
Aij depending on case. see begbie pt 2 p 193 for cases
"""
function initAij()
	α=au.α[1]
	β=au.β[1]
	γ=au.γ[1]

	m1=[α 0 0; 0 β γ; 0 γ β]
	m7=[α 0 0; 0 β γ; 0 γ β]

	m2=[α 0 0; 0 β -γ; 0 -γ β]
	m8=[α 0 0; 0 β -γ; 0 -γ β]

	m3=[β 0 γ; 0 α 0; γ 0 β]
	m9=[β 0 γ; 0 α 0; γ 0 β]

	m4=[β 0 -γ; 0 α 0; -γ 0 β]
	m10=[β 0 -γ; 0 α 0; -γ 0 β]

	m5=[β γ 0; γ β 0; 0 0 α]
	m11=[β γ 0; γ β 0; 0 0 α]

	m6=[β -γ 0; -γ β 0; 0 0 α]
	m12=[β -γ 0; -γ β 0; 0 0 α]

	fmat100=[m1,m2,m3,m4,m5,m6,m7,m8,m9,m10,m11,m12]
	u=[1/√2 -1/√6 1/√3; 0 2/√6 1/√3; -1/√2 -1/√6 1/√3] # transformation matrix. xyz unit vectors stored in rows
	fmat100u=[fmat100[i]*u for i in eachindex(fmat100)]
	fmat111=[u'*fmat100u[i] for i in eachindex(fmat100u)]
	# fmat111=u'.*fmat100.*u
	Aij=Dict(zip(1:12,fmat111))
end

"""
finds force matrix (Aij) for each nn pair case. see begbie pt 2 p 193 for matrices
"""
function initAijarray()
	aijarray=[[Aij[j] for j in nncase[i]] for i in eachindex(nncase)]
end

"""
gets force array for ith rij. called only in F_Au()
"""
function getnnifarray(mi,i)
	[-Aijarray[i][j]*mi[:,j] for j in axes(mi,2)]
end

"""
article, p7, eq 20

m is a matrix containing the current coords of all the atoms in the Au slab. 

# return f only
"""
function F_Au(m::Matrix{Float64};inVfunc::Bool=false)
	drij=getdrij(m)
	rij=drij-r0ij

	farray=[getnnifarray(rij[i],i) for i in eachindex(rij)]

	if inVfunc
		return rij,farray
	else
		f=sum.(farray)
		return f
	end
end

"""
gets v for atom i. called only in V_Au()
"""
function getnniv(mi,ai)
	varray=[dot(mi[:,j],ai[j]) for j in axes(mi,2)]
	sum(varray)
end

"""
article, p7, eq 20

m is a matrix containing the current coords of all the atoms in the Au slab. 

# return v only
"""
function V_AuAu(m::Matrix{Float64};fv::Bool=false)
	rij,farray=F_Au(m;inVfunc=true)

	varray=[getnniv(rij[i],farray[i]) for i in eachindex(rij)]
	v=-sum(varray)/2 # negative for sign flip. farray vals negative. divide 2 because of v definition

	if fv
		f=sum.(farray)
		return f,v
	else
		return v
	end
end