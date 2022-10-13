# adiabatic molecular dynamics of NO scattering off of an Au(111) surface

# system of units
# mass = g/mol = amu
# distance = Angstrom = 10^(-10) m
# time = 10 femtoseconds = 10^(-14) s
# frequency = 10^(14) s^(-1) = 100 THz
# velocity = Angstrom/(10 femtoseconds) = 10^4 m/s
# energy = 100 kJ/mol

using DataFrames
using CSV
using NearestNeighbors

#functions
function initAuParams()
	file="au params.csv"
	au=CSV.read(file,DataFrame)
	au.α[1]*=6.02e-2
	au.β[1]*=6.02e-2
	au.γ[1]*=6.02e-2
	return au
end

function initSimParams()
	file="input.csv"
	CSV.read(file,DataFrame)
end

function initPESParamGS()
	file="pes params-ground.csv"
	CSV.read(file,DataFrame)
end

function initPESParamIonic()
	file="pes params-ionic.csv"
	CSV.read(file,DataFrame)
end

function initPESParamCoup()
	file="pes params-coupled.csv"
	CSV.read(file,DataFrame)
end

function initNOParams()
	file="no params.csv"
	CSV.read(file,DataFrame)
end
# PES
# neutral
"""
article, p6, eq 9

dr=|ri-rO|
"""
function V00_AuO(dr)
	A0=PES_GS.A0[1]
	α0=PES_GS.α0[1]
	r_cutoff=PES_GS.AuOcutoff[1]
	println(A0,", ",α0,", ",r_cutoff) #check

	eterm1=exp(-α0*dr)
	eterm2=exp(-α0*r_cutoff)
	println(eterm1,", ",eterm2) #check

	v=A0*(eterm1-eterm2) 
	println(v) #check
end

"""
article, p6, eq 10

dr=|ri-rN|
"""
function V00_AuN(dr)
	B0=PES_GS.B0[1]
	β0=PES_GS.β0[1]
	r_cutoff=PES_GS.AuNcutoff[1]
	println(B0,", ",β0,", ",r_cutoff) #check

	eterm1=exp(-β0*dr)
	eterm2=exp(-β0*r_cutoff)
	println(eterm1,", ",eterm2) #check

	v=B0*(eterm1-eterm2) 
	println(v) #check
end

"""
article, p6, eq 11

dr=|rN-rO|
"""
function V00_NO(dr)
	F0=PES_GS.F0[1]
	γ0=PES_GS.γ0[1]
	r0_NO=PES_GS.r0_NO[1]
	println(F0,", ",γ0,", ",r0_NO) #check

	ddr=dr-r0_NO
	eterm=exp(-γ0*ddr)
	println(ddr,", ",eterm) #check

	v=F0*(1-eterm)^2 
	println(v) #check
end

# ionic

"""
article, p6, eq 14

dr=|ri-rO|
"""
function V11_AuO(dr)
	A1=PES_ionic.A1[1]
	α1=PES_ionic.α1[1]
	r_cutoff=PES_ionic.AuOcutoff[1]
	println(A1,", ",α1,", ",r_cutoff) #check

	eterm1=exp(-α1*dr)
	eterm2=exp(-α1*r_cutoff)
	println(eterm1,", ",eterm2) #check

	v=A1*(eterm1-eterm2) 
	println(v) #check
end

"""
article, p6, eq 15

dr=|ri-rN|

cos_th=(zO-zN)/|rN-rO|
"""
function V11_AuN(dr,cos_th)
	B1=PES_ionic.B1[1]
	β1=PES_ionic.β1[1]
	r1_AuN=PES_ionic.r1_AuN[1]
	r_cutoff=PES_ionic.AuNcutoff[1]
	println(B1,", ",β1,", ",r1_AuN,", ",r_cutoff) #check

	ddr1=dr-r1_AuN
	ddr2=r_cutoff-r1_AuN
	println(ddr1,", ",ddr2) #check

	eterm1=exp(-β1*ddr1)
	eterm2=exp(-β1*ddr2)
	println(eterm1,", ",eterm2) #check

	a=eterm1^2-eterm2^2
	b=2*cos_th^2*eterm1
	c=eterm2
	println(a,", ",b,", ",c) #check

	v=B1*(a-b-c) 
	println(v) #check
end

"""
article, p6, eq 16

dr=|rN-rO|
"""
function V11_NO(dr)
	F1=PES_ionic.F1[1]
	γ1=PES_ionic.γ1[1]
	r1_NO=PES_ionic.r1_NO[1]
	println(F1,", ",γ1,", ",r1_NO) #check

	ddr=dr-r1_NO
	eterm=exp(-γ1*ddr)
	println(ddr,", ",eterm) #check

	v=F1*(1-eterm)^2 
	println(v) #check
end

"""
article, p6, eq 13

zcom is the perpendicular distance of the center of mass of the NO molecule from the surface plane (z=0)
"""
function V11_image(zcom)
	D=PES_ionic.D[1]
	C=PES_ionic.C[1]
	zimage=PES_ionic.zimage[1]
	println(D,", ",C,", ",zimage) #check

	dz=zcom-zimage
	expr=C^2+dz^2
	sqexpr=sqrt(expr)
	println(dz,", ",expr,", ",sqexpr) #check

	v=-D/sqexpr 
	println(v) #check
end

# coupling
"""
article, p7, eq 18

variables defined differently in fortran prgm

dr=|ri-rO|
"""
function V01_AuO(dr)
	A2=PES_coup.A2[1]
	A3=PES_coup.A3[1]
	γ2=PES_coup.γ2[1]
	r_cutoff=PES_coup.AuOcutoff[1]
	println(A2,", ",A3,", ",γ2,", ",r_cutoff) #check

	eterm1=exp(γ2*dr)
	eterm2=exp(γ2*r_cutoff)
	println(eterm1,", ",eterm2) #check

	expr1=1+A3*eterm1
	expr2=1+A3*eterm2
	println(expr1,", ",expr2) #check

	v=-A2*(expr1^-1-expr2^-1) 
	println(v) #check
end

"""
article, p7, eq 19

variables defined differently in fortran prgm

dr=|ri-rN|
"""
function V01_AuN(dr)
	B2=PES_coup.B2[1]
	B3=PES_coup.B3[1]
	γ3=PES_coup.γ3[1]
	r_cutoff=PES_coup.AuNcutoff[1]
	println(B2,", ",B3,", ",γ3,", ",r_cutoff) #check

	eterm1=exp(γ3*dr)
	eterm2=exp(γ3*r_cutoff)
	println(eterm1,", ",eterm2) #check

	expr1=1+B3*eterm1
	expr2=1+B3*eterm2
	println(expr1,", ",expr2) #check

	v=-B2*(expr1^-1-expr2^-1) 
	println(v) #check
end
# left off \below
"""
input is vector of vectors. removes same atom pairs from nearest neighbor list. called in getnn() only.
"""	
function removeiipairs!(nn)
	for i in 1:length(nn)
		for j in 1:length(nn[i])
			if nn[i][j]==i
				splice!(nn[i],j)
				break
			end
		end
	end
end
# """
# throws error if idx not found
# """
# function removeiipairs2!(nn)
# 	for i in 1:length(nn)
# 		idx=findfirst(isequal(i),nn[i])
# 		deleteat!(nn[i],idx)
# 	end
# end
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
finds drij (ri-rj) for all nearest neighbor pairs. input is matrix containing coords of Au atoms

article, p7, eq 20
"""
function getdrij(rAu)
	drij=[]
	for i in 1:length(nn)
		ithrow=[]
		for j in 1:length(nn[i])
			nnj=nn[i][j]
			dx=rAu[1,nnj]-rAu[1,i]
			dy=rAu[2,nnj]-rAu[2,i]
			dz=rAu[3,nnj]-rAu[3,i]
			push!(ithrow,[dx,dy,dz])
		end
		push!(drij,ithrow)
	end
	return drij
end
"""
finds case number for each nn pair. see begbie pt 2 p 191 for cases
"""
function getcase()
	rAut=transformAuCoord(rAu)
	rijt=getdrij(rAut);
	a=au.a[1]
	rijt/=a/2
	possiblers=[[0,1,1],[0,1,-1],[1,0,1],[-1,0,1],[1,1,0],[1,-1,0],[0,-1,-1],[0,-1,1],[-1,0,-1],[1,0,-1],[-1,-1,0],[-1,1,0]]
	mapcases=Dict(zip(possiblers,1:12))
	cases=[]
	for i in 1:length(rijt)
		casesithrow=[]
		for j in 1:length(rijt[i])
			caseno=mapcases[round.(Int,rijt[i][j])]
			push!(casesithrow,caseno)
		end
		push!(cases,casesithrow)
	end
	return cases
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
	fmat100u=[fmat100[i]*u for i in 1:length(fmat100)]
	fmat111=[u'*fmat100u[i] for i in 1:length(fmat100u)]
	# fmat111=u'.*fmat100.*u
	Aij=Dict(zip(1:12,fmat111))
end
# """
# article, p7, eq 20

# m is a matrix containing the current coords of all the atoms in the Au slab. 
# """
# function V_AuAu(m)
# 	drij=getdrij(m)
# 	rij=drij-r0ij

# 	v=0
# 	for i in 1:length(rij)
# 		for j in 1:length(rij[i])
# 			case=nncase[i][j]
# 			looprij=rij[i][j]
# 			looprijt=rij[i][j]'
# 			loopA=Aij[case]
# 			v+=looprijt*loopA*looprij
# 		end
# 	end
# 	v/=2
# end
"""
article, p7, eq 20

m is a matrix containing the current coords of all the atoms in the Au slab. 
"""
function F_Au(m::Matrix{Float64})
	drij=getdrij(m)
	rij=drij-r0ij

	f=[-Aij[nncase[i][j]]*rij[i][j] for i in 1:length(rij) for j in 1:length(rij[i])]
end
function F_Au(rij::Vector{Any})
	f=[-Aij[nncase[i][j]]*rij[i][j] for i in 1:length(rij) for j in 1:length(rij[i])]
end
"""
article, p7, eq 20

m is a matrix containing the current coords of all the atoms in the Au slab. 
"""
function V_AuAu(m::Matrix{Float64})
	drij=getdrij(m)
	rij=drij-r0ij

	rarray=[rij[i][j]' for i in 1:length(rij) for j in 1:length(rij[i])]
	v=-rarray⋅F_Au(rij)/2
	println(vij[1])
end
function V_AuAu(rij::Vector{Any})
	rarray=[rij[i][j]' for i in 1:length(rij) for j in 1:length(rij[i])]
	v=-rarray⋅F_Au(rij)/2
	println(vij[1])
end

#initialize variables from input files
au=initAuParams()
param=initSimParams()
PES_GS=initPESParamGS()
PES_ionic=initPESParamIonic()
PES_coup=initPESParamCoup()
no=initNOParams()

#other variables
mAu=196.966569
mN=14.00674
mO=15.9994

Av_Et = 0
Av_theta = 0
Av_weighted_theta = 0
Norm_weighted_theta = 0
Av_Evib = 0
Av_Er = 0
Av_EAu = 0
NTRAP = 0

rAu=transpose([au.x au.y au.z]) # 3 x N matrix. xyz coords of each Au atom stored in columns
nn=getnn() # nearest neighbor array for each Au atom. ith row corresponds to Au atom i's nearest neighbors (in terms of atom number)
r0ij=getdrij(rAu) # array of vectors to each Au atom's nearest neighbors. ith row cooresponds to all vectors (in [x,y,z] format) between Au atom i and its nearest neighbors
nncase=getcase() # case # array for each Au atom. ith row cooresponds to the case #'s between Au atom i and its nearest neighbors 
Aij=initAij() # force matrix as a function of case #

run(`pwd`);
# include("test.jl")
# # gen test rAu
# testmat=copy(rAu)
# for i in 1:size(testmat,1)
# 	for j in 1:size(testmat,2)
# 		testmat[i][j]+=rand()
# 	end
# end


#allocatables
#=
nn=
gauss_num=
xi=
vi=
xiv=
F=
F_next=
F_temp=
xiN=
xiO=
rvNO=
F_neutral=
F_ion=
F_coup=
=#

#=
##########################################################################
#old


#constants
natoms=2
atom1="O"
atom2="O"

c=2.99792458e8
re=1.2075273
wn=1580.193 #cm-1
freq=wn*100*c/10^14
m1=16
m2=16
mred=m1*m2/(m1+m2)
k=4π^2*mred*freq^2

#parameters
et=1.5
dt=0.005
nstep=1000

#functions
function outputenergies()
	data=DataFrame(t=storedt,PE=storedpe,KE=storedke,Etot=storedpe+storedke)
	file="energy.csv"
	CSV.write(file,data)
end

function outputvibration()
	data=DataFrame(t=storedt,r=storedr)
	file="vibration.csv"
	CSV.write(file,data)
end

function outputtrajectory()
	file=open("trajectory.xyz","w")
	for i∈1:length(storedx1)
		println(file,natoms)
		println(file)
		curx1=storedx1[i]
		curx2=storedx2[i]
		println(file,"$atom1 $curx1 $y1 $z1")
		println(file,"$atom2 $curx2 $y2 $z2")
	end
	close(file)
end

function outputall()
	outputenergies()
	outputvibration()
	outputtrajectory()
end

function stepall()
	#step x
	global x1=x1+v1*dt+1/2*a1*dt^2
	global x2=x2+v2*dt+1/2*a2*dt^2

	#step a
	olda1=a1
	olda2=a2
	r=x2-x1
	x=r-re
	f=-k*x
	global a1=-f/m1
	global a2=f/m2

	#step v
	global v1=v1+(a1+olda1)/2*dt
	global v2=v2+(a2+olda2)/2*dt

	#step e
	global pe=1/2*k*x^2
	global ke=1/2*mred*(v2-v1)^2

	return x1,x2,v1,v2,a1,a2,pe,ke

end
#=
function stepx()
	global x1=x1+v1*dt+1/2*a1*dt^2
	global x2=x2+v2*dt+1/2*a2*dt^2
	return x1,x2
end
#=
function stepx2()
	global x2=x2+v2*dt+1/2*a2*dt^2
end
=#
function stepv()
	alla=stepa()
	nexta1=first(alla)
	global v1=v1+(a1+nexta1)/2*dt

	nexta2=last(alla)
	global v2=v2+(a2+nexta2)/2*dt
	return v1,v2
end
#=
function stepv2()
	nexta2=stepa2()
	global v2=v2+(a2+nexta2)/2*dt
end
=#
function stepa()
	nextx1,nextx2=stepx()
	r=nextx2-nextx1
	x=r-re
	f=-k*x
	global a1=-f/m1
	global a2=f/m2
	return a1,a2
end

=#
#=
function stepa2()
	tmpx1=stepx1()
	tmpx2=stepx2()
	r=tmpx2-tmpx1
	x=r-re
	f=-k*x
	global a2=f/m2
end


function curE()
	r=x2-x1
	x=r-re

	global pe=1/2*k*x^2
	global ke=et-pe
	#global te=pe+ke
	return pe,ke
end
=#
function storeall()
	push!(storedt,t)
	push!(storedx1,x1)
	push!(storedv1,v1)
	push!(storeda1,a1)
	push!(storedx2,x2)
	push!(storedv2,v2)
	push!(storeda2,a2)
	push!(storedr,x2-x1)
	push!(storedke,ke)
	push!(storedpe,pe)
	#push!(storedte,te)
end

#initial conditions
x1=-m2*re/(m1+m2)
x2=m1*re/(m1+m2)

y1=0
z1=0
y2=0
z2=0

r=x2-x1
x=r-re
pe=1/2*k*x^2
ke=et-pe
#te=0
#curE()

vrel=√(2ke/mred)
v1=m2*vrel/(m1+m2)
v2=-m1*vrel/(m1+m2)

f=-k*x
a1=-f/m1
a2=f/m2

storedx1=[]
storedv1=[]
storeda1=[]

storedx2=[]
storedv2=[]
storeda2=[]

storedr=[]

storedpe=[]
storedke=[]
#storedte=[]

storedt=[]

t=0

#main

for i∈0:nstep
	storeall()
	stepall()
	global t+=dt
end
outputall()
#=
out=open("test.txt","w")

for i∈1:nstep
	println(out,storedx1[i])
end
=#
=#