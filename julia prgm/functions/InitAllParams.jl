using DataFrames
using CSV
import PhysicalConstants.CODATA2018: N_A, k_B

function initAuParams()
	file="input files/au params.csv"
	au=CSV.read(file,DataFrame)

	#convert x,y,z,dnn,a,aPBCx,aPBCy,aPBCz to d_MD
	au.x*=u"Å"
	au.y*=u"Å"
	au.z*=u"Å"
	au.dnn*=u"Å"
	au.a*=u"Å"
	au.aPBCx*=u"Å"
	au.aPBCy*=u"Å"
	au.aPBCz*=u"Å"

	#convert alpha,beta,gamma to u"N/m" and multiply by N_A
	au.α*=u"N/m"*N_A
	au.β*=u"N/m"*N_A
	au.γ*=u"N/m"*N_A

	return au
end

function initSimParams()
	file="input files/input.csv"
	input=CSV.read(file,DataFrame)

	#convert T,dt to u"K",u"fs"
	input.T*=u"K"
	input.dt*=u"fs"

	return input
end

function initPESParamGS()
	file="input files/pes params-ground.csv"
	gs=CSV.read(file,DataFrame)

	#convert A0	α0	AuOcutoff	B0	β0	AuNcutoff	F0	γ0	r0_NO
	gs.A0*=u"kJ/mol"
	gs.α0*=u"Å"
	gs.AuOcutoff*=u"Å"
	gs.B0*=u"kJ/mol"
	gs.β0*=u"Å"
	gs.AuNcutoff*=u"Å"
	gs.F0*=u"kJ/mol"
	gs.γ0*=u"Å"
	gs.r0_NO*=u"Å"

	return gs
end

function initPESParamIonic()
	file="input files/pes params-ionic.csv"
	ion=CSV.read(file,DataFrame)

	#convert A1	α1	AuOcutoff	B1	β1	AuNcutoff	r1_AuN	D	C	zimage	F1	γ1	r1_NO	ϕ	Ea
	ion.A1*=u"kJ/mol"
	ion.α1*=u"Å"
	ion.AuOcutoff*=u"Å"
	ion.B1*=u"kJ/mol"
	ion.β1*=u"Å"
	ion.AuNcutoff*=u"Å"
	ion.r1_AuN*=u"Å"
	ion.D*=u"kJ/mol"
	ion.C*=u"Å"
	ion.zimage*=u"Å"
	ion.F1*=u"kJ/mol"
	ion.γ1*=u"Å"
	ion.r1_NO*=u"Å"
	ion.ϕ*=u"eV"
	ion.Ea*=u"eV"

	return ion
end

function initPESParamCoup()
	file="input files/pes params-coupled.csv"
	coup=CSV.read(file,DataFrame)

	#convert A2	A3	γ2	AuOcutoff	B2	B3	γ3	AuNcutoff
	coup.A2*=u"kJ/mol"
	coup.A3*=u"kJ/mol"
	coup.γ2*=u"Å"
	coup.AuOcutoff*=u"Å"
	coup.B2*=u"kJ/mol"
	coup.B3*=u"kJ/mol"
	coup.γ3*=u"Å"
	coup.AuNcutoff*=u"Å"

	return coup
end

function initNOParams()
	file="input files/no params.csv"
	no=CSV.read(file,DataFrame)

	#convert r,v,Et_i,θi to d_MD,v_MD,e_MD,u"°"
	no.r*=u"Å"
	no.v*=u"v_MD"
	no.Et_i*=u"e_MD"
	no.θi*=u"°"

	return no
end

function initAuCoords()
	t=transpose([au.x au.y au.z])
	# convert(Matrix{typeof(t)},t)
	# convert(Array{typeof(t),2},t)
end