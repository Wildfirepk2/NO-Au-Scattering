# force functions for NO/Au interactions on ground, excited, and coupled states. functions based on forms in roy art.

############################################################################################################

"""
Au/O force function for ground state. see roy art, p6, eq 9

dr=|ri-rO|
"""
function F00_AuO(dr::Unitful.Length)
	A0=PES_GS.A0[1]
	α0=PES_GS.α0[1]

	eterm1=exp(-α0*dr)

	# output F
	A0*α0*eterm1
end

############################################################################################################

"""
Au/N force function for ground state. see roy art, p6, eq 10

dr=|ri-rN|
"""
function F00_AuN(dr::Unitful.Length)
	B0=PES_GS.B0[1]
	β0=PES_GS.β0[1]

	eterm1=exp(-β0*dr)

	# output F
	B0*β0*eterm1
end

############################################################################################################

"""
N/O force function for ground state. see roy art, p6, eq 11

dr=|rN-rO|
"""
function F00_NO(dr::Unitful.Length)
	F0=PES_GS.F0[1]
	γ0=PES_GS.γ0[1]
	r0_NO=PES_GS.r0_NO[1]

	ddr=dr-r0_NO
	eterm=exp(-γ0*ddr)

	# output F
	-2F0*γ0*eterm*(1-eterm) 
end

############################################################################################################

"""
Au/O force function for excited state. see roy art, p6, eq 14

dr=|ri-rO|
"""
function F11_AuO(dr::Unitful.Length)
	A1=PES_ionic.A1[1]
	α1=PES_ionic.α1[1]

	eterm1=exp(-α1*dr)

	# output F
	A1*α1*eterm1
end

############################################################################################################

"""
Au/N force function for excited state. see roy art, p6, eq 15

dr=|ri-rN|

cos_th=(zO-zN)/|rN-rO|
"""
function F11_AuN(dr::Unitful.Length,cos_th)
	B1=PES_ionic.B1[1]
	β1=PES_ionic.β1[1]
	r1_AuN=PES_ionic.r1_AuN[1]
	
	ddr1=dr-r1_AuN

	eterm1=exp(-β1*ddr1)

	a=eterm1^2
	b=cos_th^2*eterm1

	# output F
	2B1*β1*(a-b) 
end

############################################################################################################

"""
output force (with direction) on N/O due to F11_AuN. accounts for cutoff? see FORCE_Au_N_ION in fortran

dr=|ri-rN|

cos_th=(zO-zN)/|rN-rO|
"""
function F11_AuNcutoff(dr::Unitful.Length,cos_th,rNO,uON)
	B1=PES_ionic.B1[1]
	β1=PES_ionic.β1[1]
	r1_AuN=PES_ionic.r1_AuN[1]
	r_cutoff=PES_ionic.AuNcutoff[1]

	ddr1=dr-r1_AuN
	ddr2=r_cutoff-r1_AuN

	eterm1=exp(-β1*ddr1)
	eterm2=exp(-β1*ddr2)

	fcut1=4B1*eterm1*(uON*cos_th^2/rNO) # vec
	fzcut1=4B1*eterm1*(cos_th/rNO) # scalar
	fcut2=4B1*eterm2*(uON*cos_th^2/rNO) # vec. fcut1 in fortran
	fzcut2=4B1*eterm2*(cos_th/rNO) # scalar. fcut2 in fortran

	f=fcut1-fcut2
	fz=fzcut1-fzcut2
	
	ft=f+[0u"N/mol",0u"N/mol",fz]
	push!(FNO_AuN,ft)
	return ft
end

############################################################################################################

"""
output force (with direction) on N/O due to F11_AuN. accounts for cutoff? see FORCE_Au_N_ION in fortran

\fix
O forces
"""
function F11_AuNcutoff()
	# println("FNO_AuN: $FNO_AuN")
	ft=-sum(FNO_AuN)
	# println("ft: $ft")
	empty!(FNO_AuN)
	# println("FNO_AuN: $FNO_AuN")
	return ft
end

############################################################################################################

"""
N/O force function for excited state. see roy art, p6, eq 16

dr=|rN-rO|
"""
function F11_NO(dr::Unitful.Length)
	F1=PES_ionic.F1[1]
	γ1=PES_ionic.γ1[1]
	r1_NO=PES_ionic.r1_NO[1]

	ddr=dr-r1_NO
	eterm=exp(-γ1*ddr)

	# output F
	-2F1*γ1*eterm*(1-eterm) 
end

############################################################################################################

"""
Au/NO force function for excited state at far distances. see roy art, p6, eq 13

zcom is the perpendicular distance of the center of mass of the NO molecule from the surface plane (z=0)
"""
function F11_image(zcom::Unitful.Length)
	D=PES_ionic.D[1]
	C=PES_ionic.C[1]
	zimage=PES_ionic.zimage[1]

	dz=zcom-zimage
	expr=C^2+dz^2
	sqexpr=expr^(3/2)

	# output F
	-D*dz/sqexpr 
end

############################################################################################################

"""
Au/O force function for coupled state. see roy art, p7, eq 18. variables defined differently in fortran

dr=|ri-rO|
"""
function F01_AuO(dr::Unitful.Length)
	A2=PES_coup.A2[1]
	A3=PES_coup.A3[1]
	γ2=PES_coup.γ2[1]

	eterm1=exp(γ2*dr)

	expr1=1+A3*eterm1

	# output F
	-A2*A3*γ2*eterm1*(expr1^-2) 
end

############################################################################################################

"""
Au/N force function for coupled state. see roy art, p7, eq 19. variables defined differently in fortran

dr=|ri-rN|
"""
function F01_AuN(dr::Unitful.Length)
	B2=PES_coup.B2[1]
	B3=PES_coup.B3[1]
	γ3=PES_coup.γ3[1]

	eterm1=exp(γ3*dr)

	expr1=1+B3*eterm1

	# output F
	-B2*B3*γ3*eterm1*(expr1^-2)
end

