# potential functions for NO/Au interactions on ground, excited, and coupled states. functions based on forms in roy art.

############################################################################################################

"""
Au/O interaction potential for ground state. see roy art, p6, eq 9

dr=|ri-rO|
"""
function V00_AuO(dr::Unitful.Length)
	A0=PES_GS.A0[1]
	α0=PES_GS.α0[1]
	r_cutoff=PES_GS.AuOcutoff[1]

	eterm1=exp(-α0*dr)
	eterm2=exp(-α0*r_cutoff)

	# output V
	A0*(eterm1-eterm2) 
end

############################################################################################################

"""
Au/N interaction potential for ground state. see roy art, p6, eq 10

dr=|ri-rN|
"""
function V00_AuN(dr::Unitful.Length)
	B0=PES_GS.B0[1]
	β0=PES_GS.β0[1]
	r_cutoff=PES_GS.AuNcutoff[1]

	eterm1=exp(-β0*dr)
	eterm2=exp(-β0*r_cutoff)

	# output V
	B0*(eterm1-eterm2) 
end

############################################################################################################

"""
N/O interaction potential for ground state. see roy art, p6, eq 11

dr=|rN-rO|
"""
function V00_NO(dr::Unitful.Length)
	F0=PES_GS.F0[1]
	γ0=PES_GS.γ0[1]
	r0_NO=PES_GS.r0_NO[1]

	ddr=dr-r0_NO
	eterm=exp(-γ0*ddr)

	# output V
	F0*(1-eterm)^2 
end

############################################################################################################

"""
Au/O interaction potential for excited state. see roy art, p6, eq 14

dr=|ri-rO|
"""
function V11_AuO(dr::Unitful.Length)
	A1=PES_ionic.A1[1]
	α1=PES_ionic.α1[1]
	r_cutoff=PES_ionic.AuOcutoff[1]

	eterm1=exp(-α1*dr)
	eterm2=exp(-α1*r_cutoff)

	# output V
	A1*(eterm1-eterm2) 
end

############################################################################################################

"""
Au/N interaction potential for excited state. see roy art, p6, eq 15

dr=|ri-rN|

cos_th=(zO-zN)/|rN-rO|
"""
function V11_AuN(dr::Unitful.Length,cos_th)
	B1=PES_ionic.B1[1]
	β1=PES_ionic.β1[1]
	r1_AuN=PES_ionic.r1_AuN[1]
	r_cutoff=PES_ionic.AuNcutoff[1]

	ddr1=dr-r1_AuN
	ddr2=r_cutoff-r1_AuN

	eterm1=exp(-β1*ddr1)
	eterm2=exp(-β1*ddr2)

	a=eterm1^2-eterm2^2
	b=2*cos_th^2*eterm1
	c=eterm2

	# output V
	B1*(a-b-c) 
end

############################################################################################################

"""
N/O interaction potential for excited state. see roy art, p6, eq 16

dr=|rN-rO|
"""
function V11_NO(dr::Unitful.Length)
	F1=PES_ionic.F1[1]
	γ1=PES_ionic.γ1[1]
	r1_NO=PES_ionic.r1_NO[1]

	ddr=dr-r1_NO
	eterm=exp(-γ1*ddr)

	# output V
	F1*(1-eterm)^2 
end

############################################################################################################

"""
Au/NO interaction potential for excited state at far distances. see roy art, p6, eq 13

zcom is the perpendicular distance of the center of mass of the NO molecule from the surface plane (z=0)
"""
function V11_image(zcom::Unitful.Length)
	D=PES_ionic.D[1]
	C=PES_ionic.C[1]
	zimage=PES_ionic.zimage[1]

	dz=zcom-zimage
	expr=C^2+dz^2
	sqexpr=sqrt(expr)

	# output V
	-D/sqexpr 
end

############################################################################################################

"""
Au/O interaction potential for coupled state. see roy art, p7, eq 18. variables defined differently in fortran

dr=|ri-rO|
"""
function V01_AuO(dr::Unitful.Length)
	A2=PES_coup.A2[1]
	A3=PES_coup.A3[1]
	γ2=PES_coup.γ2[1]
	r_cutoff=PES_coup.AuOcutoff[1]

	eterm1=exp(γ2*dr)
	eterm2=exp(γ2*r_cutoff)

	expr1=1+A3*eterm1
	expr2=1+A3*eterm2

	# output V
	-A2*(expr1^-1-expr2^-1) 
end

############################################################################################################

"""
Au/N interaction potential for coupled state. see roy art, p7, eq 19. variables defined differently in fortran

dr=|ri-rN|
"""
function V01_AuN(dr::Unitful.Length)
	B2=PES_coup.B2[1]
	B3=PES_coup.B3[1]
	γ3=PES_coup.γ3[1]
	r_cutoff=PES_coup.AuNcutoff[1]

	eterm1=exp(γ3*dr)
	eterm2=exp(γ3*r_cutoff)

	expr1=1+B3*eterm1
	expr2=1+B3*eterm2

	# output V
	-B2*(expr1^-1-expr2^-1)
end

