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