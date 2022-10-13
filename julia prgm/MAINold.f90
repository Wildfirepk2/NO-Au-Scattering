!general notes
!go back and check if nums exact in csv files
!PES:O->2, N->3


PROGRAM NO_Au_111_2X2_PES_MD

IMPLICIT NONE

!INTEGER     :: N            ! Number of Au atoms
INTEGER, ALLOCATABLE :: nn(:,:)      ! List of nearest-neighbour Au atoms of each Au atom in a (111) lattice, (N x 12)
!INTEGER     :: NTRAJ        ! Number of trajectories
!INTEGER     :: NSTEPS_EQ,NSTEPS_DYN  ! Number of time steps in equilibration and dynamics loops
INTEGER     :: i,j            ! counter 
!INTEGER     :: NV           ! initial vibrational state of NO
!INTEGER     :: NVIB         ! Number of entries of initial bond length and vibrational velocity of NO(v) 
INTEGER     :: NBOUNCE      ! Number of bounces 
!INTEGER     :: NTRAP        ! Number of trajectories that trapped
!REAL(8)     :: dnn          ! Nearest-neighbour distance in Au(111)
!REAl(8)     :: a            ! Lattice constant
REAL(8)     :: U(3,3)       ! Coordinate transformation matrix r --> r'
REAL(8)     :: r0(3,12)     ! Equilibrium positions of 12 nearest neighbours
!REAL(8)     :: aPBC(3)      ! Simulation box dimensions
REAL(8), ALLOCATABLE :: gauss_num(:)     ! Gaussian random number of unit width
!REAL(8)     :: dt           ! Time step (fs)
REAL(8)     :: rNO          ! Bond length of NO
REAL(8)     :: xcom,ycom,zcom,vcom,vcomsq        ! x, y and z positions of centre-of-mass of NO, 
                                                 ! velocity and velocity-squared of centre-of-mass of NO
REAL(8)     :: Et_i,theta_i,Et,Evib,Er,vrel,theta,polar,azimuthal ! translational energy, relative velocity, incidence angle, orientation angles
REAL(8)     :: v1,v2   ! z components of velocity of centre-of-mass of NO, to determine the number of bounces
REAL(8), ALLOCATABLE :: xi(:,:),vi(:,:),x(:,:),v(:,:),F(:,:),F_next(:,:),F_temp(:,:) ! Atom positions, velocities, forces (3 x N+2)
REAL(8), ALLOCATABLE :: xiN(:,:),xiO(:,:)
!REAL(8), ALLOCATABLE :: mass(:)      ! Masses of nuclei
!REAL(8), ALLOCATABLE :: rvNO(:,:)  ! Initial bond length and vibrational velocity of NO(v)
REAL(8), ALLOCATABLE :: F_NEUTRAL(:,:),F_ION(:,:),F_COUPLING(:,:) !"Neutral", "Ion" and "Coupling" forces (3 X N+2)
REAL(8)     :: EigenV(2) ! Ground-state eigen vector
REAL(8)     :: KE,PE,!T  ! Kinetic energy, potential energy, average values, surface temperature
REAL(8)     :: V_Au_Au,E_Au_i   ! Potential energy, initial total energy of Au(111)
REAL(8)     :: RAN  ! Random number between 0 and 1
REAL(8)     :: CUTOFF
REAL(8)     :: VIB_POT,VIB_KIN  ! Vibrational potential and kinetic energies of scattered NO
!REAL(8)     :: Av_Et,Av_theta,Av_weighted_theta,Norm_weighted_theta,Av_Evib,Av_Er,Av_EAu  ! Average values of translational, vibrational and rotational energies of NO
!REAL(8), PARAMETER :: pi  = 3.14159265358979d0    
!REAL(8), PARAMETER :: RT2 = 1.41421356237310d0    ! SQRT(2)
!REAL(8), PARAMETER :: RT3 = 1.73205080756888d0    ! SQRT(3)
!REAL(8), PARAMETER :: RT6 = 2.44948974278318d0    ! SQRT(6)
REAL(8), PARAMETER :: Ang = 1.0d-10       ! 1 Angstrom in m
REAL(8), PARAMETER :: eV = 1.60217646d-19    ! 1 eV in J
REAL(8), PARAMETER :: kJpermol = 1.660538863d-21 ! 1 kJ/mol in J
REAL(8), PARAMETER :: amu = 1.660538863d-27    ! 1 amu in kg
REAL(8), PARAMETER :: kB  = 1.3806503d-23     ! Boltzmann constant in J/K
REAL(8), PARAMETER :: kB_MD = 8.31447215d-5       ! Boltzmann constant in MD units: 100kJ/K/mol

! DEFINE PARAMETERS FOR "NEUTRAL" SURFACE (Units: energy = kJ/mol, distance = Angstroms)

!REAL(8)    :: A_n, alpha_n, Au_O_cutoff_n
!REAL(8)    :: B_n, beta_n, Au_N_cutoff_n
!REAL(8)    :: F_n, delta_n, rNO_e_n

! DEFINE PARAMETERS FOR "ION" SURFACE (Units: energy = kJ/mol, distance = Angstroms)
! Work function and electron affinity read in eV

!REAL(8)    :: C_i, D_i, z0_i
!REAL(8)    :: A_i, alpha_i, Au_O_cutoff_i
!REAL(8)    :: B_i, beta_i, rN_sur_e_i, Au_N_cutoff_i
!REAL(8)    :: F_i, delta_i, rNO_e_i
!REAL(8)    :: WF, EA

! DEFINE PARAMETERS FOR "COUPLING" FUNCTION (Units: energy = kJ/mol, distance = Angstroms)

!REAL(8)    :: coup_a_N, coup_b_N, coup_beta_N, Au_N_coupling_cutoff
!REAL(8)    :: coup_a_O, coup_b_O, coup_beta_O, Au_O_coupling_cutoff
!start**********************************************************************************
! INPUT AND OUTPUT FILES

!OPEN (UNIT = 10, FILE = 'Au111_coordinates.dat', STATUS = 'OLD', ACTION = 'READ')!done
!OPEN (UNIT = 20, FILE = 'simulation_input.dat', STATUS = 'OLD', ACTION = 'READ')!done
!OPEN (UNIT = 30, FILE = 'PES_PARAMETERS.dat', STATUS = 'OLD', ACTION = 'READ')
!OPEN (UNIT = 40, FILE = 'NO_incidence_E.dat', STATUS = 'OLD', ACTION = 'READ')
!OPEN (UNIT = 50, FILE = 'NO_vib_r_v.dat', STATUS = 'OLD', ACTION = 'READ')
OPEN (UNIT = 60, FILE = 'NO_init_coord.dat', STATUS = 'UNKNOWN', ACTION = 'WRITE')
OPEN (UNIT = 70, FILE = 'NO_init_vel.dat', STATUS = 'UNKNOWN', ACTION = 'WRITE') 
OPEN (UNIT = 80, FILE = 'NO_scatter_coord.dat', STATUS = 'UNKNOWN', ACTION = 'WRITE')
OPEN (UNIT = 90, FILE = 'NO_scatter_vel.dat', STATUS = 'UNKNOWN', ACTION = 'WRITE')
OPEN (UNIT = 100, FILE = 'NO_Et_Evib_Er.dat', STATUS = 'UNKNOWN', ACTION = 'WRITE')
OPEN (UNIT = 110, FILE = 'NO_t_nbounce.dat', STATUS = 'UNKNOWN', ACTION = 'WRITE')
OPEN (UNIT = 120, FILE = 'synopsis.dat', STATUS = 'UNKNOWN', ACTION = 'WRITE')

! READ INITIAL CONDITIONS OF SIMULATION !done, calcs below coded into data file

!READ(20,*)NTRAJ
!READ(20,*)T
!READ(20,*)dt
!dt = dt*0.1d0  ! Conversion to molecular dynamics units, 1E-14 s
!READ(20,*)NSTEPS_EQ
!READ(20,*)NSTEPS_DYN

! READ DETAILS OF Au SLAB !done, calcs below coded into data file

!READ (10,*) N
!READ (10,*) dnn
!a = RT2*dnn
!READ (10,*) aPBC(1)
!READ (10,*) aPBC(2)
!aPBC(3) = RT6*dnn*100.0d0 ! Set box z-dimension to 100x lattice dimensions

! READ PARAMETERS OF NO--Au(111) GROUND-STATE POTENTIAL ENERGY SURFACE!done, calcs below coded into data file

!READ(30,*)A_n, alpha_n, Au_O_cutoff_n
!READ(30,*)B_n, beta_n, Au_N_cutoff_n
!READ(30,*)F_n, delta_n, rNO_e_n

! Conversion to MD units (Energy = 100 kJ/mol, mass = g/mol, distance = Angstroms, time = 1d-14 s)
!A_n = A_n*0.01d0 ; B_n = B_n*0.01d0 ; F_n = F_n*0.01d0

!done, calcs below coded into data file

!READ(30,*) C_i, D_i, z0_i
!READ(30,*) A_i, alpha_i, Au_O_cutoff_i
!READ(30,*) B_i, beta_i, rN_sur_e_i, Au_N_cutoff_i
!READ(30,*) F_i, delta_i, rNO_e_i
!READ(30,*) WF, EA

! Conversion to MD units (Energy = 100 kJ/mol, mass = g/mol, distance = Angstroms, time = 1E-14 s)
!D_i = D_i*0.01d0 ; A_i = A_i*0.01d0 ; B_i = B_i*0.01d0 ; F_i = F_i*0.01d0
!WF = WF*(eV/kJpermol)*0.01d0 ; EA = EA*(eV/kJpermol)*0.01d0

! DEFINE PARAMETERS FOR "COUPLING" FUNCTION (Units: energy = kJ/mol, distance = Angstroms)

!READ(30,*) coup_a_N, coup_b_N, coup_beta_N, Au_N_coupling_cutoff
!READ(30,*) coup_a_O, coup_b_O, coup_beta_O, Au_O_coupling_cutoff

! Conversion to MD units (Energy = 100 kJ/mol, mass = g/mol, distance = Angstroms, time = 1E-14 s)
!coup_a_N = coup_a_N*0.01d0 ; coup_a_O = coup_a_O*0.01d0

! ALLOCATE DYNAMIC VARIABLES

ALLOCATE(gauss_num(3*N),xi(3,N+2),x(3,N+2),v(3,N+2),F(3,N+2),F_next(3,N+2),F_temp(3,N+2),mass(N+2),nn(3:N+2,12),xiN(4,3:N+2),xiO(4,3:N+2))
ALLOCATE(F_NEUTRAL(3,N+2),F_ION(3,N+2),F_COUPLING(3,N+2))

! READ Au ATOM POSITIONS

READ (10,*) xi(:,3:N+2)
xi(:,1:2) = 0.0d0  ! Initial coordinates of N and O will be input later
x = xi

! MASSES OF N, O AND Au ATOMS

!mass(1) = 14.00674d0
!mass(2) = 15.9994d0
!mass(3:N+2) = 196.966569d0

! CHANGE r0

U = RESHAPE ( (/ -1d0/RT2,0d0,1d0/RT2, 1d0/RT6,-2d0/RT6,1d0/RT6, -1d0/RT3,-1d0/RT3,-1d0/RT3 /), (/3,3/) )
r0 = MATMUL(TRANSPOSE(U),RESHAPE ( (/0,1,1,0,1,-1,1,0,1,-1,0,1,1,1,0,1,-1,0,0, &
-1,-1,0,-1,1,-1,0,-1,1,0,-1,-1,-1,0,-1,1,0/),(/3,12/) )/2.0d0*a)

! CREATE OF NEAREST NEIGHBOUR LIST

dnn = dnn * 1.01d0
CALL GetNN2(nn)

!---------------------------------------------------------------------------------------------------
! INITIALIZING AVERAGE ENERGIES AND NTRAP, READING VALUES
Av_Et = 0.0d0
Av_theta = 0.0d0
Av_weighted_theta = 0.0d0
Norm_weighted_theta = 0.0d0
Av_Evib = 0.0d0
Av_Er = 0.0d0
Av_EAu = 0.0d0
NTRAP = 0

!READ(40,*)Et_i,theta_i  ! Reading incidence energy + angle of incidence
!READ(50,*)NV
!READ(50,*)NVIB
ALLOCATE (rvNO(2,NVIB))
!READ(50,*)rvNO

WRITE(60,'(15X,A,6(18X,A))')'NTRAJ','xN','yN','zN','xO','yO','zO'
WRITE(70,'(15X,A,6(17X,A))')'NTRAJ','vxN','vyN','vzN','vxO','vyO','vzO'
WRITE(80,'(15X,A,6(18X,A))')'NTRAJ','xN','yN','zN','xO','yO','zO'
WRITE(90,'(15X,A,6(17X,A))')'NTRAJ','vxN','vyN','vzN','vxO','vyO','vzO'
WRITE(100,'(15X,A,13X,A,15X,A,15X,A,16X,A,7X,A)')'NTRAJ','E_TRANS','E_VIB','E_ROT','E_Au','SCATTER ANGLE'
WRITE(110,'(15X,A,16X,A,13X,A)')'NTRAJ','TIME','NBOUNCE'
WRITE(120,'(A50)')'SYNOPSIS'
WRITE(120,'(A50,I25)')'NUMBER OF TRAJECTORIES = ',NTRAJ
WRITE(120,'(A50,F25.16)')'LENGTH OF SIMULATION (ps) = ',DBLE(NSTEPS_DYN)*dt*0.01d0
WRITE(120,'(A50,F25.16)')'TIME STEP OF SIMULATION (fs) = ',dt*10.0d0
WRITE(120,*)''
WRITE(120,'(A50)')'INITIAL CONDITIONS:'
WRITE(120,'(A50,F25.16)')'TEMPERATURE OF Au(111) (K) = ',T
WRITE(120,*)''
WRITE(120,'(A50,F25.16)')'TRANSLATIONAL ENERGY OF NO (kJ/mol) = ',Et_i
WRITE(120,'(A50,I25)')'VIBRATIONAL STATE OF NO = ',NV                                                                                              
WRITE(120,'(A50,F25.16)')'ROTATIONAL ENERGY OF NO (kJ/mol) = ',0.0d0
WRITE(120,'(A50,F25.16)')'ANGLE OF INCIDENCE OF NO (degrees) = ',theta_i
!---------------------------------------------------------------------------------------------------

! BEGIN TRAJECTORIES:

CALL RANDOM_SEED
j = 0

DO j = 1, NTRAJ

x = xi ; v = 0.0d0 ; F = 0.0d0 ; F_next = 0.0d0 ; F_temp = 0.0d0 ; xiN = 0.0d0 ; xiO = 0.0d0 ; gauss_num = 0.0d0

!----------------------------------------------------------------------------------------------------
! BEGIN EQUILIBRATION OF Au(111)

! GENERATE INITIAL VELOCITIES OF Au ATOMS USING MAXWELL-BOLTZMANN DISTRIBUTION AT TEMPERATURE T

v(:,1:2) = 0.0d0 ! Set velocities of N and O equal to 0
CALL GAUSS(gauss_num)
v(:,3:N+2) = RESHAPE(gauss_num,(/3,N/)) * DSQRT(2.0d0*kB_MD*T/SPREAD(mass(3:N+2),1,3))
v(:,399:530) = 0.0d0 ! Set velocities of last layer of Au atoms equal to 0

! INITIALIZE ENERGIES AND FORCES OF Au ATOMS

CALL GetFNN2(F)
F(:,399:530) = 0.0d0  ! Set forces on of last layer of Au atoms equal to 0 i.e. fix last layer

i = 0

! VELOCITY VERLET ALGORITHM

DO i = 1, NSTEPS_EQ

  x(:,3:N+2) = x(:,3:N+2) + v(:,3:N+2)*dt + 0.5d0*(-F(:,3:N+2)/SPREAD(mass(3:N+2),1,3))*dt*dt
  CALL GetFNN2(F_next)
  F_next(:,399:530) = 0.0d0
  v(:,3:N+2) = v(:,3:N+2) + 0.5d0*dt*(-(F(:,3:N+2) + F_next(:,3:N+2))/SPREAD(mass(3:N+2),1,3))
  F = F_next

END DO

CALL GetVNN2(V_Au_Au)
KE = 0.5d0*SUM(SPREAD(mass(3:N+2),1,3)*v(:,3:N+2)*v(:,3:N+2))
E_Au_i = KE+V_Au_Au
!---------------------------------------------------------------------------------------------------

! READ INITIAL COORDINATES AND VELOCITIES OF NO

zcom = 12.0d0
Et = Et_i ; theta = theta_i
Et = Et*0.01d0  ! Converting incidence kinetic energy into MD units
CALL RANDOM_NUMBER(RAN) ; xcom = aPBC(1)*RAN ; CALL RANDOM_NUMBER(RAN) ; ycom = aPBC(2)*RAN
CALL RANDOM_NUMBER(RAN) ; polar = pi*RAN ; CALL RANDOM_NUMBER(RAN) ; azimuthal = 2.0d0*pi*RAN
CALL RANDOM_NUMBER(RAN) ; rNO = rvNO(1,IDINT(NVIB*RAN + 1.0d0)) ; vrel = rvNO(2,IDINT(NVIB*RAN + 1.0d0))
x(1,1) = xcom - (mass(2)*rNO*DSIN(polar)*DCOS(azimuthal)/(mass(1)+mass(2)))
x(2,1) = ycom - (mass(2)*rNO*DSIN(polar)*DSIN(azimuthal)/(mass(1)+mass(2)))
x(3,1) = zcom - (mass(2)*rNO*DCOS(polar)/(mass(1)+mass(2)))
x(1,2) = xcom + (mass(1)*rNO*DSIN(polar)*DCOS(azimuthal)/(mass(1)+mass(2)))
x(2,2) = ycom + (mass(1)*rNO*DSIN(polar)*DSIN(azimuthal)/(mass(1)+mass(2)))
x(3,2) = zcom + (mass(1)*rNO*DCOS(polar)/(mass(1)+mass(2)))

vcom = -DSQRT(2.0d0*Et/(mass(1)+mass(2)))
v(1,1) = vcom*DSIN(2.0d0*pi*theta/360.0d0) ; v(2,1) = 0.0d0 ; v(3,1) = vcom*DCOS(2.0d0*pi*theta/360.0d0)
v(:,2) = v(:,1)
v(1,1) = v(1,1) - (mass(2)*vrel*DSIN(polar)*DCOS(azimuthal)/(mass(1)+mass(2)))
v(2,1) = v(2,1) - (mass(2)*vrel*DSIN(polar)*DSIN(azimuthal)/(mass(1)+mass(2)))
v(3,1) = v(3,1) - (mass(2)*vrel*DCOS(polar)/(mass(1)+mass(2)))
v(1,2) = v(1,2) + (mass(1)*vrel*DSIN(polar)*DCOS(azimuthal)/(mass(1)+mass(2)))
v(2,2) = v(2,2) + (mass(1)*vrel*DSIN(polar)*DSIN(azimuthal)/(mass(1)+mass(2)))
v(3,2) = v(3,2) + (mass(1)*vrel*DCOS(polar)/(mass(1)+mass(2)))

! CALCULATE TRANSLATIONAL ENERGY AND EXIT ANGLE OF SCATTERED NO
  vcomsq = ((v(1,1)*mass(1)+v(1,2)*mass(2))/(mass(1)+mass(2)))**2&
  +((v(2,1)*mass(1)+v(2,2)*mass(2))/(mass(1)+mass(2)))**2+((v(3,1)*mass(1)+v(3,2)*mass(2))/(mass(1)+mass(2)))**2
  Et = 0.5d0*(mass(1)+mass(2))*vcomsq
  theta = DACOS(((v(3,1)*mass(1)+v(3,2)*mass(2))/(mass(1)+mass(2)))/DSQRT(vcomsq))
! CALCULATE VIBRATIONAL ENERGY OF SCATTERED NO
  CALL ENERGY_N_O_NEUTRAL(VIB_POT)
  VIB_KIN = 0.5d0*(mass(1)*mass(2)/(mass(1)+mass(2)))*(DOT_PRODUCT(x(:,2)-x(:,1),v(:,2)-v(:,1)))**2/(rNO*rNO)
  polar = DACOS((x(3,2)-x(3,1))/rNO)
  azimuthal = DACOS((x(1,2)-x(1,1))/(rNO*DSIN(polar)))
  Evib = VIB_KIN + VIB_POT
! CALCULATE ROTATIONAL ENERGY BY SUBTRACTION OF TRANSLATION AND VIBRATIONAL KINETIC FROM TOTAL KINETIC
  Er = (0.5d0*SUM(SPREAD(mass(1:2),1,3)*v(:,1:2)*v(:,1:2)) - Et - VIB_KIN)*100.0d0
! CALCULATE ENERGY OF Au(111)
  KE = 0.5d0*SUM(SPREAD(mass(3:N+2),1,3)*v(:,3:N+2)*v(:,3:N+2))

!---------------------------------------------------------------------------------------------------

! BEGIN SCATTERING OF NO ON Au(111)

i = 0

xiN(1:3,:) = x(:,3:N+2) - SPREAD(x(1:3,1),2,N)  ! Store instantaneous components of Au(i)--N distances in array xiN
xiO(1:3,:) = x(:,3:N+2) - SPREAD(x(1:3,2),2,N)  ! Store instantaneous components of Au(i)--O distances in array xiO
xiN(1:2,:) = xiN(1:2,:) - SPREAD(aPBC(1:2),2,N)*DNINT(xiN(1:2,:)/SPREAD(aPBC(1:2),2,N)) ! Apply PBC to x and y components of xiN
xiO(1:2,:) = xiO(1:2,:) - SPREAD(aPBC(1:2),2,N)*DNINT(xiO(1:2,:)/SPREAD(aPBC(1:2),2,N)) ! Apply PBC to x and y components of xiO
xiN(4,:) = DSQRT(SUM(xiN(1:3,:)*xiN(1:3,:),1))  ! Store instantaneous Au(i)--N distances in array xiN
xiO(4,:) = DSQRT(SUM(xiO(1:3,:)*xiO(1:3,:),1))  ! Store instantaneous Au(i)--O distances in array xiO
rNO = DSQRT(SUM((x(:,1)-x(:,2))*(x(:,1)-x(:,2)))) ! Calculate bond length of NO
zcom = ((x(3,1)*mass(1))+(x(3,2)*mass(2)))/(mass(1)+mass(2)) ! Calculate z distance of centre-of-mass of NO from Au(111) surface

CALL ENERGY_MATRIX(PE,EigenV) ! Calculate ground-state energy of NO--Au(111)
KE = 0.5d0*SUM(SPREAD(mass,1,3)*v*v)
CALL FORCE_MATRIX(F,EigenV)  ! Calculate ground-state forces on Au, N and O atoms
F(:,399:530) = 0.0d0
WRITE(60,'(I20,6F20.8)')j,x(1,1),x(2,1),x(3,1),x(1,2),x(2,2),x(3,2)
WRITE(70,'(I20,6F20.8)')j,v(1,1),v(2,1),v(3,1),v(1,2),v(2,2),v(3,2)                                                                                                                   

! VELOCITY VERLET ALGORITHM

CUTOFF = Au_N_cutoff_n  ! Cutoff distance of NO from Au(111) above which the trajectory will be stopped
NBOUNCE = 0             ! Counter for number of bounces initialized to 0
v1 = (v(3,1)*mass(1)+v(3,2)*mass(2))/(mass(1)+mass(2))  ! Instantaneous z-component of velocity of centre-of-mass of NO

DO i = 1, NSTEPS_DYN

  x = x + v*dt + 0.5d0*(F/SPREAD(mass,1,3))*dt*dt

  xiN(1:3,:) = x(:,3:N+2) - SPREAD(x(1:3,1),2,N)  
  xiO(1:3,:) = x(:,3:N+2) - SPREAD(x(1:3,2),2,N)  
  xiN(1:2,:) = xiN(1:2,:) - SPREAD(aPBC(1:2),2,N)*DNINT(xiN(1:2,:)/SPREAD(aPBC(1:2),2,N)) 
  xiO(1:2,:) = xiO(1:2,:) - SPREAD(aPBC(1:2),2,N)*DNINT(xiO(1:2,:)/SPREAD(aPBC(1:2),2,N)) 
  xiN(4,:) = DSQRT(SUM(xiN(1:3,:)*xiN(1:3,:),1))  
  xiO(4,:) = DSQRT(SUM(xiO(1:3,:)*xiO(1:3,:),1))  
  rNO = DSQRT(SUM((x(:,1)-x(:,2))*(x(:,1)-x(:,2)))) 
  zcom = ((x(3,1)*mass(1))+(x(3,2)*mass(2)))/(mass(1)+mass(2)) 

  CALL ENERGY_MATRIX(PE,EigenV)
  CALL FORCE_MATRIX(F_next,EigenV)
  F_next(:,399:530) = 0.0d0

  v = v + 0.5d0*dt*(F + F_next)/SPREAD(mass,1,3)
  F = F_next

  v2 = (v(3,1)*mass(1)+v(3,2)*mass(2))/(mass(1)+mass(2)) ! Calculating the number of times NO collides with Au(111)
  IF (v1*v2 < 0.0d0) NBOUNCE = NBOUNCE + 1
  v1 = v2

  IF ((v1 > 0.0d0).AND.(x(3,1) >= CUTOFF).AND.(x(3,2) >= CUTOFF)) EXIT

END DO

IF ((x(3,1) >= CUTOFF).AND.(x(3,2) >= CUTOFF)) THEN

  WRITE(80,'(I20,6F20.8)')j,x(1,1),x(2,1),x(3,1),x(1,2),x(2,2),x(3,2)
  WRITE(90,'(I20,6F20.8)')j,v(1,1),v(2,1),v(3,1),v(1,2),v(2,2),v(3,2)
  WRITE(110,'(I20,F20.8,I20)')j,DBLE(i)*dt*0.01d0,(NBOUNCE+1)/2

! CALCULATE TRANSLATIONAL ENERGY AND EXIT ANGLE OF SCATTERED NO
  vcomsq = ((v(1,1)*mass(1)+v(1,2)*mass(2))/(mass(1)+mass(2)))**2&
  +((v(2,1)*mass(1)+v(2,2)*mass(2))/(mass(1)+mass(2)))**2+((v(3,1)*mass(1)+v(3,2)*mass(2))/(mass(1)+mass(2)))**2
  Et = 0.5d0*(mass(1)+mass(2))*vcomsq
  theta = DACOS(((v(3,1)*mass(1)+v(3,2)*mass(2))/(mass(1)+mass(2)))/DSQRT(vcomsq))
  Av_Et = Av_Et + Et
  Av_theta = Av_theta + theta
  Av_weighted_theta = Av_weighted_theta + (theta/DSIN(theta))
  Norm_weighted_theta = Norm_weighted_theta + (1.0d0/DSIN(theta))
! CALCULATE VIBRATIONAL ENERGY OF SCATTERED NO
  CALL ENERGY_N_O_NEUTRAL(VIB_POT)
  VIB_KIN = 0.5d0*(mass(1)*mass(2)/(mass(1)+mass(2)))*(DOT_PRODUCT(x(:,2)-x(:,1),v(:,2)-v(:,1)))**2/(rNO*rNO)
  polar = DACOS((x(3,2)-x(3,1))/rNO)
  azimuthal = DACOS((x(1,2)-x(1,1))/(rNO*DSIN(polar)))
  Evib = VIB_KIN + VIB_POT
  Av_Evib = Av_Evib + Evib
! CALCULATE ROTATIONAL ENERGY BY SUBTRACTION OF TRANSLATION AND VIBRATIONAL KINETIC FROM TOTAL KINETIC
  Er = (0.5d0*SUM(SPREAD(mass(1:2),1,3)*v(:,1:2)*v(:,1:2)) - Et - VIB_KIN)
  Av_Er = Av_Er + Er
! CALCULATE ENERGY OF Au(111)
  CALL GetVNN2(V_Au_Au)
  KE = 0.5d0*SUM(SPREAD(mass(3:N+2),1,3)*v(:,3:N+2)*v(:,3:N+2))
  Av_EAu = Av_EAu + V_Au_Au + KE - E_Au_i

  WRITE(100,'(I20,5F20.8)')j,Et*100.0d0,Evib*100.0d0,Er*100.0d0,(V_Au_Au + KE - E_Au_i)*100.0d0,theta*(360.0d0/(2.0d0*pi))

ELSE

  NTRAP = NTRAP + 1

END IF

END DO

IF (NTRAP /= NTRAJ) THEN
Av_Et = Av_Et/DBLE(NTRAJ-NTRAP)
Av_Evib = Av_Evib/DBLE(NTRAJ-NTRAP)
Av_Er = Av_Er/DBLE(NTRAJ-NTRAP)
Av_EAu = Av_EAu/DBLE(NTRAJ-NTRAP)
Av_theta = Av_theta/DBLE(NTRAJ-NTRAP)
Av_weighted_theta = Av_weighted_theta/Norm_weighted_theta
END IF

WRITE(120,'(A50)')'AFTER SCATTERING'
WRITE(120,'(A50,I25)')'NUMBER TRAPPED = ',NTRAP
WRITE(120,'(A50,F25.16)')'TRAPPING PROBABILITY = ',DBLE(NTRAP)/DBLE(NTRAJ)
WRITE(120,'(A50,F25.16)')'AVERAGE TRANSLATIONAL ENERGY (kJ/mol) = ',Av_Et*100.0d0
WRITE(120,'(A50,F25.16)')'AVERAGE SCATTERING ANGLE (degrees) = ',Av_theta*(360.0d0/(2.0d0*pi))
WRITE(120,'(A50,F25.16)')'AVERAGE SCATTERING ANGLE (weighted, degrees) = ',Av_weighted_theta*(360.0d0/(2.0d0*pi))
WRITE(120,'(A50,F25.16)')'AVERAGE VIBRATIONAL ENERGY (kJ/mol) = ',Av_Evib*100.0d0
WRITE(120,'(A50,F25.16)')'AVERAGE ROTATIONAL ENERGY (kJ/mol) = ',Av_Er*100.0d0
WRITE(120,'(A50,F25.16)')'AVERAGE ENERGY TRANSFERRED TO Au(111) (kJ/mol) = ',Av_EAu*100.0d0

!---------------------------------------------------------------------------------------------------

DEALLOCATE(gauss_num,x,v,F,F_next,F_temp,mass,nn,xiN,xiO,F_NEUTRAL,F_ION,F_COUPLING,rvNO)

CLOSE(10)
CLOSE(20)
CLOSE(30)
CLOSE(40)
CLOSE(50)
CLOSE(60)
CLOSE(70)
CLOSE(80)
CLOSE(90)
CLOSE(100)

CONTAINS

!-----------------------------------------------------------

SUBROUTINE ENERGY_MATRIX(PE,EigenV)

IMPLICIT NONE

REAL(8), INTENT(OUT) :: PE, EigenV(2)  
REAL(8)     :: V_Au_N_n,V_Au_O_n,V_N_O_n   ! Components OF "Neutral" Au--NO energy
REAL(8)     :: V_Au_N_i,V_Au_O_i,V_N_O_i,V_image_pot   ! Components of "Ion" Au--NO energy
REAL(8)     :: En,Ei,Ec,Eg  ! "Neutral", "Ion", "Coupling", ground-state energies

!--- ENERGY OF THE "NEUTRAL" Au--NO SURFACE

CALL ENERGY_Au_N_NEUTRAL(V_Au_N_n)
CALL ENERGY_Au_O_NEUTRAL(V_Au_O_n)
CALL ENERGY_N_O_NEUTRAL(V_N_O_n)

En = V_Au_N_n + V_Au_O_n + V_N_O_n 

!--- ENERGY OF THE "ION" SURFACE

CALL ENERGY_Au_N_ION(V_Au_N_i)
CALL ENERGY_Au_O_ION(V_Au_O_i)
CALL ENERGY_N_O_ION(V_N_O_i)
CALL ENERGY_IMAGE_POTENTIAL(V_image_pot)

Ei = V_Au_N_i + V_Au_O_i + V_N_O_i + V_image_pot + WF - EA

!--- COUPLING FUNCTION

CALL ENERGY_COUPLING(Ec)

Eg = 0.5d0*(En+Ei-DSQRT((En-Ei)*(En-Ei)+(4.0d0*Ec*Ec)))

EigenV(1) = (Ei-Eg)/DSQRT((Ei-Eg)*(Ei-Eg)+(Ec*Ec))

EigenV(2) = -Ec/DSQRT((Ei-Eg)*(Ei-Eg)+(Ec*Ec))

PE = Eg 

END SUBROUTINE ENERGY_MATRIX

!---------------------------------------------------------

SUBROUTINE FORCE_MATRIX(F,EigenV)

IMPLICIT NONE

REAL(8), INTENT(OUT) :: F(:,:)
REAL(8), INTENT(IN)  :: EigenV(2)

F = 0.0d0 ; F_temp = 0.0d0
F_NEUTRAL = 0.0d0 ; F_ION = 0.0d0 ; F_COUPLING = 0.0d0

!--- FORCES FROM "NEUTRAL" SURFACE:

CALL FORCE_Au_N_NEUTRAL(F_temp)
F_NEUTRAL = F_temp

CALL FORCE_Au_O_NEUTRAL(F_temp)
F_NEUTRAL = F_NEUTRAL+F_temp

CALL FORCE_N_O_NEUTRAL(F_temp)
F_NEUTRAL = F_NEUTRAL+F_temp

!--- FORCES FROM "ION" SURFACE:

CALL FORCE_Au_N_ION(F_temp) ! *****Expensive
F_ION = F_temp

CALL FORCE_Au_O_ION(F_temp) ! *****Expensive
F_ION = F_ION + F_temp

CALL FORCE_N_O_ION(F_temp)
F_ION = F_ION + F_temp

CALL FORCE_IMAGE_POTENTIAL(F_temp)
F_ION = F_ION + F_temp

!--- FORCES DUE TO COUPLING:

CALL FORCE_COUPLING(F_COUPLING) ! *****Expensive

!--- FORCES DUE TO Au--Au INTERACTION:

CALL GetFNN2(F_temp)

F = EigenV(1)*EigenV(1)*F_NEUTRAL + EigenV(2)*EigenV(2)*F_ION + 2.0d0*EigenV(1)*EigenV(2)*F_COUPLING - F_temp

END SUBROUTINE FORCE_MATRIX

!---------------------------------------------------------------------

SUBROUTINE ENERGY_Au_N_NEUTRAL(V_Au_N_n)

IMPLICIT NONE

REAL(8), INTENT(OUT)  :: V_Au_N_n
REAL(8)  :: VCUT  ! Cutoff energy

VCUT = B_n*DEXP(-beta_n*Au_N_cutoff_n)

V_Au_N_n = SUM(B_n*DEXP(-beta_n*xiN(4,:)) - VCUT, MASK = xiN(4,:) < 10.0d0)

END SUBROUTINE ENERGY_Au_N_NEUTRAL

!------------------------------------------------------------

SUBROUTINE ENERGY_Au_O_NEUTRAL(V_Au_O_n)

IMPLICIT NONE

REAL(8), INTENT(OUT)  :: V_Au_O_n
REAL(8)  :: VCUT  ! Cutoff energy

VCUT = A_n*DEXP(-alpha_n*Au_O_cutoff_n)

V_Au_O_n = SUM(A_n*DEXP(-alpha_n*xiO(4,:)) - VCUT, MASK = xiO(4,:) < 10.0d0)

END SUBROUTINE ENERGY_Au_O_NEUTRAL

!---------------------------------------------------------

SUBROUTINE ENERGY_Au_N_ION(V_Au_N_i)

IMPLICIT NONE

REAL(8), INTENT(OUT)  :: V_Au_N_i
REAL(8)               :: cos_sq ! Cosine of the angle between the bond axis and the surface normal
REAL(8)  :: VCUT  ! Cutoff energy

cos_sq = (x(3,2)-x(3,1))*(x(3,2)-x(3,1))/(rNO*rNO)
VCUT = B_i*(DEXP(-2.0d0*beta_i*(Au_N_cutoff_i-rN_sur_e_i)) - 2.0d0*cos_sq*DEXP(-beta_i*(Au_N_cutoff_i-rN_sur_e_i)))

V_Au_N_i = SUM(B_i*(DEXP(-2.0d0*beta_i*(xiN(4,:)-rN_sur_e_i)) - 2.0d0*cos_sq*DEXP(-beta_i*(xiN(4,:)-rN_sur_e_i))) - VCUT,&
MASK = xiN(4,:) < 10.0d0)

END SUBROUTINE ENERGY_Au_N_ION

!-------------------------------------------------------------

SUBROUTINE ENERGY_Au_O_ION(V_Au_O_i)

IMPLICIT NONE

REAL(8), INTENT(OUT)  :: V_Au_O_i
REAL(8)  :: VCUT  ! Cutoff energy

VCUT = A_i*DEXP(-alpha_i*Au_O_cutoff_i)

V_Au_O_i = SUM(A_i*DEXP(-alpha_i*xiO(4,:)) - VCUT, MASK = xiO(4,:) < 10.0d0)

END SUBROUTINE ENERGY_Au_O_ION

!-----------------------------------------------------------------

SUBROUTINE ENERGY_IMAGE_POTENTIAL(V_image_pot)

IMPLICIT NONE

REAL(8), INTENT(OUT)  :: V_image_pot

V_image_pot = -D_i/DSQRT(C_i*C_i + (zcom - z0_i)*(zcom - z0_i))

END SUBROUTINE ENERGY_IMAGE_POTENTIAL

!-------------------------------------------------------------------

SUBROUTINE ENERGY_N_O_NEUTRAL(V_N_O_n)

IMPLICIT NONE

REAL(8), INTENT(OUT)  :: V_N_O_n

V_N_O_n = F_n*(1.0d0-DEXP(-delta_n*(rNO - rNO_e_n)))*(1.0d0-DEXP(-delta_n*(rNO - rNO_e_n)))

END SUBROUTINE ENERGY_N_O_NEUTRAL

!---------------------------------------------------------------------

SUBROUTINE ENERGY_N_O_ION(V_N_O_i)

IMPLICIT NONE

REAL(8), INTENT(OUT)  :: V_N_O_i

V_N_O_i = F_i*(1.0d0-DEXP(-delta_i*(rNO - rNO_e_i)))*(1.0d0-DEXP(-delta_i*(rNO - rNO_e_i)))

END SUBROUTINE ENERGY_N_O_ION

!-----------------------------------------------------------------------

SUBROUTINE ENERGY_COUPLING(V_coupling)

IMPLICIT NONE

REAL(8), INTENT(OUT)  :: V_coupling
REAL(8) :: FCUTN,FCUTO

FCUTN = coup_a_N/(1.0d0+coup_b_N*DEXP(coup_beta_N*Au_N_coupling_cutoff))
FCUTO = coup_a_O/(1.0d0+coup_b_O*DEXP(coup_beta_O*Au_O_coupling_cutoff))

V_coupling = SUM(coup_a_N/(1.0d0+coup_b_N*DEXP(coup_beta_N*xiN(4,:)))-FCUTN, MASK = xiN(4,:) < 10.0d0)
V_coupling = V_coupling + SUM(coup_a_O/(1.0d0+coup_b_O*DEXP(coup_beta_O*xiO(4,:)))-FCUTO, MASK = xiO(4,:) < 10.0d0)

END SUBROUTINE ENERGY_COUPLING

!--------------------------------------------------------------------------------

SUBROUTINE FORCE_Au_N_NEUTRAL(F_Au_N_n)

IMPLICIT NONE

REAL(8), INTENT(OUT)  :: F_Au_N_n(:,:)

F_Au_N_n = 0.0d0
         
F_Au_N_n(:,3:N+2) = -xiN(1:3,:)*SPREAD(-B_n*beta_n*DEXP(-beta_n*xiN(4,:))/xiN(4,:),1,3)

WHERE (SPREAD(xiN(4,:),1,3) >= 10.0d0)
 F_Au_N_n(:,3:N+2) = 0.0d0
END WHERE
 
F_Au_N_n(:,1) = SUM(-xiN(1:3,:)*B_n*beta_n*DEXP(-beta_n*SPREAD(xiN(4,:),1,3))/SPREAD(xiN(4,:),1,3),2,&
MASK = SPREAD(xiN(4,:),1,3) < 10.0d0)

END SUBROUTINE FORCE_Au_N_NEUTRAL 

!--------------------------------------------------------------------------

SUBROUTINE FORCE_Au_O_NEUTRAL(F_Au_O_n)

IMPLICIT NONE

REAL(8), INTENT(OUT)  :: F_Au_O_n(:,:)

F_Au_O_n = 0.0d0

F_Au_O_n(:,3:N+2) = -xiO(1:3,:)*SPREAD(-A_n*alpha_n*DEXP(-alpha_n*xiO(4,:))/xiO(4,:),1,3)

WHERE (SPREAD(xiO(4,:),1,3) >= 10.0d0)
 F_Au_O_n(:,3:N+2) = 0.0d0
END WHERE

F_Au_O_n(:,2) = SUM(-xiO(1:3,:)*A_n*alpha_n*DEXP(-alpha_n*SPREAD(xiO(4,:),1,3))/SPREAD(xiO(4,:),1,3),2,&
MASK = SPREAD(xiO(4,:),1,3) < 10.0d0)
         
END SUBROUTINE FORCE_Au_O_NEUTRAL 

!-----------------------------------------------------------------------------------

SUBROUTINE FORCE_Au_N_ION(F_Au_N_i)

IMPLICIT NONE

REAL(8), INTENT(OUT)  :: F_Au_N_i(:,:)
REAL(8), ALLOCATABLE  :: dr1(:,:) ! Difference between Au--N distance and equilibrium Au--N distance
REAL(8) :: dr2  ! Difference between Au--N cutoff distance and equilibrium Au--N distance
REAL(8) :: cos_sq ! Cosine of the angle between the bond axis and the surface normal
REAL(8) :: FCUT1(3) ! Component 1 of force due to cutoff in energy
REAL(8) :: FCUT2  ! Compoenent 2 of force due to cutoff in energy
         
ALLOCATE (dr1(3,3:N+2))

F_Au_N_i = 0.0d0
cos_sq = (x(3,2)-x(3,1))*(x(3,2)-x(3,1))/(rNO*rNO)

dr1 = 0.0d0
WHERE (SPREAD(xiN(4,:),1,3) < 10.0d0)
dr1 = SPREAD(xiN(4,:)-rN_sur_e_i,1,3)
END WHERE

dr2 = Au_N_cutoff_i - rN_sur_e_i
FCUT1(:) = B_i*4.0d0*DEXP(-beta_i*dr2)*(x(:,1)-x(:,2))*cos_sq/(rNO*rNO)
FCUT2 = B_i*4.0d0*DEXP(-beta_i*dr2)*(x(3,2) - x(3,1))/(rNO*rNO)

F_Au_N_i(:,1)=SUM(B_i*((4.0d0*DEXP(-beta_i*dr1)*SPREAD(x(:,1)-x(:,2),2,N)*cos_sq/(rNO*rNO))&
+(2.0d0*beta_i*xiN(1:3,:)*DEXP(-2.0d0*beta_i*dr1)/SPREAD(xiN(4,:),1,3))&
-(2.0d0*beta_i*xiN(1:3,:)*DEXP(-beta_i*dr1)*cos_sq/SPREAD(xiN(4,:),1,3)))-SPREAD(FCUT1,2,N),2,MASK = dr1 /= 0.0d0)

F_Au_N_i(3,1) = F_Au_N_i(3,1) + SUM((B_i*4.0d0*DEXP(-beta_i*dr1(3,:))*(x(3,2)-x(3,1))/(rNO*rNO))-FCUT2,MASK = dr1(3,:) /= 0.0d0)

F_Au_N_i(:,2) = SUM((-B_i*4.0d0*DEXP(-beta_i*dr1)*SPREAD(x(:,1)-x(:,2),2,N)*cos_sq/(rNO*rNO))+SPREAD(FCUT1,2,N),2,MASK = dr1 /= 0.0d0)

F_Au_N_i(3,2) = F_Au_N_i(3,2) + SUM((-B_i*4.0d0*DEXP(-beta_i*dr1(3,:))*(x(3,2)-x(3,1))/(rNO*rNO))+FCUT2,MASK = dr1(3,:) /= 0.0d0)

WHERE (dr1 /= 0.0d0)
F_Au_N_i(:,3:N+2) = B_i*((-2.0d0*beta_i*DEXP(-2.0d0*beta_i*dr1)*xiN(1:3,:)/SPREAD(xiN(4,:),1,3))&
+(2.0d0*beta_i*DEXP(-beta_i*dr1)*xiN(1:3,:)*cos_sq/SPREAD(xiN(4,:),1,3)))
END WHERE

F_Au_N_i = -F_Au_N_i

DEALLOCATE (dr1)

END SUBROUTINE FORCE_Au_N_ION

!------------------------------------------------------------------

SUBROUTINE FORCE_Au_O_ION(F_Au_O_i)

IMPLICIT NONE

REAL(8), INTENT(OUT)  :: F_Au_O_i(:,:)

F_Au_O_i = 0.0d0

F_Au_O_i(:,3:N+2) = -xiO(1:3,:)*SPREAD(-A_i*alpha_i*DEXP(-alpha_i*xiO(4,:))/xiO(4,:),1,3)

WHERE (SPREAD(xiO(4,:),1,3) >= 10.0d0)
 F_Au_O_i(:,3:N+2) = 0.0d0
END WHERE

F_Au_O_i(:,2) = SUM(-xiO(1:3,:)*A_i*alpha_i*DEXP(-alpha_i*SPREAD(xiO(4,:),1,3))/SPREAD(xiO(4,:),1,3),2,&
MASK = SPREAD(xiO(4,:),1,3) < 10.0d0)

END SUBROUTINE FORCE_Au_O_ION 

!-------------------------------------------------------------------
         
SUBROUTINE FORCE_IMAGE_POTENTIAL(F_image_pot) 

IMPLICIT NONE 

REAL(8), INTENT(OUT)  :: F_image_pot(:,:)
REAL(8)               :: FIMP  ! Component of force due to image potential

F_image_pot = 0.0d0

FIMP = D_i*(zcom - z0_i)/(C_i*C_i+(zcom - z0_i)*(zcom - z0_i))**1.5d0
F_image_pot(3,1) = -FIMP * mass(1)/(mass(1) + mass(2))
F_image_pot(3,2) = -FIMP * mass(2)/(mass(1) + mass(2)) 

END SUBROUTINE FORCE_IMAGE_POTENTIAL
 
!-------------------------------------------------------------------

SUBROUTINE FORCE_N_O_NEUTRAL(F_N_O_n)

IMPLICIT NONE

REAL(8), INTENT(OUT)  :: F_N_O_n(:,:)
REAL(8)               :: FNO ! Compenent of force due to N--O vibration on "neutral" surface

F_N_O_n = 0.0d0

FNO = 2.0d0*F_n*delta_n*(DEXP(-2.0d0*delta_n*(rNO - rNO_e_n))-DEXP(-delta_n*(rNO - rNO_e_n)))
F_N_O_n(:,1) = -FNO*(x(:,2)-x(:,1))/rNO
F_N_O_n(:,2) = -F_N_O_n(:,1)

END SUBROUTINE FORCE_N_O_NEUTRAL

!--------------------------------------------------------------------

SUBROUTINE FORCE_N_O_ION(F_N_O_i)

IMPLICIT NONE

REAL(8), INTENT(OUT)  :: F_N_O_i(:,:)
REAL(8)               :: FNO ! Compenent of force due to N--O vibration on "ion" surface

F_N_O_i = 0.0d0

FNO = 2.0d0*F_i*delta_i*(DEXP(-2.0d0*delta_i*(rNO - rNO_e_i))-DEXP(-delta_i*(rNO - rNO_e_i)))
F_N_O_i(:,1) = -FNO*(x(:,2)-x(:,1))/rNO
F_N_O_i(:,2) = -F_N_O_i(:,1)

END SUBROUTINE FORCE_N_O_ION

!---------------------------------------------------------------------

SUBROUTINE FORCE_COUPLING(F_COUPLING)

IMPLICIT NONE

REAL(8), INTENT(OUT)  :: F_COUPLING(:,:)
REAL(8), ALLOCATABLE  :: FCN(:,:),FCO(:,:)

ALLOCATE (FCN(3,3:N+2),FCO(3,3:N+2))

F_COUPLING = 0.0d0
FCN = 0.0d0
FCO = 0.0d0

WHERE (SPREAD(xiN(4,:),1,3) < 10.0d0)
FCN = coup_a_N*coup_b_N*coup_beta_N*DEXP(coup_beta_N*SPREAD(xiN(4,:),1,3))*xiN(1:3,:)/&
(SPREAD(xiN(4,:),1,3)*(1+coup_b_N*DEXP(coup_beta_N*SPREAD(xiN(4,:),1,3)))**2)
END WHERE
WHERE (SPREAD(xiO(4,:),1,3) < 10.0d0)
FCO = coup_a_O*coup_b_O*coup_beta_O*DEXP(coup_beta_O*SPREAD(xiO(4,:),1,3))*xiO(1:3,:)/&
(SPREAD(xiO(4,:),1,3)*(1+coup_b_O*DEXP(coup_beta_O*SPREAD(xiO(4,:),1,3)))**2)
END WHERE

F_COUPLING(:,1) = SUM(-FCN,2,MASK = FCN /= 0.0d0)
F_COUPLING(:,2) = SUM(-FCO,2,MASK = FCO /= 0.0d0)
F_COUPLING(:,3:N+2) = FCN + FCO

DEALLOCATE (FCN,FCO)

END SUBROUTINE FORCE_COUPLING

!---------------------------------------------------------------------------

SUBROUTINE GetNN2(nn)
! Calculates 1st nearest neighbor array for a set of atom positions in an fcc lattice
! Revised code based on coordinate transformation from 
! ([1 0 0],[0 1 0],[0 0 1]) --> ([1 0 -1]/sqrt(2)],[1 -4 1]/sqrt(6),[1 1 1]/sqrt(3))
! for use with [1 1 1] surface scattering simulations

IMPLICIT NONE

INTEGER, INTENT(OUT) :: nn(3:N+2,12) 		! nearest neighbors atoms, N x 12

REAL(8) :: r(3) 		! vector to nearest neighbor
REAL(8) :: rb(3)		! vector to nearest neighbor in Cartesian coordinates
REAL(8) :: TEMP(3)		
REAL(8)	:: U(3,3)	! Coordinate transormation matrix r(111) --> r(001)
INTEGER :: i,j,k,sint,jind
INTEGER :: s(3)

U = RESHAPE ( (/ -1d0/RT2,0d0,1d0/RT2, 1d0/RT6,-2d0/RT6,1d0/RT6, -1d0/RT3,-1d0/RT3,-1d0/RT3 /), (/3,3/) )
DO i = 3,N+2
DO j = 1,12
  nn(i,j) = 0
END DO
DO j = 3,N+2
  TEMP = (x(:,j) - x(:,i))/aPBC
  TEMP = TEMP - floor(TEMP + 0.5)
  r = TEMP*aPBC
IF (SUM(r*r) < dnn*dnn) THEN
rb = MATMUL(U,r)
s = NINT(rb/dnn)
sint = s(1)*9+s(2)*3+s(3)
SELECT CASE (sint)
case (4)
  nn(i,1) = j
case (2)
  nn(i,2) = j
case (10)
  nn(i,3) = j
case (-8)
  nn(i,4) = j
case (12)
  nn(i,5) = j
case (6)
  nn(i,6) = j
case (-4)
  nn(i,7) = j
case (-2)
  nn(i,8) = j
case (-10)
  nn(i,9) = j
case (8)
  nn(i,10) = j
case (-12)
  nn(i,11) = j
case (-6)
  nn(i,12) = j
END SELECT
END IF
END DO
END DO

END SUBROUTINE GetNN2
!------------------------------------------------------------------------------------
     
SUBROUTINE GetVNN2(Vme)

! Calculates potential energy for generalized force model
! Revised code based on coordinate transformation from 
! ([1 0 0],[0 1 0],[0 0 1]) --> ([1 0 -1]/sqrt(2)],[1 -4 1]/sqrt(6),[1 1 1]/sqrt(3))
! for use with [1 1 1] surface scattering simulations

IMPLICIT NONE

REAL(8), INTENT(OUT)	:: Vme		! Potential energy

REAL(8)			:: alpha(3) = (/-4.94d0, 17.15d0, 19.4d0/)     ! generalized force parameters
REAL(8)			:: r(3)		! Vector to nearest neighbors
REAL(8) 		:: TEMP(3)
REAL(8)			:: FmatOld(3,3)
REAL(8)			:: Fmat1(3,3),Fmat2(3,3),Fmat3(3,3),Fmat4(3,3),Fmat5(3,3),Fmat6(3,3)
					! Force matrices for each nearest neighbor position
REAL(8)			:: U(3,3),TU(3,3)	! Coordinate transormation matrix r(111) --> r(001)
INTEGER 		:: i,j,k
INTEGER			:: init = 0

SAVE init,Fmat1,Fmat2,Fmat3,Fmat4,Fmat5,Fmat6

alpha = alpha*6.02214199d-2

IF (init .EQ. 0) THEN
U = RESHAPE ( (/ -1d0/RT2,0d0,1d0/RT2, 1d0/RT6,-2d0/RT6,1d0/RT6, -1d0/RT3,-1d0/RT3,-1d0/RT3 /), (/3,3/) )
TU = TRANSPOSE(U)
FmatOld = RESHAPE ( (/ alpha(1),0d0,0d0,0d0,alpha(2),alpha(3),0d0,alpha(3),alpha(2) /), (/3,3/) )
Fmat1 = MATMUL(TU,MATMUL(FmatOld,U))
FmatOld = RESHAPE ( (/ alpha(1),0d0,0d0,0d0,alpha(2),-alpha(3),0d0,-alpha(3),alpha(2) /), (/3,3/) )
Fmat2 = MATMUL(TU,MATMUL(FmatOld,U))
FmatOld = RESHAPE ( (/ alpha(2),0d0,alpha(3),0d0,alpha(1),0d0,alpha(3),0d0,alpha(2) /), (/3,3/) )
Fmat3 = MATMUL(TU,MATMUL(FmatOld,U))
FmatOld = RESHAPE ( (/ alpha(2),0d0,-alpha(3),0d0,alpha(1),0d0,-alpha(3),0d0,alpha(2) /), (/3,3/) )
Fmat4 = MATMUL(TU,MATMUL(FmatOld,U))
FmatOld = RESHAPE ( (/ alpha(2),alpha(3),0d0,alpha(3),alpha(2),0d0,0d0,0d0,alpha(1) /), (/3,3/) )
Fmat5 = MATMUL(TU,MATMUL(FmatOld,U))
FmatOld = RESHAPE ( (/ alpha(2),-alpha(3),0d0,-alpha(3),alpha(2),0d0,0d0,0d0,alpha(1) /), (/3,3/) )
Fmat6 = MATMUL(TU,MATMUL(FmatOld,U))
init = 1
END IF


Vme = 0
DO i = 3,N+2
   DO j = 1,12
     IF (nn(i,j) /= 0) THEN
      TEMP = (x(:,nn(i,j))-x(:,i))/aPBC
      TEMP = TEMP - floor(TEMP + .5d0)
      r = TEMP*aPBC-r0(:,j)
      SELECT CASE (j)
         CASE (1)
	 Vme = Vme + DOT_PRODUCT(r,MATMUL(Fmat1,r))
         CASE (7)
	 Vme = Vme + DOT_PRODUCT(r,MATMUL(Fmat1,r))
         CASE (2)
	 Vme = Vme + DOT_PRODUCT(r,MATMUL(Fmat2,r))
	 CASE (8)
	 Vme = Vme + DOT_PRODUCT(r,MATMUL(Fmat2,r))
	 CASE (3)
	 Vme = Vme + DOT_PRODUCT(r,MATMUL(Fmat3,r))
	 CASE (9)
	 Vme = Vme + DOT_PRODUCT(r,MATMUL(Fmat3,r))
	 CASE (4)
	 Vme = Vme + DOT_PRODUCT(r,MATMUL(Fmat4,r))
	 CASE (10)
	 Vme = Vme + DOT_PRODUCT(r,MATMUL(Fmat4,r))
	 CASE (5)
	 Vme = Vme + DOT_PRODUCT(r,MATMUL(Fmat5,r))
	 CASE (11)
	 Vme = Vme + DOT_PRODUCT(r,MATMUL(Fmat5,r))
	 CASE (6)
	 Vme = Vme + DOT_PRODUCT(r,MATMUL(Fmat6,r))
	 CASE (12)
	 Vme = Vme + DOT_PRODUCT(r,MATMUL(Fmat6,r))
      END SELECT
      END IF
   END DO
END DO
Vme = Vme/4.0d0
END SUBROUTINE GetVNN2

!-------------------------------------------------------------------

SUBROUTINE GetFNN2(F_Au)

! Calculates forces for generalized force model
! Revised code based on coordinate transformation from 
! ([1 0 0],[0 1 0],[0 0 1]) --> ([1 0 -1]/sqrt(2)],[1 -4 1]/sqrt(6),[1 1 1]/sqrt(3))
! for use with [1 1 1] surface scattering simulations
! Definitions for nearest neighbors, forces, etc. can be found in G.H. Begbie, Proc. Roy. Soc. Lon.,
! (1947) 180-208.

IMPLICIT NONE

REAL(8), INTENT(OUT)	:: F_Au(:,:)	! forces on atoms

REAL(8)			:: alpha(3) = (/-4.94d0, 17.15d0, 19.4d0/)    ! generalized force parameters
REAL(8)			:: r(3)		! vector to nearest neighbors
REAL(8) 		:: TEMP(3)
REAL(8)			:: FmatOld(3,3)
REAL(8)			:: Fmat1(3,3),Fmat2(3,3),Fmat3(3,3),Fmat4(3,3),Fmat5(3,3),Fmat6(3,3)
					! Force matrices for each nearest neighbor position
REAL(8)			:: U(3,3),TU(3,3)	! Coordinate transormation matrix r(111) --> r(001)
INTEGER 		:: i,j,k
INTEGER			:: init = 0

SAVE init,Fmat1,Fmat2,Fmat3,Fmat4,Fmat5,Fmat6

alpha = alpha*6.02214199d-2

IF (init .EQ. 0) THEN
U = RESHAPE ( (/ -1d0/RT2,0d0,1d0/RT2, 1d0/RT6,-2d0/RT6,1d0/RT6, -1d0/RT3,-1d0/RT3,-1d0/RT3 /), (/3,3/) )
TU = TRANSPOSE(U)
FmatOld = RESHAPE ( (/ alpha(1),0d0,0d0,0d0,alpha(2),alpha(3),0d0,alpha(3),alpha(2) /), (/3,3/) )
Fmat1 = MATMUL(TU,MATMUL(FmatOld,U))
FmatOld = RESHAPE ( (/ alpha(1),0d0,0d0,0d0,alpha(2),-alpha(3),0d0,-alpha(3),alpha(2) /), (/3,3/) )
Fmat2 = MATMUL(TU,MATMUL(FmatOld,U))
FmatOld = RESHAPE ( (/ alpha(2),0d0,alpha(3),0d0,alpha(1),0d0,alpha(3),0d0,alpha(2) /), (/3,3/) )
Fmat3 = MATMUL(TU,MATMUL(FmatOld,U))
FmatOld = RESHAPE ( (/ alpha(2),0d0,-alpha(3),0d0,alpha(1),0d0,-alpha(3),0d0,alpha(2) /), (/3,3/) )
Fmat4 = MATMUL(TU,MATMUL(FmatOld,U))
FmatOld = RESHAPE ( (/ alpha(2),alpha(3),0d0,alpha(3),alpha(2),0d0,0d0,0d0,alpha(1) /), (/3,3/) )
Fmat5 = MATMUL(TU,MATMUL(FmatOld,U))
FmatOld = RESHAPE ( (/ alpha(2),-alpha(3),0d0,-alpha(3),alpha(2),0d0,0d0,0d0,alpha(1) /), (/3,3/) )
Fmat6 = MATMUL(TU,MATMUL(FmatOld,U))
init = 1
END IF

F_Au = 0.0d0

DO i = 3,N+2
   F_Au(:,i) = 0
   DO j = 1,12
     IF (nn(i,j) /= 0) THEN
      TEMP = (x(:,nn(i,j))-x(:,i))/aPBC
      TEMP = TEMP - floor(TEMP + .5d0)
      r = TEMP*aPBC-r0(:,j)
      SELECT CASE (j)
         CASE (1)
	 F_Au(:,i) = F_Au(:,i) - MATMUL(Fmat1,r)
         CASE (7)
	 F_Au(:,i) = F_Au(:,i) - MATMUL(Fmat1,r)
         CASE (2)
	 F_Au(:,i) = F_Au(:,i) - MATMUL(Fmat2,r)
	 CASE (8)
	 F_Au(:,i) = F_Au(:,i) - MATMUL(Fmat2,r)
	 CASE (3)
	 F_Au(:,i) = F_Au(:,i) - MATMUL(Fmat3,r)
	 CASE (9)
	 F_Au(:,i) = F_Au(:,i) - MATMUL(Fmat3,r)
	 CASE (4)
	 F_Au(:,i) = F_Au(:,i) - MATMUL(Fmat4,r)
	 CASE (10)
	 F_Au(:,i) = F_Au(:,i) - MATMUL(Fmat4,r)
	 CASE (5)
	 F_Au(:,i) = F_Au(:,i) - MATMUL(Fmat5,r)
	 CASE (11)
	 F_Au(:,i) = F_Au(:,i) - MATMUL(Fmat5,r)
	 CASE (6)
	 F_Au(:,i) = F_Au(:,i) - MATMUL(Fmat6,r)
	 CASE (12)
	 F_Au(:,i) = F_Au(:,i) - MATMUL(Fmat6,r)
      END SELECT
      END IF
   END DO
END DO

END SUBROUTINE GetFNN2

!-----------------------------------------------------------------------

SUBROUTINE GAUSS(gauss_num)

REAL(4) :: a0,a1,b0,b1
REAL(4), ALLOCATABLE :: smm(:),tmm(:),rmm(:)
REAL(8), ALLOCATABLE :: eps(:)
REAL(8), INTENT(OUT) :: gauss_num(:)

ALLOCATE (smm(3*N),tmm(3*N),rmm(3*N),eps(3*N))

DATA a0,a1,b0,b1/2.30753,.27061,.99229,.04481/

CALL RANDOM_NUMBER(eps)
smm = eps-0.5
tmm = sqrt(-alog(smm**2.0))
rmm = tmm-(a0+a1*tmm)/(1.+b0*tmm+b1*tmm*tmm)
gauss_num = sign(rmm,smm)

DEALLOCATE (smm,tmm,rmm,eps)

END SUBROUTINE GAUSS

!---------------------------------------------------------------------------

END PROGRAM NO_Au_111_2X2_PES_MD
