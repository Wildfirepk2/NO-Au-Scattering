# adiabatic molecular dynamics of NO scattering off of an Au(111) surface

# index
# roy art: Roy, S.; Shenvi, N. A.; Tully, J. C. Model Hamiltonian for the Interaction of NO with the Au(111) Surface. J. Chem. Phys. 2009, 130 (17), 174716. https://doi.org/10.1063/1.3122989.
# begbie: Begbie, G. H.; Born, M. Thermal Scattering of X-Rays by Crystals II. The Thermal Scattering of the Face-Centred Cubic and the Close-Packed Hexagonal Lattices. Proceedings of the Royal Society of London. Series A. Mathematical and Physical Sciences 1947, 188 (1013), 189–208. https://doi.org/10.1098/rspa.1947.0004.
# fortran: Dr. Roy's fortran program for roy art
# notebook: personal notebook, no-au scattering notebook.docx

# checked: 10/27/22

############################################################################################################

# initialize all functions
include("functions/InitMDUnits.jl")
include("functions/InitAllParams.jl")
include("functions/InitAllAuNOVFunc.jl")
include("functions/InitAllAuVFunc.jl")

############################################################################################################

# parameter variables (type: DataFrame) from input files
au=initAuParams()
param=initSimParams()
PES_GS=initPESParamGS()
PES_ionic=initPESParamIonic()
PES_coup=initPESParamCoup()
no=initNOParams()

############################################################################################################

# main variables

# 3 x N matrix. xyz coords of each Au atom stored in columns
rAu=initAuCoords()

# array of arrays. nearest neighbors for each Au atom. ith row corresponds to Au atom i's nearest neighbors (in terms of atom number)
nn=getnn()

# array of matrices. initial distances to nearest neighbors for each atom. nn xyz coords stored in columns
r0ij=getdrij(rAu)

# array of arrays. case number for nearest neighbors for each atom. 
nncase=getcase()

# Dict mapping nn case number to force matrix.
Aij=initAij()

# array of arrays of matrices. force matrix for nearest neighbors for each atom. 
Aijarray=initAijarray()

############################################################################################################

# other variables
mAu=196.966569u"u"
mN=14.00674u"u"
mO=15.9994u"u"

Av_Et = 0
Av_theta = 0
Av_weighted_theta = 0
Norm_weighted_theta = 0
Av_Evib = 0
Av_Er = 0
Av_EAu = 0
NTRAP = 0

############################################################################################################

# custom Molly interactions

struct AuSlabInteraction <: PairwiseInteraction
    nl_only::Bool
    # Any other properties, e.g. constants for the interaction or cutoff parameters
end

function Molly.force(inter::AuSlabInteraction,
                        vec_ij,
                        coord_i,
                        coord_j,
                        atom_i,
                        atom_j,
                        boundary)
    v=getdrij(atom_i,coord_i)
    rij=[v[i]-r0ij[:,i] for i in 1:axes(r0ij,2)]

    farray=[-Aijarray[i]*rij[i] for i in 1:eachindex(rij)]


    # Replace this with your force calculation
    # A positive force causes the atoms to move apart
    f = 0.0

    # Obtain a vector for the force
    fdr = f * normalize(vec_ij)
    return fdr

    # get current nn distances for each atom. array of matrices
	drij=getdrij(rAu)

	# find distances away from equilibrium positions. array of matrices
	rij=drij-r0ij

	# temp array with forces on each atom due to nn's. array of arrays of arrays
	farray=[getnnifarray(rij[i],i) for i in eachindex(rij)]

    function getnnifarray(rij_i,i)
        # F=−Aij(nn case no)*rij
        [-Aijarray[i][j]*rij_i[:,j] for j in axes(rij_i,2)]
    end

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

# start Molly simulation for Au slab equilibration

n_atoms = au.N[1]
boundary = CubicBoundary(au.aPBCx[1], au.aPBCy[1], Inf*u"Å")
temp = param.T[1]
atom_mass = mAu
nsteps=param.Nsteps_eq[1]

atoms = [Atom(index=i, mass=atom_mass) for i in 1:n_atoms]
coords = [SA[au.x[i],au.y[i],-au.z[i]] for i in 1:n_atoms]
velocities = [velocity(atom_mass, temp) for i in 1:n_atoms]
pairwise_inters = ()
simulator = VelocityVerlet(
    dt=param.dt[1],
    coupling=AndersenThermostat(temp, 1.0u"10fs"),
)

initcoords=copy(coords)

sys = System(
    atoms=atoms,
    pairwise_inters=pairwise_inters,
    coords=coords,
    velocities=velocities,
    boundary=boundary,
    loggers=(
        # temp=TemperatureLogger(100),
        coords=CoordinateLogger(10),
    ),
)

simulate!(sys, simulator, nsteps)
visualize(sys.loggers.coords, boundary, "au slab equilibration.mp4")
