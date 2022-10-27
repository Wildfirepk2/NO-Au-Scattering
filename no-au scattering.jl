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