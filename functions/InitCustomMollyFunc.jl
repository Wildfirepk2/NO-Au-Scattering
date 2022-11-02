# custom Molly interactions/force functions

############################################################################################################

using Molly

############################################################################################################

# interaction for Au slab equilibration
struct AuSlabInteraction <: PairwiseInteraction
    nl_only::Bool
    # Any other properties, e.g. constants for the interaction or cutoff parameters
end

############################################################################################################

# force function for Au slab equilibration
function Molly.force(inter::AuSlabInteraction,
                        vec_ij,
                        coord_i,
                        coord_j,
                        atom_i,
                        atom_j,
                        boundary)
    i=atom_i.index
    if i<auatomcutoff
        # normal operation
        v=getdrij(atom_i,coord_i)
    else
        # no force on Au atoms
        return zeros(3)u"N/mol"
    end
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