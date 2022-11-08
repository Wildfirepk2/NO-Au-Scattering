# custom Molly interactions/force functions

############################################################################################################

using Molly
using GLMakie

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
    # fetch atom_i's index
    i=atom_i.index

    # fetch system's force units. needed for Molly compatibility. may change later
    sysunits=sys_Au.force_units

    # no forces on atoms in last layer
    if i<auatomcutoff
        # normal operation

        # nn coords for atom i
        nncoords=[sys_Au.coords[j] for j in nn[i]]

        # nn coords dist away from atom i
        drij=[nncoords[j]-coord_i for j in eachindex(nncoords)]

        # dist of nn pairs away from equilibrium
        rij=[drij[j]-r0ij[i][:,j] for j in eachindex(drij)]

        # force each nn pair exerts on atom i
        Fij=[-Aijarray[i][j]*rij[j] for j in eachindex(rij)]

        # total force on atom i
        F=sum(Fij)

        # return force in system units.
        F .|> sysunits

        # v=getdrij(atom_i,coord_i)
    else
        # no force on Au atoms
        zeros(3)sysunits
    end
    # rij=[v[i]-r0ij[:,i] for i in 1:axes(r0ij,2)]

    # farray=[-Aijarray[i]*rij[i] for i in 1:eachindex(rij)]


    # # Replace this with your force calculation
    # # A positive force causes the atoms to move apart
    # f = 0.0

    # # Obtain a vector for the force
    # fdr = f * normalize(vec_ij)
    # return fdr

    # # get current nn distances for each atom. array of matrices
	# drij=getdrij(rAu)

	# # find distances away from equilibrium positions. array of matrices
	# rij=drij-r0ij

	# # temp array with forces on each atom due to nn's. array of arrays of arrays
	# farray=[getnnifarray(rij[i],i) for i in eachindex(rij)]

    # function getnnifarray(rij_i,i)
    #     # F=−Aij(nn case no)*rij
    #     [-Aijarray[i][j]*rij_i[:,j] for j in axes(rij_i,2)]
    # end

	# if inVfunc
	# 	# only if called in V_AuAu function. option added for speed
	# 	return rij,farray
	# else
	# 	# normal behavior. return only forces

	# 	# sum forces on each atom
	# 	f=sum.(farray)
	# 	return f
	# end
end

############################################################################################################

"""
output atom coords w time to csv
"""
function outputsyscoords()
    # coords from molly's loggers
    datasrc=sys_Au.loggers.coords.history

    # write csv files for ea time step
    for i in eachindex(datasrc)
        # xyz positions stored as temp vars
        xs=[datasrc[i][j][1] for j in eachindex(datasrc[i])]
        ys=[datasrc[i][j][2] for j in eachindex(datasrc[i])]
        zs=[datasrc[i][j][3] for j in eachindex(datasrc[i])]

        # write to csv
        data=DataFrame(x=xs,y=ys,z=zs)
        file="results/au slab equilibration/coords/syscoords-step $i.csv"
        CSV.write(file,data)
    end
end
