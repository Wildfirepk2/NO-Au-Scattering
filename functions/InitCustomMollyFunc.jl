# custom Molly interactions/force functions

############################################################################################################

using Molly

############################################################################################################

# interaction for Au slab equilibration
struct AuSlabInteraction <: PairwiseInteraction
    nl_only::Bool
end

############################################################################################################

# force function for Au slab equilibration. calculating F for each atom. inline/inbounds: copied from Lennard Jones force function. may also consider @fastmath
#= @inline @inbounds =# function Molly.force(inter::AuSlabInteraction,
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

        # nn coords for atom i. array of arrays
        nncoords=[sys_Au.coords[j] for j in nn[i]]

        # nn coords dist away from atom i. array of arrays
        drij=[vector(coord_i,nncoords[j],boundary) for j in eachindex(nncoords)]

        # dist of nn pairs away from equilibrium. array of arrays
        rij=[vector(r0ij[i][:,j],drij[j],boundary) for j in eachindex(drij)]

        # force each nn pair exerts on atom i. array of arrays
        Fij=[SVector{3}(-Aijarray[i][j]*rij[j]) for j in eachindex(rij)]

        # total force on atom i. array
        F=sum(Fij)

        # return force in system units.
        F .|> sysunits
    else
        # tmp step counter
        if i==au.N[1]
            println("Step $step_no")
            global step_no+=1
        end

        # no force on Au atoms
        SVector{3}(zeros(3)sysunits)
    end
end

############################################################################################################

# potential energy (V) for Au slab equilibration. calculating V for each atom. molly will sum all Vs. similar to force function. inline/inbounds: copied from Lennard Jones force function. may also consider @fastmath
#= @inline @inbounds =# function Molly.potential_energy(inter::AuSlabInteraction,
                          dr,
                          coord_i,
                          coord_j,
                          atom_i,
                          atom_j,
                          boundary)
    # fetch atom_i's index
    i=atom_i.index

    # # fetch system's force units. needed for Molly compatibility. may change later
    # sysunits=sys_Au.force_units

    # nn coords for atom i. array of arrays
    nncoords=[sys_Au.coords[j] for j in nn[i]]

    # nn coords dist away from atom i. array of arrays
    drij=[vector(coord_i,nncoords[j],boundary) for j in eachindex(nncoords)]

    # dist of nn pairs away from equilibrium. array of arrays
    rij=[vector(r0ij[i][:,j],drij[j],boundary) for j in eachindex(drij)]

    # V each nn pair exerts on atom i. array of numbers
    Vij=[dot(rij[j],Aijarray[i][j],rij[j]) for j in eachindex(rij)]

    # total V on atom i. divide by 2 because of v definition
    V=sum(Vij)/2
end
