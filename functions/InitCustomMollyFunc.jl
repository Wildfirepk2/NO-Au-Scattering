# custom Molly interactions/force functions

############################################################################################################

# interaction for Au slab equilibration
struct AuSlabInteraction <: PairwiseInteraction
    nl_only::Bool
end

############################################################################################################

# mostly a copy of molly force! (level above force). needed to calculate force correctly for my system
@inline @inbounds function Molly.force!(fs, inter::AuSlabInteraction, s::System, i::Integer, j::Integer, force_units, weight_14::Bool=false)
    dr = vector(s.coords[i], s.coords[j], s.boundary)
    fdr = force(inter, dr, s.coords[i], s.coords[j], s.atoms[i], s.atoms[j], s.boundary, weight_14)
    Molly.check_force_units(fdr, force_units)
    fdr_ustrip = ustrip.(fdr)
    fs[i] += fdr_ustrip
    # println("using custom force function")
    return nothing
end

############################################################################################################


# force function for Au slab equilibration. calculating F for each atom due to NN. 
# inline/inbounds: copied from Lennard Jones force function. may also consider @fastmath
# see roy art, p7, eq 20. see fortran GetFNN2
@inline @inbounds function Molly.force(inter::AuSlabInteraction,
                        vec_ij,
                        coord_i,
                        coord_j,
                        atom_i,
                        atom_j,
                        boundary)
    # fetch atom indices
    i=atom_i.index
    j=atom_j.index

    # find j in nn[i] and return index. needed to access items in variables based off nn
    idx_j=findfirst(isequal(j), nn[i])

    # fetch system's force units. needed for Molly compatibility. may change later
    sysunits=sys_Au.force_units

    # dist of nn pair away from equilibrium. static array
    rij=vector(r0ij[i][:,idx_j],vec_ij,boundary)

    # force nn pair exerts on atom i. static array
    Fij=SVector{3}(-Aijarray[i][idx_j]*rij)

    # tmp step counter. updates when loop reaches cutoff atom/394 pair
    if i==auatomcutoff-1 && j==394
        println("Step $step_no")
        global step_no+=1
    end

    # return force in system units
    Fij .|> sysunits

    # # no forces on atoms in last layer
    # if i<auatomcutoff
    #     # normal operation

    #     # dist of nn pair away from equilibrium. static array
    #     rij=vector(r0ij[i][:,idx_j],vec_ij,boundary)

    #     # force nn pair exerts on atom i. static array
    #     Fij=SVector{3}(-Aijarray[i][idx_j]*rij)

    #     # return force in system units
    #     Fij .|> sysunits
    # else
    #     # tmp step counter. updates when loop reaches atom N/526 pair
    #     if i==au.N[1] && j==526
    #         println("Step $step_no")
    #         global step_no+=1
    #     end

    #     # no force on Au atoms
    #     SVector{3}(zeros(3)sysunits)
    # end
end

############################################################################################################


# potential energy (V) for Au slab equilibration. calculating V for each atom due to NN. molly will sum all Vs. similar to force function. 
# inline/inbounds: copied from Lennard Jones force function. may also consider @fastmath. 
# due to truncated NN list, doesnt calc V for back layer. ok since V would be constant (layer doesnt move)
# see roy art, p7, eq 20. see fortran GetVNN2
@inline @inbounds function Molly.potential_energy(inter::AuSlabInteraction,
                          dr,
                          coord_i,
                          coord_j,
                          atom_i,
                          atom_j,
                          boundary)
    # fetch atom indices
    i=atom_i.index
    j=atom_j.index

    # find j in nn[i] and return index. needed to access items in variables based off nn
    idx_j=findfirst(isequal(j), nn[i])

    # dist of nn pair away from equilibrium. static array
    rij=vector(r0ij[i][:,idx_j],dr,boundary)

    # V nn pair exerts on atom i. number. divide by 2 because of v definition
    dot(rij,Aijarray[i][idx_j],rij)/2
end

############################################################################################################

# custom neighbor finder for molly compatibility
struct AuNeighborFinder
end

############################################################################################################


# neighbor finding function for custom neighbor finder. form copied from molly doc. just returns same neighbors each time.
function Molly.find_neighbors(s,
                        nf::AuNeighborFinder,
                        current_neighbors=nothing,
                        step_n::Integer=0;
                        n_threads::Integer=Threads.nthreads())
    nn_molly
end
