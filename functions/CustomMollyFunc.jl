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
    fs[i] -= fdr_ustrip # negative force just works. investigate later
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

    # iterating through all NNs. may change later

    # # dist of nn pair away from equilibrium. static array
    # rij=vector(r0ij[i][:,idx_j],vec_ij,boundary)

    # # force nn pair exerts on atom i. static array
    # Fij=SVector{3}(-Aijarray[i][idx_j]*rij)

    # # tmp step counter. updates when loop reaches cutoff atom/394 pair
    # if i==auatomcutoff-1 && j==394
    #     println("Step $step_no")
    #     global step_no+=1
    # end

    # # return force in system units
    # Fij .|> sysunits

    # no forces on atoms in last layer
    if i<auatomcutoff
        # normal operation

        # dist of nn pair away from equilibrium. static array
        rij=vector(r0ij[i][:,idx_j],vec_ij,boundary)

        # force nn pair exerts on atom i. static array
        Fij=SVector{3}(-Aijarray[i][idx_j]*rij)

        # return force in system units
        Fij .|> sysunits
    else
        # tmp step counter. updates when loop reaches atom N/526 pair
        if i==au.N[1] && j==526
            println("Step $step_no")
            global step_no+=1
        end

        # no force on Au atoms
        SVector{3}(zeros(3)sysunits)
    end
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

# custom molly neighbor finder for Au slab
struct AuNeighborFinder
end

############################################################################################################

# neighbor finding function for Au neighbor finder. form copied from molly doc. just returns same neighbors each time.
function Molly.find_neighbors(s,
                        nf::AuNeighborFinder,
                        current_neighbors=nothing,
                        step_n::Integer=0;
                        n_threads::Integer=Threads.nthreads())
    nn_molly
end

############################################################################################################

# interaction for NO/Au scattering
struct NOAuInteraction <: PairwiseInteraction
    nl_only::Bool
end

############################################################################################################

# custom molly neighbor finder for NO/Au scattering
struct NONeighborFinder{D1,D2}
    N_cutoff::D1
    O_cutoff::D2
end

############################################################################################################

# helper function to construct NONeighborFinder w named parameters
function NONeighborFinder(;
                        N_cutoff,
                        O_cutoff)
    return NONeighborFinder{typeof(N_cutoff),typeof(O_cutoff)}(N_cutoff, O_cutoff)
end

############################################################################################################

# neighbor finding function for NO/Au neighbor finder. form copied from molly doc. just returns same neighbors each time.
function Molly.find_neighbors(s,
                        nf::NONeighborFinder,
                        current_neighbors=nothing,
                        step_n::Integer=0;
                        n_threads::Integer=Threads.nthreads())
    # stripping units for NearestNeighbors compatibility. standardizing units to Angstroms
    Ncut_nounit=ustrip(u"Å",nf.N_cutoff)
    Ocut_nounit=ustrip(u"Å",nf.O_cutoff)
    PBC_nounit=ustrip.(u"Å",s.boundary)
	PBC=PeriodicEuclidean(PBC_nounit) # simulation box with periodic boundary conditions
    coords_nounit=ustrip_vec.(u"Å",s.coords)

    # initializing coordinates for nearest neighbor search. BallTree for PBC compatibility. tested to be same as BruteTree (purely distance searching)
	tree=BallTree(coords_nounit, PBC)

    # search for all points within distances
	nn_N=inrange(tree,coords_nounit[1],Ncut_nounit)
    nn_O=inrange(tree,coords_nounit[2],Ocut_nounit)
    nn=[nn_N,nn_O]

    # remove same atom pairs from nn list
	removeiipairs!(nn)

    # nn list in tuple form for molly compatibility. (atom i, atom j, weight 14=false (0)). then pass to molly neighbor object
	nntuples=[(i,nn[i][j],false) for i in eachindex(nn) for j in eachindex(nn[i])]

    # molly neighbor object
	nn_molly=NeighborList(length(nntuples), nntuples)

    # above looks ok

############################################################################################################
dist_unit = unit(first(first(s.coords)))
    bv = ustrip.(dist_unit, s.boundary)
    btree = BallTree(ustrip_vec.(s.coords), PeriodicEuclidean(bv))
    dist_cutoff = ustrip(dist_unit, nf.dist_cutoff)

    # @floop ThreadedEx(basesize = length(s) ÷ n_threads) for i in 1:length(s)
    #     ci = ustrip.(s.coords[i])
    #     nbi = @view nf.nb_matrix[:, i]
    #     w14i = @view nf.matrix_14[:, i]
    #     idxs = inrange(btree, ci, dist_cutoff, true)
    #     for j in idxs
    #         if nbi[j] && i > j
    #             nn = (i, j, w14i[j])
    #             @reduce(neighbors_list = append!(Tuple{Int, Int, Bool}[], (nn,)))
    #         end
    #     end
    # end

    return NeighborList(length(neighbors_list), neighbors_list)

    ############################################################################################################









    rAu_nounit=ustrip.(u"Å",rAu)
	dnn_nounit=ustrip(u"Å",au.dnn[1])
	PBC_nounit=ustrip.(u"Å",simboxdims)
	PBC=PeriodicEuclidean(PBC_nounit) # simulation box with periodic boundary conditions

	# initializing coordinates for nearest neighbor search. BallTree for PBC compatibility. tested to be same as BruteTree (purely distance searching)
	tree=BallTree(rAu_nounit, PBC)

	# search for all points within distance r. +1 for rounding effects
	r=dnn_nounit+1
	nn=inrange(tree,rAu_nounit,r)

	# remove same atom pairs from nn list
	removeiipairs!(nn)

	# nn list in tuple form for molly compatibility. (atom i, atom j, weight 14=false (0)). then pass to molly neighbor object
	nntuples=[(i,nn[i][j],false) for i in eachindex(nn) for j in eachindex(nn[i])]
	nn_molly=NeighborList(length(nntuples), nntuples)

	# output as 528 length array of arrays and molly neighbor object
	return nn,nn_molly
end
