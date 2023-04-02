# custom Molly interactions/force functions for O/Au

############################################################################################################

# interaction for O/Au scattering
struct OAuInteraction <: PairwiseInteraction
    nl_only::Bool
end

############################################################################################################

# custom molly neighbor finder for O/Au scattering
struct ONeighborFinder{D1}
    O_cutoff::D1
end

############################################################################################################

# helper function to construct ONeighborFinder w named parameters
function ONeighborFinder(;
                        O_cutoff)
    return ONeighborFinder{typeof(O_cutoff)}(O_cutoff)
end

############################################################################################################

# neighbor finding function for O/Au neighbor finder. form copied from molly doc. just returns same neighbors each time.
function Molly.find_neighbors(s,
                        nf::ONeighborFinder,
                        current_neighbors=nothing,
                        step_n::Integer=0;
                        n_threads::Integer=Threads.nthreads())
    # stripping units for NearestNeighbors compatibility. standardizing units to Angstroms
    Ocut_nounit=ustrip(u"Å",nf.O_cutoff)
    PBC_nounit=ustrip.(u"Å",s.boundary)
	PBC=PeriodicEuclidean(PBC_nounit) # simulation box with periodic boundary conditions
    coords_nounit=ustrip_vec.(u"Å",s.coords)

    # initializing coordinates for nearest neighbor search. BallTree for PBC compatibility. tested to be same as BruteTree (purely distance searching)
	tree=BallTree(coords_nounit, PBC)

    # search for all points within distances. type: Vector{Int64}
    nn_O=inrange(tree,coords_nounit[1],Ocut_nounit)

    # nn list in tuple form for molly compatibility. (atom i, atom j, weight 14=false (0)). then pass to molly neighbor object
	nntuples_o=[(1,nn_O[i],false) for i in eachindex(nn_O)]
    nntuples=vcat(nntuples_o,nntuples_oau)

    # molly neighbor object
	NeighborList(length(nntuples), nntuples)

    # above looks ok

    # parallelization
    # dist_unit = unit(first(first(s.coords)))
    # bv = ustrip.(dist_unit, s.boundary)
    # btree = BallTree(ustrip_vec.(s.coords), PeriodicEuclidean(bv))
    # dist_cutoff = ustrip(dist_unit, nf.dist_cutoff)

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

    # return NeighborList(length(neighbors_list), neighbors_list)
end

############################################################################################################

# alter molly PE function above for O/Au. need to calculate E all at once
function Molly.potential_energy(s::System{D, false, T, CU, A, AD, PI} where {D,T,CU,A,AD,PI<:Tuple{OAuInteraction}}, neighbors=nothing)
    En=0 * s.energy_units
    Ei=0 * s.energy_units
    Ec=0 * s.energy_units
    EAu=0 * s.energy_units

    Ocoord=s.coords[1]
    dz=Ocoord[3]-maximum(au.z)

    @inbounds for ni in 1:neighbors.n
        i, j, weight_14 = neighbors.list[ni]
        # i<3: NOT factoring in V_AuAu.\fix?
        # if i<3
            # setup for E calc
            bound=s.boundary
            dr = vector(s.coords[i], s.coords[j], bound)
            distbtwn=euclidean(dr,zeros(3)u"Å")
            
            t_En,t_Ei,t_Ec,t_EAu=getVij_OAu(i,j,distbtwn,dr,dz,bound)
            En+=t_En
            Ei+=t_Ei
            Ec+=t_Ec
            EAu+=t_EAu
        # end
    end

    # final ground state energy/eigenvalues
    Eg=1/2*(En+Ei-√((En-Ei)^2+4Ec^2))
    λ1=(Ei-Eg)/√((Ei-Eg)^2+Ec^2)
    λ2=-Ec/√((Ei-Eg)^2+Ec^2)

    # store eigenvalues for force calculation
    df=DataFrame(Eg=Eg,λ1=λ1,λ2=λ2)
    append!(storeEs,df)

    # case for neu+ion?
    if neutral_PES_active && ionic_PES_active && coupled_PES_active
        return uconvert(s.energy_units, Eg+EAu)
    elseif ionic_PES_active
        return uconvert(s.energy_units, Ei+EAu)
    else
        return uconvert(s.energy_units, En+EAu)
    end
end

############################################################################################################

# mostly a copy of molly force! (level above force). needed to freeze Au atoms in scattering FOR NOW
@inline @inbounds function Molly.force!(fs, inter::OAuInteraction, s::System, i::Integer, j::Integer, force_units, weight_14::Bool=false)
    dr = vector(s.coords[i], s.coords[j], s.boundary)
    fdr = force(s, inter, dr, s.coords[i], s.coords[j], s.atoms[i], s.atoms[j], s.boundary)
    Molly.check_force_units(fdr, force_units)
    fdr_ustrip = ustrip.(fdr)
    fs[i] -= fdr_ustrip
    # reactive force on Au atoms
    if j<auatomcutoff
        fs[j] += fdr_ustrip
    end
    return nothing
end

############################################################################################################

# force function for O/Au scattering. calculating F for each atom due to NN. 
# inline/inbounds: copied from Lennard Jones force function. may also consider @fastmath
# see fortran FORCE_MATRIX
@inline @inbounds function Molly.force(s::System,
                            inter::OAuInteraction,
                            vec_ij,
                            coord_i,
                            coord_j,
                            atom_i,
                            atom_j,
                            boundary)
    # fetch atom indices
    i=atom_i.index
    j=atom_j.index
    sysunits=s.force_units

    # force calculation
    Ocoord=s.coords[1]
    dz=Ocoord[3]-maximum(au.z)
    F=getFij_OAu(i,j,vec_ij,dz,boundary)
    F .|> sysunits
end

############################################################################################################

# mod molly simulate! to fix force fn lagging behind eigenvalues
function Molly.simulate!(sys::System{D, false, T, CU, A, AD, PI} where {D,T,CU,A,AD,PI<:Tuple{OAuInteraction}},
                    sim::VelocityVerlet,
                    n_steps::Integer;
                    n_threads::Integer=Threads.nthreads())
    sys.coords = wrap_coords.(sys.coords, (sys.boundary,))
    sim.remove_CM_motion && remove_CM_motion!(sys)
    neighbors = find_neighbors(sys, sys.neighbor_finder; n_threads=n_threads)
    run_loggers!(sys, neighbors, 0; n_threads=n_threads)
    accels_t = accelerations(sys, neighbors; n_threads=n_threads)
    accels_t_dt = zero(accels_t)

    for step_n in 1:n_steps
        old_coords = copy(sys.coords)
        sys.coords += sys.velocities .* sim.dt .+ (Molly.remove_molar.(accels_t) .* sim.dt ^ 2) ./ 2

        apply_constraints!(sys, old_coords, sim.dt)
        sys.coords = wrap_coords.(sys.coords, (sys.boundary,))

        potential_energy(sys, neighbors) # modded part. forces \lambda to be pushed to storeEs
        accels_t_dt = accelerations(sys, neighbors; n_threads=n_threads)

        sys.velocities += Molly.remove_molar.(accels_t .+ accels_t_dt) .* sim.dt / 2

        sim.remove_CM_motion && remove_CM_motion!(sys)
        apply_coupling!(sys, sim.coupling, sim)

        run_loggers!(sys, neighbors, step_n; n_threads=n_threads)

        if step_n != n_steps
        neighbors = find_neighbors(sys, sys.neighbor_finder, neighbors, step_n;
                                n_threads=n_threads)
        accels_t = accels_t_dt
        end
    end
    return sys
end