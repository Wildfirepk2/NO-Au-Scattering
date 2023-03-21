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
    fdr = force(s, inter, dr, s.coords[i], s.coords[j], s.atoms[i], s.atoms[j], s.boundary)
    Molly.check_force_units(fdr, force_units)
    fdr_ustrip = ustrip.(fdr)
    fs[i] -= fdr_ustrip # negative force just works. investigate later
    if j<auatomcutoff
        fs[j] += fdr_ustrip
    end
    return nothing
end

############################################################################################################


# force function for Au slab equilibration. calculating F for each atom due to NN. 
# inline/inbounds: copied from Lennard Jones force function. may also consider @fastmath
# see roy art, p7, eq 20. see fortran GetFNN2
@inline @inbounds function Molly.force(s::System,
                        inter::AuSlabInteraction,
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
    sysunits=s.force_units

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
	nntuples_no=[(i,nn[i][j],false) for i in eachindex(nn) for j in eachindex(nn[i])]
    nntuples=vcat(nntuples_no,nntuples_au)

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

# alter molly PE function above for NO/Au. need to calculate E all at once. doesnt have V_AuAu for speed. plus not in fortran
function Molly.potential_energy(s::System{D, false, T, CU, A, AD, PI} where {D,T,CU,A,AD,PI<:Tuple{NOAuInteraction}}, neighbors=nothing)
    En=0 * s.energy_units
    Ei=0 * s.energy_units
    Ec=0 * s.energy_units
    EAu=0 * s.energy_units

    @inbounds for ni in 1:neighbors.n
        i, j, weight_14 = neighbors.list[ni]
        # i<3: NOT factoring in V_AuAu.\fix?
        # if i<3
            # setup for E calc
            bound=s.boundary
            dr = vector(s.coords[i], s.coords[j], bound)
            distbtwn=euclidean(dr,zeros(3)u"Å")

            Ncoords=s.coords[1]
            Ocoords=s.coords[2]
            mN=s.atoms[1].mass
            mO=s.atoms[2].mass
            cosθ=getcosth(Ncoords,Ocoords,bound)
            dz=getzcom(mN,mO,Ncoords,Ocoords)

            t_En,t_Ei,t_Ec,t_EAu=getVij_NOAu(i,j,distbtwn,cosθ,dz,dr,bound)
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
    push!(storeEs,[Eg,λ1,λ2])

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
@inline @inbounds function Molly.force!(fs, inter::NOAuInteraction, s::System, i::Integer, j::Integer, force_units, weight_14::Bool=false)
    dr = vector(s.coords[i], s.coords[j], s.boundary)
    fdr = force(s, inter, dr, s.coords[i], s.coords[j], s.atoms[i], s.atoms[j], s.boundary)
    Molly.check_force_units(fdr, force_units)
    fdr_ustrip = ustrip.(fdr)
    fs[i] -= fdr_ustrip
    # reactive force on Au atoms/NO reactive force for (1,2) or (2,1) pairs. need those calculated individually
    if j>2 && j<auatomcutoff
        fs[j] += fdr_ustrip
    end
    return nothing
end

############################################################################################################

# force function for NO/Au scattering. calculating F for each atom due to NN. 
# inline/inbounds: copied from Lennard Jones force function. may also consider @fastmath
# see fortran FORCE_MATRIX
@inline @inbounds function Molly.force(s::System,
                            inter::NOAuInteraction,
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

    # get eigenvalues
    λ1=storeEs[end,2]
    λ2=storeEs[end,3]
    a=λ1^2
    b=λ2^2
    c=2*λ1*λ2

    # setup for force calculation
    distbtwn=euclidean(vec_ij,zeros(3)u"Å")
    u=normalize(vec_ij)

    Ncoords=s.coords[1]
    Ocoords=s.coords[2]
    vecON=vector(Ocoords,Ncoords,boundary)
    uON=normalize(vecON)
    rNO=peuclidean(Ncoords,Ocoords,boundary.side_lengths)
    mN=s.atoms[1].mass
    mO=s.atoms[2].mass
    cosθ=getcosth(Ncoords,Ocoords,boundary)
    dz=getzcom(mN,mO,Ncoords,Ocoords)

    F=getFij_NOAu(i,j,vec_ij,distbtwn,rNO,uON,u,cosθ,dz,a,b,c,boundary)
    F .|> sysunits
end

############################################################################################################

# mod molly simulate! to fix force fn lagging behind eigenvalues
function Molly.simulate!(sys::System{D, false, T, CU, A, AD, PI} where {D,T,CU,A,AD,PI<:Tuple{NOAuInteraction}},
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