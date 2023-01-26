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
	nntuples=[(i,nn[i][j],false) for i in eachindex(nn) for j in eachindex(nn[i])]

    # molly neighbor object
	NeighborList(length(nntuples), nntuples)

    # above looks ok

    # parallelization
    # dist_unit = unit(first(first(s.coords)))
#     bv = ustrip.(dist_unit, s.boundary)
#     btree = BallTree(ustrip_vec.(s.coords), PeriodicEuclidean(bv))
#     dist_cutoff = ustrip(dist_unit, nf.dist_cutoff)

#     @floop ThreadedEx(basesize = length(s) ÷ n_threads) for i in 1:length(s)
#         ci = ustrip.(s.coords[i])
#         nbi = @view nf.nb_matrix[:, i]
#         w14i = @view nf.matrix_14[:, i]
#         idxs = inrange(btree, ci, dist_cutoff, true)
#         for j in idxs
#             if nbi[j] && i > j
#                 nn = (i, j, w14i[j])
#                 @reduce(neighbors_list = append!(Tuple{Int, Int, Bool}[], (nn,)))
#             end
#         end
#     end

#     return NeighborList(length(neighbors_list), neighbors_list)
end

############################################################################################################

# alter molly PE function above for NO/Au. need to calculate E all at once
function Molly.potential_energy(s::System{D, false, T, CU, A, AD, PI} where {D,T,CU,A,AD,PI<:Tuple{NOAuInteraction}}, neighbors=nothing)
    En=0 * s.energy_units
    Ei=0 * s.energy_units
    Ec=0 * s.energy_units

    @inbounds for ni in 1:neighbors.n
        # setup for E calc
        i, j, weight_14 = neighbors.list[ni]
        dr = vector(s.coords[i], s.coords[j], s.boundary)
        distbtwn=euclidean(dr,zeros(3)u"Å")

        Ncoords=s.coords[1]
        Ocoords=s.coords[2]
        mN=s.atoms[1].mass
        mO=s.atoms[2].mass
        cosθ=getcosth(Ncoords,Ocoords)
        dz=getzcom(mN,mO,Ncoords,Ocoords)

        t_En,t_Ei,t_Ec=getVij_NOAu(i,j,distbtwn,cosθ,dz)
        En+=t_En
        Ei+=t_Ei
        Ec+=t_Ec
    end

    # final ground state energy/eigenvalues
    Eg=1/2*(En+Ei-√((En-Ei)^2+4Ec^2))
    λ1=(Ei-Eg)/√((Ei-Eg)^2+Ec^2)
    λ2=-Ec/√((Ei-Eg)^2+Ec^2)

    # store eigenvalues for force calculation
    push!(storeEs,[Eg,λ1,λ2])
    
    return uconvert(s.energy_units, Eg)
end

############################################################################################################

function getVij_NOAu(i,j,distbtwn,cosθ,dz)
    En=0u"e_MD"
    Ei=0u"e_MD"
    Ec=0u"e_MD"
    if i==1
        # NO interaction + 1 time V_image, WF, Ea calc
        if j==2
            if neutral_PES_active
                En=V00_NO(distbtwn)
            end
            if ionic_PES_active
                Ei=V11_NO(distbtwn)+V11_image(dz)+PES_ionic.ϕ[1]-PES_ionic.Ea[1]
            end
        # N-Au interaction
        else
            if neutral_PES_active
                En=V00_AuN(distbtwn)
            end
            if ionic_PES_active
                Ei=V11_AuN(distbtwn,cosθ)
            end
            if coupled_PES_active
                Ec=V01_AuN(distbtwn)
            end
        end
    else
        # O-Au interactions
        if j>1
            if neutral_PES_active
                En=V00_AuO(distbtwn)
            end
            if ionic_PES_active
                Ei=V11_AuO(distbtwn)
            end
            if coupled_PES_active
                Ec=V01_AuO(distbtwn)
            end
        end
    end
    return En,Ei,Ec
end

############################################################################################################

# force function for NO/Au scattering. calculating F for each atom due to NN. 
# inline/inbounds: copied from Lennard Jones force function. may also consider @fastmath
# see fortran FORCE_MATRIX
@inline @inbounds function Molly.force(inter::NOAuInteraction,
                            vec_ij,
                            coord_i,
                            coord_j,
                            atom_i,
                            atom_j,
                            boundary)
    # fetch atom indices
    i=atom_i.index
    j=atom_j.index
    sysunits=sys_NOAu.force_units

    # get eigenvalues
    λ1=storeEs[end,2]
    λ2=storeEs[end,3]
    a=λ1^2
    b=λ2^2
    c=2*λ1*λ2

    # setup for force calculation
    distbtwn=euclidean(vec_ij,zeros(3)u"Å")
    u=normalize(vec_ij)

    Ncoords=sys_NOAu.coords[1]
    Ocoords=sys_NOAu.coords[2]
    mN=sys_NOAu.atoms[1].mass
    mO=sys_NOAu.atoms[2].mass
    cosθ=getcosth(Ncoords,Ocoords)
    dz=getzcom(mN,mO,Ncoords,Ocoords)

    F=getFij_NOAu(i,j,distbtwn,cosθ,dz,a,b,c)
    F*u .|> sysunits
end

############################################################################################################

function getFij_NOAu(i,j,distbtwn,cosθ,dz,a,b,c)
    Fn=0u"N/mol"
    Fi=0u"N/mol"
    Fc=0u"N/mol"

    if i==1
        # NO interaction
        if j==2
            if neutral_PES_active
                Fn=F00_NO(distbtwn)
            end
            if ionic_PES_active
                Fi=F11_NO(distbtwn)
            end
            if coupled_PES_active
                Fc=F11_image(dz)
            end
        # N-Au interaction
        else
            if neutral_PES_active
                Fn=F00_AuN(distbtwn)
            end
            if ionic_PES_active
                Fi=F11_AuN(distbtwn,cosθ)
            end
            if coupled_PES_active
                Fc=F01_AuN(distbtwn)
            end
        end
    else
        # O-Au interactions 
        if j>1
            if neutral_PES_active
                Fn=F00_AuO(distbtwn)
            end
            if ionic_PES_active
                Fi=F11_AuO(distbtwn)
            end
            if coupled_PES_active
                Fc=F01_AuO(distbtwn)
            end
        end
    end
    Fn*a+Fi*b+Fc*c
end
