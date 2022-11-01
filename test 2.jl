

using Molly
using GLMakie



###########################################################################################

# 8/20/22
# julia> @btime include("test.jl");
#   27.319 s (1473790 allocations: 413.36 MiB)

n_atoms2 = 100
boundary2 = CubicBoundary(2.0u"nm", 2.0u"nm", 2.0u"nm")
temp2 = 298.0u"K"
atom_mass2 = 10.0u"u"

atoms2 = [Atom(mass=atom_mass2, σ=0.3u"nm", ϵ=0.2u"kJ * mol^-1") for i in 1:n_atoms2]
coords2 = place_atoms(n_atoms2, boundary2; min_dist=0.3u"nm")
velocities2 = [velocity(atom_mass2, temp2) for i in 1:n_atoms2]
pairwise_inters2 = (LennardJones(),)
simulator2 = VelocityVerlet(
    dt=0.002u"ps",
    coupling=AndersenThermostat(temp2, 1.0u"ps"),
)

sys2 = System(
    atoms=atoms2,
    pairwise_inters=pairwise_inters2,
    coords=coords2,
    velocities=velocities2,
    boundary=boundary2,
    loggers=(
        temp=TemperatureLogger(100),
        coords=CoordinateLogger(10),
    ),
)

simulate!(sys2, simulator2, 10_000)
# visualize(sys2.loggers.coords, boundary2, "sim_lj2.mp4")

###########################################################################################

# julia> @btime include("test.jl");
#   20.482 s (1743160 allocations: 265.73 MiB)

@enum Status susceptible infected recovered

# Custom atom type
mutable struct Person
    i::Int
    status::Status
    mass::Float64
    σ::Float64
    ϵ::Float64
end

Molly.mass(person::Person) = person.mass

# Custom PairwiseInteraction
struct SIRInteraction <: PairwiseInteraction
    nl_only::Bool
    dist_infection::Float64
    prob_infection::Float64
    prob_recovery::Float64
end

# Custom force function
function Molly.force(inter::SIRInteraction,
                        vec_ij,
                        coord_i,
                        coord_j,
                        atom_i,
                        atom_j,
                        boundary)
    if (atom_i.status == infected && atom_j.status == susceptible) ||
                (atom_i.status == susceptible && atom_j.status == infected)
        # Infect close people randomly
        r2 = sum(abs2, vec_ij)
        if r2 < inter.dist_infection^2 && rand() < inter.prob_infection
            atom_i.status = infected
            atom_j.status = infected
        end
    end
    # Workaround to obtain a self-interaction
    if atom_i.i == (atom_j.i + 1)
        # Recover randomly
        if atom_i.status == infected && rand() < inter.prob_recovery
            atom_i.status = recovered
        end
    end
    return zero(coord_i)
end

# Custom Logger
function fracs_SIR(s::System, neighbors=nothing; n_threads::Integer=Threads.nthreads())
    counts_sir = [
        count(p -> p.status == susceptible, s.atoms),
        count(p -> p.status == infected   , s.atoms),
        count(p -> p.status == recovered  , s.atoms)
    ]
    return counts_sir ./ length(s)
end

SIRLogger(n_steps) = GeneralObservableLogger(fracs_SIR, Vector{Float64}, n_steps)

temp = 1.0
boundary = RectangularBoundary(10.0, 10.0)
n_steps = 1_000
n_people = 500
n_starting = 2
atoms = [Person(i, i <= n_starting ? infected : susceptible, 1.0, 0.1, 0.02) for i in 1:n_people]
coords = place_atoms(n_people, boundary; min_dist=0.1)
velocities = [velocity(1.0, temp; dims=2) for i in 1:n_people]
pairwise_inters = (
    LennardJones=LennardJones(nl_only=true),
    SIR=SIRInteraction(false, 0.5, 0.06, 0.01),
)
neighbor_finder = DistanceNeighborFinder(
    nb_matrix=trues(n_people, n_people),
    n_steps=10,
    dist_cutoff=2.0,
)
simulator = VelocityVerlet(
    dt=0.02,
    coupling=AndersenThermostat(temp, 5.0),
)

sys = System(
    atoms=atoms,
    pairwise_inters=pairwise_inters,
    coords=coords,
    velocities=velocities,
    boundary=boundary,
    neighbor_finder=neighbor_finder,
    loggers=(
        coords=CoordinateLogger(Float64, 10; dims=2),
        SIR=SIRLogger(10),
    ),
    force_units=NoUnits,
    energy_units=NoUnits,
)

simulate!(sys, simulator, n_steps)

# visualize(sys.loggers.coords, boundary, "sim_agent.mp4"; markersize=0.1)