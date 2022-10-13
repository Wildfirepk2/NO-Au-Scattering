# using DataFrames
# using CSV
# using NearestNeighbors
# using LinearAlgebra
# using Molly
# using GLMakie



#=
n_atoms = au.N[1]
boundary = CubicBoundary(au.aPBCx[1]*u"Å", au.aPBCy[1]*u"Å", Inf*u"Å")
temp = param.T[1]*u"K"
atom_mass = mAu*u"u"
nsteps=param.Nsteps_eq[1]

atoms = [Atom(index=i, mass=atom_mass) for i in 1:n_atoms]
coords = [SA[au.x[i]*u"Å",au.y[i]*u"Å",-au.z[i]*u"Å"] for i in 1:n_atoms]
velocities = [velocity(atom_mass, temp) for i in 1:n_atoms]
pairwise_inters = ()
simulator = VelocityVerlet(
    dt=param.dt[1]*u"10fs",
    coupling=AndersenThermostat(temp, 1.0u"10fs"),
)

sys = System(
    atoms=atoms,
    pairwise_inters=pairwise_inters,
    coords=coords,
    velocities=velocities,
    boundary=boundary,
    loggers=(
        temp=TemperatureLogger(100),
        coords=CoordinateLogger(10),
    ),
)

simulate!(sys, simulator, nsteps)
visualize(sys.loggers.coords, boundary, "au no int_inf.mp4")


=#