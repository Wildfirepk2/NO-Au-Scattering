# custom Molly interactions/force functions

############################################################################################################

using Molly
using GLMakie
using Dates

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
        drij=[nncoords[j]-coord_i for j in eachindex(nncoords)]

        # dist of nn pairs away from equilibrium. array of arrays
        rij=[drij[j]-r0ij[i][:,j] for j in eachindex(drij)]

        # force each nn pair exerts on atom i. array of arrays
        Fij=[-Aijarray[i][j]*rij[j] for j in eachindex(rij)]

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
        zeros(3)sysunits
    end
end

############################################################################################################

"""
output atom coords w time to csv
"""
function outputsyscoords(sys,path)
    # create folder for coords at $path
    coordpath=mkpath("$path/coords")

    # coords from molly's loggers
    datasrc=sys.loggers.coords.history

    # write csv files for ea time step
    for i in eachindex(datasrc)
        # simulation step no
        step=(i-1)*actsteplog

        # xyz positions stored as temp vars
        xs=[datasrc[i][j][1] for j in eachindex(datasrc[i])]
        ys=[datasrc[i][j][2] for j in eachindex(datasrc[i])]
        zs=[datasrc[i][j][3] for j in eachindex(datasrc[i])]

        # write to csv
        data=DataFrame(x=xs,y=ys,z=zs)
        file="$coordpath/syscoords-step $step.csv"
        CSV.write(file,data)
    end
end

############################################################################################################

"""
output system energies w time to csv
"""
function outputsysE(sys,path)
    # create folder for energies at $path
    Epath=mkpath("$path/energies")

    # number steps in run
    nsteps=length(sys.loggers.et.history)-1

    # time step of run (dt)
    dt=simulator.dt

    # tmp vars for time, kinetic, potential, and total energies
    time=[i*dt for i in 0:nsteps]
    KE=sys.loggers.ke.history
    PE=sys.loggers.pe.history
    TE=sys.loggers.et.history

    # write to csv
    data=DataFrame(t=time,KE=KE,PE=PE,TE=TE)
    file="$Epath/sysE.csv"
    CSV.write(file,data)
end

############################################################################################################

"""
output last forces on atoms to csv. called only in outputsysinfo
"""
function outputlastforces(sys,path)
    # create folder for last forces at $path
    forcepath=mkpath("$path/last step/forces")

    # number steps in run
    nsteps=sys.loggers.forces.n_steps

    # time step of run (dt)
    dt=simulator.dt

    # coords from molly's loggers
    datasrc=sys.loggers.forces.history

    # Fx,Fy,Fz stored as temp vars
    xs=[datasrc[end][j][1] for j in eachindex(datasrc[end])]
    ys=[datasrc[end][j][2] for j in eachindex(datasrc[end])]
    zs=[datasrc[end][j][3] for j in eachindex(datasrc[end])]

    # write to csv
    data=DataFrame(x=xs,y=ys,z=zs)
    file="$forcepath/forces-step $nsteps-dt $dt.csv"
    CSV.write(file,data)
end

############################################################################################################

"""
output last velocities on atoms to csv. called only in outputsysinfo
"""
function outputlastvelocities(sys,path)
    # create folder for last velocities at $path
    velocitypath=mkpath("$path/last step/velocities")

    # number steps in run
    nsteps=sys.loggers.velocities.n_steps

    # time step of run (dt)
    dt=simulator.dt

    # coords from molly's loggers
    datasrc=sys.loggers.velocities.history

    # Vx,Vy,Vz stored as temp vars
    xs=[datasrc[end][j][1] for j in eachindex(datasrc[end])]
    ys=[datasrc[end][j][2] for j in eachindex(datasrc[end])]
    zs=[datasrc[end][j][3] for j in eachindex(datasrc[end])]

    # write to csv
    data=DataFrame(x=xs,y=ys,z=zs)
    file="$velocitypath/velocities-step $nsteps-dt $dt.csv"
    CSV.write(file,data)
end

############################################################################################################

"""
output all system info. animation. logger stuff. forces/velocities last step

need to give a system description as string. 
"""
function outputsysinfo(sys,desc)
    # make directories to store results of current run
    date=Dates.format(now(), "yyyy-mm-ddTHMS")
    mainpath="results/$desc--$date"
    mkpath(mainpath)

    # output animation of simulation in $mainpath
    visualize(sys.loggers.coords, sys.boundary, "$mainpath/animation.mp4")

    # output coords/energies to csv in separate folder
    outputsyscoords(sys,mainpath)
    outputsysE(sys,mainpath)

    # output final forces/velocities to csv in separate folder
    outputlastforces(sys,mainpath)
    outputlastvelocities(sys,mainpath)
end

############################################################################################################

# potential energy (V) for Au slab equilibration. calculating V for each atom. molly will sum all Vs. similar to force function. inline/inbounds: copied from Lennard Jones force function. may also consider @fastmath
function Molly.potential_energy(inter::AuSlabInteraction,
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
    drij=[nncoords[j]-coord_i for j in eachindex(nncoords)]

    # dist of nn pairs away from equilibrium. array of arrays
    rij=[drij[j]-r0ij[i][:,j] for j in eachindex(drij)]

    # force each nn pair exerts on atom i. array of arrays
    Fij=[-Aijarray[i][j]*rij[j] for j in eachindex(rij)]

    # V each nn pair exerts on atom i. array of numbers
    Vij=[dot(rij[j],-Fij[j]) for j in eachindex(Fij)]

    # total V on atom i. divide by 2 because of v definition
    V=sum(Vij)/2
end
