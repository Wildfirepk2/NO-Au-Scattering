# custom Molly interactions/force functions

############################################################################################################

using Molly
using GLMakie

############################################################################################################

# interaction for Au slab equilibration
struct AuSlabInteraction <: PairwiseInteraction
    nl_only::Bool
end

############################################################################################################

# force function for Au slab equilibration. inline/inbounds: copied from Lennard Jones force function. may also consider @fastmath
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
        step=(i-1)*stepslogging

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

    # output final coords to csv in separate folder
    outputsyscoords(sys,mainpath)

    # output final velocities to csv
    data=DataFrame(vx=sys.loggers.velocity.history[end][:,1],vy=sys.loggers.velocity.history[end][:,2],vz=sys.loggers.velocity.history[end][:,3])
    CSV.write("results/au slab equilibration/velocities/velocities-last step.csv",data)

    # output final forces to csv
    data=DataFrame(Fx=sys.loggers.forces.history[end][:,1],Fy=sys.loggers.forces.history[end][:,2],Fz=sys.loggers.forces.history[end][:,3])
    CSV.write("results/au slab equilibration/forces/forces-last step.csv",data)
end