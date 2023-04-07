# functions to analyze Molly Au simulations after running

############################################################################################################

function getatomilabel(sys::System{D, false, T, CU, A, AD, PI} where {D,T,CU,A,AD,PI<:Tuple{AuSlabInteraction}},atom_i::Int64)
    "zAu$atom_i"
end

############################################################################################################

function getEtitle(sys::System{D, false, T, CU, A, AD, PI} where {D,T,CU,A,AD,PI<:Tuple{AuSlabInteraction}})
    "Au slab equilibration: Checking energy conservation"
end

############################################################################################################

function getsystitle(s::System{D, false, T, CU, A, AD, PI} where {D,T,CU,A,AD,PI<:Tuple{AuSlabInteraction}})
    "Au(111) Slab Equilibration"
end

############################################################################################################

function outputsummary(sys::System{D, false, T, CU, A, AD, PI} where {D,T,CU,A,AD,PI<:Tuple{AuSlabInteraction}},dt,T,simsteps=NaN,runtime=NaN,runpath=".")
    # description of run
    title="Au(111) Slab Equilibration"
    
    # time of run
    daterun=Dates.format(now(), "yyyy-mm-dd HH:MM:SS")

    # number of atoms in system
    natoms=length(sys.atoms)

    # time step of run (dt) in MD units. rounding due to floating point errors
    dtmd=round(u"t_MD",dt;digits=2)

    # total time of run in given/MD units. rounding due to floating point errors
    ttotal=simsteps*dt
    ttotalmd=round(u"t_MD",ttotal;digits=2)

    # number of steps between logging points for quantities
    stepsbtwnlogsCoords=sys.loggers.coords.n_steps
    stepsbtwnlogsE=sys.loggers.et.n_steps
    stepsbtwnlogsV=sys.loggers.velocities.n_steps

    # time step between logging points for quantities in given/MD units. rounding due to floating point errors
    logdtCoords=stepsbtwnlogsCoords*dt
    logdtE=stepsbtwnlogsE*dt
    logdtV=stepsbtwnlogsV*dt
    logdtCoordsmd=round(u"t_MD",logdtCoords;digits=2)
    logdtEmd=round(u"t_MD",logdtE;digits=2)
    logdtVmd=round(u"t_MD",logdtV;digits=2)

    # stepsbtwnlogsF=sys.loggers.forces.n_steps
    #     logdtF=stepsbtwnlogsF*dt
    #     logdtFmd=round(u"t_MD",logdtF;digits=2)

    if !simplerun
        stepsbtwnlogsF=sys.loggers.forces.n_steps
        logdtF=stepsbtwnlogsF*dt
        logdtFmd=round(u"t_MD",logdtF;digits=2)
    else
        stepsbtwnlogsF=NaN
        logdtF=NaN
        logdtFmd=NaN
    end

    # final slab energy in given/MD units. rounding due to floating point errors
    finalE=sys.loggers.et.history[end]
    finalEmd=round(u"e_MD",finalE;digits=2)

    # write to txt
    file="$runpath/summary.txt"
    open(file,"w") do io
        println(io,title)
        println(io)
        println(io,"Time of run: $daterun")
        println(io,"Number of atoms: $natoms")
        println(io,"Number of steps: $simsteps")
        println(io,"Time step: $dt ($dtmd)")
        println(io,"Total time: $ttotal ($ttotalmd)")
        println(io,"Simulation runtime: $runtime")
        println(io,"Ran on ISAAC: $isaac")
        println(io,"T: $T")
        println(io)
        println(io,"Steps between logging quantities")
        println(io,"    Coords: $stepsbtwnlogsCoords")
        println(io,"    Energies: $stepsbtwnlogsE")
        println(io,"    Forces: $stepsbtwnlogsF")
        println(io,"    Velocities: $stepsbtwnlogsV")
        println(io)
        println(io,"Time period between logging quantities")
        println(io,"    Coords: $logdtCoords ($logdtCoordsmd)")
        println(io,"    Energies: $logdtE ($logdtEmd)")
        println(io,"    Forces: $logdtF ($logdtFmd)")
        println(io,"    Velocities: $logdtV ($logdtVmd)")
        println(io)
        println(io,"Final System Total Energy: $finalE ($finalEmd)")
    end
end

############################################################################################################

function outputanimation(sys::System{D, false, T, CU, A, AD, PI} where {D,T,CU,A,AD,PI<:Tuple{AuSlabInteraction}},path=".")
    # output animation of simulation in $path
#\fix z label not showing
    kwargs = [:show_boundary => false,
            :color => :gold,]
    if isaac
        myvisualize(sys.loggers.coords, sys.boundary, "$path/animation.mp4";kwargs...)
    else
        visualize(sys.loggers.coords, sys.boundary, "$path/animation.mp4";kwargs...)
    end
end

############################################################################################################

function outputsysinfo(sys::System{D, false, T, CU, A, AD, PI} where {D,T,CU,A,AD,PI<:Tuple{AuSlabInteraction}},dt,path::String=".")
    # output quantities to excel file in separate folder
    if isaac
        outputallsyscoords(sys,dt,path,false,false)
    else
        outputallsyscoords(sys,dt,path)
    end
    outputallsysvelocities(sys,dt,path)
    if !simplerun
        outputanimation(sys,path)
        outputsysE(sys,dt,path)
        outputallsysforces(sys,dt,path)
    end
end
