# functions to analyze Molly simulations after running

############################################################################################################

"""
output atom coords at ith logging point to csv. default path is current directory. used only in outputallsyscoords but can be called independently
"""
function outputsyscoords(sys,i,path=".")
    # logging interval. ie log coords after every $stepslog steps
    stepslog=sys.loggers.coords.n_steps

    # step no of simulation
    step=(i-1)*stepslog

    # coords from molly's loggers
    datasrc=sys.loggers.coords.history[i]

    # xyz positions stored as temp vars. converted to MD units then unit stripped
    xs=ustrip([datasrc[j][1] for j in eachindex(datasrc)] .|> u"d_MD")
    ys=ustrip([datasrc[j][2] for j in eachindex(datasrc)] .|> u"d_MD")
    zs=ustrip([datasrc[j][3] for j in eachindex(datasrc)] .|> u"d_MD")

    # write to csv
    data=DataFrame(x=xs,y=ys,z=zs)
    file="$path/syscoords-step $step.csv"
    CSV.write(file,data)
end

############################################################################################################

"""
output atom coords w time to csv at path. default path is current directory. used only in outputsysinfo but can be called independently
"""
function outputallsyscoords(sys,path=".")
    # create folder for coords at $path
    coordpath=mkpath("$path/coords")

    # coords from molly's loggers
    datasrc=sys.loggers.coords.history

    # write csv files for ea time step. store in coords folder
    for i in eachindex(datasrc)
        outputsyscoords(sys,i,coordpath)
    end
end

############################################################################################################

"""
output system energies w time to csv + make graph
"""
function outputsysE(sys,systype,path=".")
    # create folder for energies at $path
    Epath=mkpath("$path/energies")

    # number steps in run
    nsteps=length(sys.loggers.et.history)-1

    # time step of run (dt)
    dt=simulator.dt

    # tmp vars for time, kinetic, potential, and total energies. converted to MD units then unit stripped
    time=ustrip([i*dt for i in 0:nsteps] .|> u"t_MD")
    KE=ustrip(sys.loggers.ke.history .|> u"e_MD")
    PE=ustrip(sys.loggers.pe.history .|> u"e_MD")
    TE=ustrip(sys.loggers.et.history .|> u"e_MD")
    
    # write to csv at $Epath
    data=DataFrame(t=time,KE=KE,PE=PE,TE=TE)
    file="$Epath/sysE.csv"
    CSV.write(file,data)

    ## make E vs t graph. possibly make more generic in future?

    # graph description based on system type
    graphdesc = systype=="Au" ? "Au slab equilibration: Checking energy conservation" : "NO/Au Scattering: Checking energy conservation"

    # create plot obj w params. add KE series. 
    fig,ax,KEseries=scatterlines(time,KE,label="KE";
    axis = (; title = graphdesc, xlabel = "Time (t_MD)", ylabel = "Energy (e_MD)"))

    # add PE/TE series
    PEseries=scatterlines!(time,PE,label="PE")
    TEseries=scatterlines!(time,TE,label="TE")
    axislegend()
    
    # save plot to $Epath
    save("$Epath/Evt.png", fig)
end

############################################################################################################

"""
output forces on atoms at ith logging point to csv. default path is current directory. used only in outputsysinfo but can be called independently
"""
function outputsysforces(sys,i,path=".")
    # logging interval. ie log forces after every $stepslog steps
    stepslog=sys.loggers.forces.n_steps

    # step no of simulation
    step=(i-1)*stepslog

    # forces from molly's loggers
    datasrc=sys.loggers.forces.history[i]

    # xyz forces stored as temp vars. converted to MD units then unit stripped
    fx=ustrip([datasrc[j][1] for j in eachindex(datasrc)] .|> u"e_MD/d_MD")
    fy=ustrip([datasrc[j][2] for j in eachindex(datasrc)] .|> u"e_MD/d_MD")
    fz=ustrip([datasrc[j][3] for j in eachindex(datasrc)] .|> u"e_MD/d_MD")

    # write to csv
    data=DataFrame(Fx=fx,Fy=fy,Fz=fz)
    file="$path/sysforces-step $step.csv"
    CSV.write(file,data)
end

############################################################################################################

"""
output velocities on atoms at ith logging point to csv. default path is current directory. used only in outputsysinfo but can be called independently
"""
function outputsysvelocities(sys,i,path=".")
    # logging interval. ie log velocities after every $stepslog steps
    stepslog=sys.loggers.velocities.n_steps

    # step no of simulation
    step=(i-1)*stepslog

    # velocities from molly's loggers
    datasrc=sys.loggers.velocities.history[i]

    # xyz velocities stored as temp vars. converted to MD units then unit stripped
    vx=ustrip([datasrc[j][1] for j in eachindex(datasrc)] .|> u"v_MD")
    vy=ustrip([datasrc[j][2] for j in eachindex(datasrc)] .|> u"v_MD")
    vz=ustrip([datasrc[j][3] for j in eachindex(datasrc)] .|> u"v_MD")

    # write to csv
    data=DataFrame(vx=vx,vy=vy,vz=vz)
    file="$path/sysvelocities-step $step.csv"
    CSV.write(file,data)
end

############################################################################################################

"""
print txt file with summary of run.
"""
function outputsummary(sys,runpath)
    # create folder for summary at $path
    summarypath=mkpath("$runpath/summary")

    # number steps in run
    nsteps=sys.loggers.forces.n_steps

    # time step of run (dt)
    dt=simulator.dt

    # number of atoms in system
    natoms=length(sys.atoms)

    # number of atoms in system
    nsteps=length(sys.loggers.forces.history)

    # time step of run (dt)
    dt=simulator.dt

    # total time of run
    ttotal=nsteps*dt

    # write to txt
    file="$summarypath/summary.txt"
    open(file,"w") do io
        println(io,"Number of atoms: $natoms")
        println(io,"Number of steps: $nsteps")
        println(io,"Time step: $dt")
        println(io,"Total time: $ttotal")
    end
end

############################################################################################################

"""
output all system info. animation. logger stuff. forces/velocities last step

need to give description of run as string and specify type of system ("Au" or "NO/Au"). 
"""
function outputsysinfo(sys,rundesc,systype)
    # make folder to store results of current run
    date=Dates.format(now(), "yyyy-mm-ddTHHMMSS")
    mainpath=mkpath("results/$rundesc--$date")

    # output animation of simulation in $mainpath
    visualize(sys.loggers.coords, sys.boundary, "$mainpath/animation.mp4";show_boundary=false)

    # output coords/energies to csv in separate folder
    outputallsyscoords(sys,mainpath)
    outputsysE(sys,systype,mainpath)

    # output final forces/velocities to csv in separate folder
    laststeppath=mkpath("$mainpath/last step")
    i_lastforce=length(sys.loggers.forces.history)
    i_lastvel=length(sys.loggers.velocities.history)
    outputsysforces(sys,i_lastforce,laststeppath)
    outputsysvelocities(sys,i_lastvel,laststeppath)
end
