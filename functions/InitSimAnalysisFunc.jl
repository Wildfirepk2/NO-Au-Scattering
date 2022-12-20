# functions to analyze Molly simulations after running

############################################################################################################

# extend length method for XLSXFile to give number of sheets. may revert later.
Base.length(xf::XLSX.XLSXFile) = XLSX.sheetcount(xf)

"""
write data to excel file. intermediate function for output functions
"""
function write_xlsx(file, i, stepslog, df::DataFrame)
    # step no of simulation
    step=(i-1)*stepslog

    # data to go in excel file
    data = collect(eachcol(df))
    cols = names(df)

    # if excel file not in directory, create it
    if !isfile(file)
        # xf: excel file object
        XLSX.openxlsx(file, mode="w") do xf
            sheet = xf[length(xf)]
            XLSX.writetable!(sheet, data, cols)
            XLSX.rename!(sheet, "step $step")
        end
    else
        # add sheet to existing file and write to it
        XLSX.openxlsx(file, mode="rw") do xf
            XLSX.addsheet!(xf)
            sheet = xf[length(xf)]
            XLSX.writetable!(sheet, data, cols)
            XLSX.rename!(sheet, "step $step")
        end
    end
end

############################################################################################################

"""
simple version of write_xlsx for outputing energies. called in outputsysE
"""
function write_xlsx(file, df::DataFrame)
    dfnounits=ustrip.(df)
    data = collect(eachcol(dfnounits))
    cols = names(dfnounits)
    XLSX.writetable(file, data, cols)
end

############################################################################################################

"""
output atom coords at ith logging point to excel file. default path is current directory. used only in outputallsyscoords but can be called independently
"""
function outputsyscoords(sys,i,path=".")
    # logging interval. ie log coords after every $stepslog steps
    stepslog=sys.loggers.coords.n_steps

    # coords from molly's loggers
    datasrc=sys.loggers.coords.history[i]

    # xyz positions stored as temp vars. converted to MD units then unit stripped
    xs=ustrip([datasrc[j][1] for j in eachindex(datasrc)] .|> u"d_MD")
    ys=ustrip([datasrc[j][2] for j in eachindex(datasrc)] .|> u"d_MD")
    zs=ustrip([datasrc[j][3] for j in eachindex(datasrc)] .|> u"d_MD")

    # write to xlsx
    data=DataFrame(x=xs,y=ys,z=zs)
    file="$path/syscoords.xlsx"
    write_xlsx(file,i,stepslog,data)
end

############################################################################################################

"""
output atom coords w time to excel file at path. default path is current directory. used only in outputsysinfo but can be called independently

firstlastonly for speed. XLSX slow.
"""
function outputallsyscoords(sys,path=".",firstlastonly=true)
    # coords from molly's loggers
    datasrc=sys.loggers.coords.history

    if firstlastonly
        # write excel file for first/last time step only.
        outputsyscoords(sys,1,path)
        outputsyscoords(sys,length(datasrc),path)
    else
        # write excel file for ea time step
        for i in eachindex(datasrc)
            outputsyscoords(sys,i,path)
        end
    end
end

############################################################################################################
"""
output graph given DataFrame. called in outputsysE
"""
function outputgraph(desc,df,path=".")
    # grabbing key data from df
    colnames=names(df)
    numcols=ncol(df)
    dffirstrow=first.(eachcol(df))
    colunittypes=qtytolabel.(dffirstrow)
    colunits=unit.(dffirstrow)
    dfnounits=ustrip.(df)

    # create plot obj
    xunittype=colunittypes[begin]
    yunittype=colunittypes[begin+1]
    xunit=colunits[begin]
    yunit=colunits[begin+1]
    fig=scatterlines(dfnounits[!,begin],dfnounits[!,begin+1],label=colnames[begin+1];
        axis = (; title = desc, xlabel = "$xunittype ($xunit)", ylabel = "$yunittype ($yunit)"))

    # add to plot obj
    for i in 2:numcols-1
        scatterlines!(dfnounits[!,begin],dfnounits[!,i+1],label=colnames[i+1])
    end
    axislegend()

    # save plot to $path
    save("$path/Evt.png", fig)
end

############################################################################################################

"""
output system energies w time to excel file + make graph
"""
function outputsysE(sys,systype,path=".")
    # number steps in run
    nsteps=length(sys.loggers.et.history)-1

    # time step of run (dt)
    dt=simulator.dt

    # tmp vars for time, kinetic, potential, and total energies + converted to MD units
    time=[i*dt for i in 0:nsteps] .|> u"t_MD"
    KE=sys.loggers.ke.history .|> u"e_MD"
    PE=sys.loggers.pe.history .|> u"e_MD"
    TE=sys.loggers.et.history .|> u"e_MD"
    
    # write to excel file at $path
    data=DataFrame(t=time,KE=KE,PE=PE,TE=TE)
    file="$path/sysE.xlsx"
    write_xlsx(file,data)

    # graph description based on system type
    graphdesc = systype=="Au" ? "Au slab equilibration: Checking energy conservation" : "NO/Au Scattering: Checking energy conservation"

    # make energy vs time graph
    outputgraph(graphdesc,data,path)
end

############################################################################################################

"""
output atom forces at ith logging point to excel file. default path is current directory. used only in outputallsysforces but can be called independently
"""
function outputsysforces(sys,i,path=".")
    # logging interval. ie log forces after every $stepslog steps
    stepslog=sys.loggers.forces.n_steps

    # forces from molly's loggers
    datasrc=sys.loggers.forces.history[i]

    # xyz forces stored as temp vars. converted to MD units then unit stripped
    fx=ustrip([datasrc[j][1] for j in eachindex(datasrc)] .|> u"e_MD/d_MD")
    fy=ustrip([datasrc[j][2] for j in eachindex(datasrc)] .|> u"e_MD/d_MD")
    fz=ustrip([datasrc[j][3] for j in eachindex(datasrc)] .|> u"e_MD/d_MD")

    # write to excel file
    data=DataFrame(Fx=fx,Fy=fy,Fz=fz)
    file="$path/sysforces.xlsx"
    write_xlsx(file,i,stepslog,data)
end

############################################################################################################
"""
output atom forces w time to excel file at path. default path is current directory. used only in outputsysinfo but can be called independently
"""
function outputallsysforces(sys,path=".")
    # forces from molly's loggers
    datasrc=sys.loggers.forces.history

    # write excel files for ea time step. store in forces folder
    for i in eachindex(datasrc)
        outputsysforces(sys,i,path)
    end
end

############################################################################################################
"""
output atom velocities at ith logging point to excel file. default path is current directory. used only in outputallsysvelocities but can be called independently
"""
function outputsysvelocities(sys,i,path=".")
    # logging interval. ie log velocities after every $stepslog steps
    stepslog=sys.loggers.velocities.n_steps

    # velocities from molly's loggers
    datasrc=sys.loggers.velocities.history[i]

    # xyz velocities stored as temp vars. converted to MD units then unit stripped
    vx=ustrip([datasrc[j][1] for j in eachindex(datasrc)] .|> u"v_MD")
    vy=ustrip([datasrc[j][2] for j in eachindex(datasrc)] .|> u"v_MD")
    vz=ustrip([datasrc[j][3] for j in eachindex(datasrc)] .|> u"v_MD")

    # write to excel file
    data=DataFrame(Vx=vx,Vy=vy,Vz=vz)
    file="$path/sysvelocities.xlsx"
    write_xlsx(file,i,stepslog,data)
end

############################################################################################################
"""
output atom velocities w time to excel file at path. default path is current directory. used only in outputsysinfo but can be called independently
"""
function outputallsysvelocities(sys,path=".")
    # velocities from molly's loggers
    datasrc=sys.loggers.velocities.history

    # write excel files for ea time step. store in velocities folder
    for i in eachindex(datasrc)
        outputsysvelocities(sys,i,path)
    end
end

############################################################################################################

"""
print txt file with summary of run.
"""
function outputsummary(sys,simsteps=NaN,runtime=NaN,runpath=".")
    # number of atoms in system
    natoms=length(sys.atoms)

    # time step of run (dt) in given/MD units. rounding due to floating point errors
    dt=simulator.dt
    dtmd=round(u"t_MD",dt;digits=2)

    # total time of run in given/MD units. rounding due to floating point errors
    ttotal=simsteps*dt
    ttotalmd=round(u"t_MD",ttotal;digits=2)

    # number of steps between logging points for quantities
    stepsbtwnlogsCoords=sys.loggers.coords.n_steps
    stepsbtwnlogsE=sys.loggers.et.n_steps
    stepsbtwnlogsF=sys.loggers.forces.n_steps
    stepsbtwnlogsV=sys.loggers.velocities.n_steps

    # time step between logging points for quantities in given/MD units. rounding due to floating point errors
    logdtCoords=stepsbtwnlogsCoords*dt
    logdtE=stepsbtwnlogsE*dt
    logdtF=stepsbtwnlogsF*dt
    logdtV=stepsbtwnlogsV*dt
    logdtCoordsmd=round(u"t_MD",logdtCoords;digits=2)
    logdtEmd=round(u"t_MD",logdtE;digits=2)
    logdtFmd=round(u"t_MD",logdtF;digits=2)
    logdtVmd=round(u"t_MD",logdtV;digits=2)

    # write to txt
    file="$runpath/summary.txt"
    open(file,"w") do io
        println(io,"Number of atoms: $natoms")
        println(io,"Number of steps: $simsteps")
        println(io,"Time step: $dt ($dtmd)")
        println(io,"Total time: $ttotal ($ttotalmd)")
        println(io,"Simulation runtime: $runtime")
        println(io)
        println(io,"Time period between logging quantities")
        println(io,"    Coords: $logdtCoords ($logdtCoordsmd)")
        println(io,"    Energies: $logdtE ($logdtEmd)")
        println(io,"    Forces: $logdtF ($logdtFmd)")
        println(io,"    Velocities: $logdtV ($logdtVmd)")
    end
end

############################################################################################################

"""
output all system info. animation. logger quantities

need to give description of run as string and specify type of system ("Au" or "NO/Au"). 
"""
function outputsysinfo(sys,systype,path=".")
    # output animation of simulation in $path
    visualize(sys.loggers.coords, sys.boundary, "$path/animation.mp4";show_boundary=false)

    # output quantities to excel file in separate folder
    outputallsyscoords(sys,path)
    outputsysE(sys,systype,path)
    outputallsysforces(sys,path)
    outputallsysvelocities(sys,path)
end

############################################################################################################

"""
run Au slab equilibration and output run info to results folder
"""
function AuSlabEquilibration(sys,simulator,steps)
    # make folder to store results of current run
    date=Dates.format(now(), "yyyy-mm-ddTHHMMSS")
    path=mkpath("results/$rundesc--$date")

    # run MD+give run time
    runtimeAu=@elapsed simulate!(sys, simulator, steps)
    runtimeAu*=u"s"

    # output all system data: animation, coords, last velocities/forces
    outputsysinfo(sys,"Au",path)

    # output summary of run
    outputsummary(sys,steps,runtimeAu,path)
end
