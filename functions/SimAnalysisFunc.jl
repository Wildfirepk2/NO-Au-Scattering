# functions to analyze Molly simulations after running

############################################################################################################

# extend length method for XLSXFile to give number of sheets. may revert later.
Base.length(xf::XLSX.XLSXFile) = XLSX.sheetcount(xf)

"""
write data to excel file. intermediate function for output functions
"""
function write_xlsx(file, dt, i, stepslog, df::DataFrame)
    # step no of simulation
    step=(i-1)*stepslog

    # time at step no
    time=round(unit(dt),dt*step;digits=2)
    
    # data to go in excel file
    data = collect(eachcol(df))
    cols = names(df)

    # if excel file not in directory, create it
    if !isfile(file)
        # xf: excel file object
        XLSX.openxlsx(file, mode="w") do xf
            sheet = xf[length(xf)]
            XLSX.writetable!(sheet, data, cols)
            XLSX.rename!(sheet, "step $step (t=$time)")
        end
    else
        # add sheet to existing file and write to it
        XLSX.openxlsx(file, mode="rw") do xf
            XLSX.addsheet!(xf)
            sheet = xf[length(xf)]
            XLSX.writetable!(sheet, data, cols)
            XLSX.rename!(sheet, "step $step (t=$time)")
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
function outputsyscoords(sys,dt,i,path=".",asexcel=true)
    # logging interval. ie log coords after every $stepslog steps
    stepslog=sys.loggers.coords.n_steps

    # coords from molly's loggers
    datasrc=sys.loggers.coords.history[i]

    # xyz positions stored as temp vars. converted to MD units then unit stripped
    xs=ustrip([datasrc[j][1] for j in eachindex(datasrc)] .|> u"d_MD")
    ys=ustrip([datasrc[j][2] for j in eachindex(datasrc)] .|> u"d_MD")
    zs=ustrip([datasrc[j][3] for j in eachindex(datasrc)] .|> u"d_MD")

    if asexcel
        # write to excel
        data=DataFrame(x=xs,y=ys,z=zs)
        file="$path/syscoords.xlsx"
        write_xlsx(file,dt,i,stepslog,data)
    else
        # write to csv
        data=DataFrame(x=xs,y=ys,z=zs)
        file="$path/syscoords $i.csv"
        CSV.write(file,data)
    end
end

############################################################################################################

"""
output atom coords w time to excel file at path. default path is current directory. used only in outputsysinfo but can be called independently

firstlastonly for speed. XLSX slow.
"""
function outputallsyscoords(sys,dt,path=".",asexcel=true,firstlastonly=true)
    # coords from molly's loggers
    datasrc=sys.loggers.coords.history

    if firstlastonly
        if asexcel
            # write excel file for first/last time step only.
            outputsyscoords(sys,dt,1,path)
            outputsyscoords(sys,dt,length(datasrc),path)
        else
            coordpath=mkpath("$path/syscoords")
            # write csv file for first/last time step only.
            outputsyscoords(sys,dt,1,coordpath,false)
            outputsyscoords(sys,dt,length(datasrc),coordpath,false)
        end
    else
        if asexcel
            # write excel file for ea time step
            for i in eachindex(datasrc)
                outputsyscoords(sys,dt,i,path)
            end
        else
            coordpath=mkpath("$path/syscoords")
            # write csv file for ea time step
            for i in eachindex(datasrc)
                outputsyscoords(sys,dt,i,coordpath,false)
            end
        end
    end
end

############################################################################################################

"""
output atom i z coords w time to excel file + make graph
"""
function outputallatomizcoords(sys,dt,atom_i::Int64,path=".")
    # logging interval. ie log coords after every $stepslog steps
    stepslog=sys.loggers.coords.n_steps

    # number of logs
    nsteps=length(sys.loggers.coords.history)

    # time between logs
    dtlog=stepslog*dt

    # tmp vars for time/coords + converted to MD units
    z_surf=maximum(au.z)
    time=[(i-1)*dtlog for i in 1:nsteps] .|> u"t_MD" # i-1 because first log is at step 0
    zcoords=[sys.loggers.coords.history[i][atom_i][3]-z_surf for i in eachindex(sys.loggers.coords.history)] .|> u"d_MD"
    
    # write to excel file at $path
    data=DataFrame(t=time,z=zcoords)
    file="$path/z$atom_i-z_surf v time.xlsx"
    write_xlsx(file,data)
    
    # make zcoord vs time graph
    graphdesc = "Atom $atom_i Height Above Surface with Time"
    outputgraph(graphdesc,data,path)
end

"""
output atom z coords (atoms #'s specified by vec) w time to excel file + make graph
"""
function outputallatomizcoords(sys,dt,atomrange,path=".")
    # logging interval. ie log coords after every $stepslog steps
    stepslog=sys.loggers.coords.n_steps

    # number of logs
    nsteps=length(sys.loggers.coords.history)

    # time between logs
    dtlog=stepslog*dt

    # tmp vars for time/coords + converted to MD units
    z_surf=maximum(au.z)
    time=[(i-1)*dtlog for i in 1:nsteps] .|> u"t_MD" # i-1 because first log is at step 0
    
    dfheader=["t"]
    dfdata=Vector[time]
    
    for atom_i in atomrange
        push!(dfheader,"z$atom_i")
        zcoords=[sys.loggers.coords.history[i][atom_i][3]-z_surf for i in eachindex(sys.loggers.coords.history)] .|> u"d_MD"
        push!(dfdata,zcoords)
    end
    
    # write to excel file at $path
    data=DataFrame(Dict(zip(dfheader,dfdata)))
    file="$path/z-z_surf v time.xlsx"
    write_xlsx(file,data)
    
    # make zcoord vs time graph
    graphdesc = "Atom Height Above Surface with Time"
    outputgraph(graphdesc,data,path)
end

############################################################################################################

"""
output graph given DataFrame. called in outputsysE
"""
function outputgraph(desc,df::DataFrame,path=".")
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
    save("$path/$yunittype v $xunittype.png", fig)
end

############################################################################################################

"""
output system energies w time to excel file + make graph
"""
function outputsysE(sys,dt,systype,path=".")
    # logging interval. ie log energies after every $stepslog steps
    stepslog=sys.loggers.et.n_steps

    # number of logs
    nsteps=length(sys.loggers.et.history)

    # time between logs
    dtlog=stepslog*dt

    # tmp vars for time, kinetic, potential, and total energies + converted to MD units
    time=[(i-1)*dtlog for i in 1:nsteps] .|> u"t_MD" # i-1 because first log is at step 0
    KE=sys.loggers.ke.history .|> u"e_MD"
    PE=sys.loggers.pe.history .|> u"e_MD"
    TE=sys.loggers.et.history .|> u"e_MD"
    
    # write to excel file at $path
    data=DataFrame(t=time,KE=KE,PE=PE,TE=TE)
    file="$path/sysE.xlsx"
    write_xlsx(file,data)

    # graph description based on system type
    graphdesc = contains(systype,"Au slab") ? "Au slab equilibration: Checking energy conservation" : "NO/Au Scattering: Checking energy conservation"

    # make energy vs time graph
    outputgraph(graphdesc,data,path)
end

############################################################################################################

"""
output atom forces at ith logging point to excel file. default path is current directory. used only in outputallsysforces but can be called independently
"""
function outputsysforces(sys,dt,i,path=".")
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
    write_xlsx(file,dt,i,stepslog,data)
end

############################################################################################################

"""
output atom forces w time to excel file at path. default path is current directory. used only in outputsysinfo but can be called independently
"""
function outputallsysforces(sys,dt,path=".")
    # forces from molly's loggers
    datasrc=sys.loggers.forces.history

    # write excel files for ea time step. store in forces folder
    for i in eachindex(datasrc)
        outputsysforces(sys,dt,i,path)
    end
end

############################################################################################################

"""
output atom velocities at ith logging point to excel file. default path is current directory. used only in outputallsysvelocities but can be called independently
"""
function outputsysvelocities(sys,dt,i,path=".")
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
    write_xlsx(file,dt,i,stepslog,data)
end

############################################################################################################

"""
output atom velocities w time to excel file at path. default path is current directory. used only in outputsysinfo but can be called independently
"""
function outputallsysvelocities(sys,dt,path=".")
    # velocities from molly's loggers
    datasrc=sys.loggers.velocities.history

    # write excel files for ea time step. store in velocities folder
    for i in eachindex(datasrc)
        outputsysvelocities(sys,dt,i,path)
    end
end

############################################################################################################

"""
print txt file with summary of run.
"""
function outputsummary(sys,dt,desc,simsteps=NaN,runtime=NaN,runpath=".")
    # description of run
    title = desc==aurundesc ? "Au(111) Slab Equilibration" : "NO/Au(111) Scattering"
    
    # time of run
    daterun=Dates.format(now(), "yyyy-mm-dd HH:MM:SS")

    # number of atoms in system
    natoms=length(sys.atoms)

    # time step of run (dt) in MD units. rounding due to floating point errors
    dtmd=round(u"t_MD",dt;digits=2)

    # total time of run in given/MD units. rounding due to floating point errors
    ttotal=simsteps*dt
    ttotalmd=round(u"t_MD",ttotal;digits=2)

    # mult runs
    ei=no.Et_i[1]
    eimd=round(u"e_MD",ei;digits=2)

    # xy position of N
    xNi=sys.loggers.coords.history[1][1][1] |> u"Å"
    yNi=sys.loggers.coords.history[1][1][2] |> u"Å"

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

    if !multirun
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
        println(io)
        println(io,"T: $(param.T[1])")
        println(io,"Incident energy of molecule: $ei ($eimd)")
        println(io,"Initial x position of molecule: $xNi")
        println(io,"Initial y position of molecule: $yNi")
        println(io)
        if desc==noaurundesc
            println(io,"PESs:")
            println(io,"    neutral_PES_active = $neutral_PES_active")
            println(io,"    ionic_PES_active = $ionic_PES_active")
            println(io,"    coupled_PES_active = $coupled_PES_active")
            println(io)
        end
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

"""
output animation of trajectory for system
"""
function outputanimation(sys,path=".")
    # type of interaction. au or no/au
    inter_type=typeof(first(sys.pairwise_inters))
    
    if inter_type==AuSlabInteraction
        # output animation of simulation in $path
        if isaac
            myvisualize(sys.loggers.coords, sys.boundary, "$path/animation.mp4";show_boundary=false)
        else
            visualize(sys.loggers.coords, sys.boundary, "$path/animation.mp4";show_boundary=false)
        end
    elseif inter_type==NOAuInteraction
        # connecting N and O atoms with purple line
        connections=[(1,2)]
        connection_color=[:purple]

        # colors of all atoms
        NOcolors=[:blue,:red]
        Aucolors=[:gold for _ in 1:au.N[1]]
        syscolors=vcat(NOcolors,Aucolors)

        # output animation of simulation in $path
        if isaac
            myvisualize(sys.loggers.coords, sys.boundary, "$path/animation.mp4";show_boundary=false,connections=connections,connection_color=connection_color,color=syscolors,)
        else
            visualize(sys.loggers.coords, sys.boundary, "$path/animation.mp4";show_boundary=false,connections=connections,connection_color=connection_color,color=syscolors,)
        end
    end
end

############################################################################################################

"""
output all system info. animation. logger quantities

run description needs to contain either "Au slab" or "NO/Au"
"""
function outputsysinfo(sys,dt,systype,path=".",simplerun=false)
    # output quantities to excel file in separate folder
    if isaac
        outputallsyscoords(sys,dt,path,false,false)
    else
        outputallsyscoords(sys,dt,path)
    end
    if systype==noaurundesc
        outputallatomizcoords(sys,dt,1:2,path)
    end

    if !simplerun
        # output animation of simulation in $path
        outputanimation(sys,path)

        outputsysE(sys,dt,systype,path)
        outputallsysforces(sys,dt,path)
        outputallsysvelocities(sys,dt,path)
    end
end

############################################################################################################

"""
helper function to make folder for results
"""
function makeresultsfolder(desc::String,steps)
   date=Dates.format(now(), "yyyy-mm-ddTHHMMSS")
   mkpath("results/$desc-$steps--$date")
end

function makeresultsfolder(desc::String)
    date=Dates.format(now(), "yyyy-mm-ddTHHMMSS")
    mkpath("$desc--$date")
 end

############################################################################################################

"""
run MD on a system and output run info to path
"""
function runMDprintresults(sys::System,desc::String,simulator,steps::Int64,path::String=makeresultsfolder(desc,steps))
   # time step in simulation
   dt=simulator.dt

   # run MD+give run time
   runtime=@elapsed simulate!(sys, simulator, steps)
   runtime*=u"s"

   # output all system data: animation, coords, last velocities/forces
   outputsysinfo(sys,dt,desc,path,multirun)

   # output summary of run
   outputsummary(sys,dt,desc,steps,runtime,path)
end

############################################################################################################

function checkscattering(s::System)
    cutoff=10u"Å"+maximum(au.z)
    finalzN=s.loggers.coords.history[end][1][3]

    if finalzN>=cutoff
        true
    else
        false
    end
end

############################################################################################################

function checkEconserved(s::System)
    initialE=s.loggers.et.history[1]
    finalE=s.loggers.et.history[end]
    percentdif=abs(initialE-finalE)/initialE

    if percentdif<0.01
        true
    else
        false
    end
end

############################################################################################################

"""
get final NO coords/vel + energies. output as vec in MD units
"""
function finalE_NO(s::System)
    # NO masses
    mN=s.atoms[1].mass
    mO=s.atoms[2].mass
    mt=mN+mO
    μNO=mN*mO/mt

    # NO last coords
    finalrN=s.loggers.coords.history[end][1] .|> u"d_MD"
    finalrO=s.loggers.coords.history[end][2] .|> u"d_MD"
    rvec_NO=vector(finalrN,finalrO,s.boundary)
    ur_NO=normalize(rvec_NO)

    # NO last velocities
    finalvN=s.loggers.velocities.history[end][1] .|> u"v_MD"
    finalvO=s.loggers.velocities.history[end][2] .|> u"v_MD"
    vvec_NO=finalvO-finalvN

    # NO final translational energy
    vcom_sq=sum([getCOM(mN,mO,finalvN[i],finalvO[i])^2 for i in eachindex(finalvN)])
    finalEtrans=0.5*mt*vcom_sq*N_A |> u"e_MD"

    # NO final vibrational KE
    vvib_sq=dot(ur_NO,vvec_NO)^2
    finalKEvib=0.5*μNO*vvib_sq*N_A |> u"e_MD"

    # NO final total KE/rotational energy
    finalKEtot=0.5*sum(mN.*finalvN.^2+mO.*finalvO.^2)*N_A |> u"e_MD"
    finalErot=finalKEtot-finalEtrans-finalKEvib

    # output final NO coords, velocities, energies
    Evec=[finalKEtot,finalEtrans,finalErot,finalKEvib]
    ustrip_vec(vcat(finalrN,finalrO,finalvN,finalvO,Evec))
end

############################################################################################################

"""helper function: zip folders in dir"""
function zipfolders(dirpath::String)
    curdir=pwd()
    cd(dirpath)
    dircontents=readdir()
    archivename="traj_details.zip"
    run(`zip -r $archivename $dircontents`)

    for i in eachindex(dircontents)
        rm(dircontents[i],recursive=true)
    end
    cd(curdir)
end

############################################################################################################

"""
print txt file with summary of the multiple runs.
"""
function outputmultirunsummary(vary,counttraj,runpath=".")
    file="$runpath/summary.txt"
    open(file,"w") do io
        println(io,"Summary of multiple runs")
        println(io)
        println(io,"Variables varied: $vary")
        println(io,"Ran on ISAAC: $isaac")
        println(io)
        println(io,counttraj)
    end
end

############################################################################################################

"""
print txt file with summary of the multiple runs.
"""
function outputtrajinfo(trajscatter::DataFrame,trajtrap::DataFrame,runpath=".")
    name_sc="$runpath/traj_scattered.xlsx"
    name_tr="$runpath/traj_trapped.xlsx"
    write_xlsx(name_sc,trajscatter)
    write_xlsx(name_tr,trajtrap)
end
