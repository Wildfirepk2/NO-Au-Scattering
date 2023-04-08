# functions to analyze Molly simulations after running

############################################################################################################

# extend length method for XLSXFile to give number of sheets. may revert later.
Base.length(xf::XLSX.XLSXFile) = XLSX.sheetcount(xf)

"""
write data to excel file. intermediate function for output functions
"""
function write_xlsx(file::String, dt, i, stepslog, df::DataFrame)
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
function write_xlsx(file::String, df::DataFrame)
    df=ustrip.(df)
    data = collect(eachcol(df))
    cols = names(df)
    XLSX.writetable(file, data, cols)
end

############################################################################################################

"""
output atom coords at ith logging point to excel file. default path is current directory. used only in outputallsyscoords but can be called independently
"""
function outputsyscoords(sys::System,dt,i,path=".",asexcel=true)
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
function outputallsyscoords(sys::System,dt,path=".",asexcel=true,firstlastonly=true)
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
output atom z coords (atoms #'s specified by vec) w time to excel file + make graph
"""
function outputallatomizcoords(sys::System,dt,atomrange,path=".")
    # logging interval. ie log coords after every $stepslog steps
    stepslog=sys.loggers.coords.n_steps

    # number of logs
    nsteps=length(sys.loggers.coords.history)

    # time between logs
    dtlog=stepslog*dt

    # tmp vars for time/coords + converted to MD units
    z_surf=maximum(au.z)
    time=[(i-1)*dtlog for i in 1:nsteps] .|> u"t_MD" # i-1 because first log is at step 0
    
    data=DataFrame(t=time)
    for atom_i in atomrange
        label_i=getatomilabel(sys,atom_i)
        zcoords=[sys.loggers.coords.history[i][atom_i][3]-z_surf for i in eachindex(sys.loggers.coords.history)] .|> u"d_MD"
        data=hcat(data,DataFrame(label_i=>zcoords))
    end
    
    # write to excel file at $path
    file="$path/z-z_surf v time.xlsx"
    write_xlsx(file,data)
    
    # make zcoord vs time graph
    graphdesc = getzgraphdesc(sys)
    outputgraph(graphdesc,data,path)
end

############################################################################################################

"""
output graph given DataFrame. called in outputsysE
"""
function outputgraph(desc,df::DataFrame,path=".";stackplots=false)
    # grabbing key data from df
    colnames=names(df)
    numcols=ncol(df)
    dffirstrow=first.(eachcol(df)) # \fix try eltype later?
    colunittypes=qtytolabel.(dffirstrow)
    colunits=getgraphunitlabel.(dffirstrow)
    xunittype=colunittypes[begin]
    yunittype=colunittypes[begin+1]
    xunit=colunits[begin]
    yunit=colunits[begin+1]

    # convert units then unit strip all df cols
    for i in eachindex(colunits)
        df[!,i]=ustrip.(colunits[i],df[!,i])
    end
    
    # create plot objs
    f = Figure(resolution = (900, 600)) # default 800x600
    Axis(f[1, 1], title = desc, xlabel = "$xunittype ($xunit)", ylabel = "$yunittype ($yunit)")
    plotobjs=[]

    # add to plot obj
    if stackplots
        # \fix
    else
        for i in 1:numcols-1
            obj=scatterlines!(df[!,begin],df[!,i+1],label=colnames[i+1],markersize=5,linewidth=0.5)
            push!(plotobjs,obj)
        end
    end
    if length(plotobjs)>1;Legend(f[1, 2],plotobjs,colnames[2:end]);end
    
    # save plot to $path
    save("$path/$yunittype v $xunittype.png", f)
    return f
end

############################################################################################################
    
"""helper function: get description for charge graph"""
function getchargegraphdesc(s::(System{D, false, T, CU, A, AD, PI} where {D, T, CU, A, AD, PI <: Tuple{NOAuInteraction}}))
    "Charge on NO vs Time"
end

function getchargegraphdesc(s::(System{D, false, T, CU, A, AD, PI} where {D, T, CU, A, AD, PI <: Tuple{OAuInteraction}}))
    "Charge on O vs Time"
end

############################################################################################################
#\fix
"""
graft charge to existing figure
"""
function outputcharge(f::Figure,sys::System,dt,path=".")
    fc=outputcharge(sys,dt,path)
    # mod fc axis to be on right

    #stitch f./fc graphs 

    # example:
    f = Figure()

    ax1 = Axis(f[1, 1], yticklabelcolor = :blue)
    ax2 = Axis(f[1, 1], yticklabelcolor = :red, yaxisposition = :right)

    lines!(ax1, 0..10, sin, color = :blue)
    lines!(ax2, 0..10, x -> 100 * cos(x), color = :red)

    f
end

"""
make standalone charge graph
"""
function outputcharge(sys::System,dt,path=".")
    charges=sys.loggers.charge.history .|> u"e⁻"

    # logging interval. ie log coords after every $stepslog steps
    stepslog=sys.loggers.charge.n_steps

    # number of logs
    nsteps=length(charges)

    # time between logs
    dtlog=stepslog*dt

    # tmp var for time + converted to MD units
    time=[(i-1)*dtlog for i in 1:nsteps] .|> u"t_MD" # i-1 because first log is at step 0
    
    data=DataFrame(t=time,Charge=charges)
    
    # write to excel file at $path
    file="$path/charge v time.xlsx"
    write_xlsx(file,data)
    
    # make zcoord vs time graph
    graphdesc = getchargegraphdesc(sys)
    outputgraph(graphdesc,data,path)
end

############################################################################################################
#\fix
"""
output z coords vs charge on 1 graph
"""
function outputzcoordvscharge(sys::System,dt,atomrange,path=".")
    f=outputallatomizcoords(sys,dt,atomrange,path)
    outputcharge(f,sys,dt,path)
end

############################################################################################################

"""
output system energies w time to excel file + make graph
"""
function outputsysE(sys,dt,path=".")
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
    graphdesc = getEtitle(sys)

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
helper function to make folder for results
"""
function makeresultsfolder(desc::String,steps)
    date=Dates.format(now(), "yyyy-mm-ddTHHMMSS")
    path=mkpath("results/$desc-$steps--$date")
    println("Made folder: $path")
    return path
end

function makeresultsfolder(desc::String)
    date=Dates.format(now(), "yyyy-mm-ddTHHMMSS")
    path=mkpath("$desc--$date")
    println("Made folder: $path")
    return path
 end

############################################################################################################

"""
run MD on a system and output run info to path
"""
function runMDprintresults(sys::System,desc::String,simulator,steps::Int64,path::String=makeresultsfolder(desc,steps),TEi...)
    # time step in simulation
    dt=simulator.dt

    # run MD+give run time
    runtime=@elapsed simulate!(sys, simulator, steps)
    runtime*=u"s"

    #\fix. if no au, error: no velocities
    # # output all system data: animation, coords, last velocities/forces
    # outputsysinfo(sys,dt,path)

    # output summary of run
    outputsummary(sys,dt,TEi...,steps,runtime,path)
end

# #\fix later
# function runMDprintresults(sys::System{D, false, T, CU, A, AD, PI} where {D,T,CU,A,AD,PI<:Tuple{AuSlabInteraction}},desc::String,simulator,steps::Int64,T,path::String=makeresultsfolder(desc,steps))
#     Ei=NaN
#     runMDprintresults(sys,desc,simulator,steps,path,T,Ei)
# end

############################################################################################################

function checkEconserved(s::System)
    initialE=s.loggers.et.history[1]
    finalE=s.loggers.et.history[end]
    percentdif=abs(initialE-finalE)/initialE

    percentdif<0.01 ? println("Energy conserved") : error("Energy not conserved. Initial E: $initialE, Final E: $finalE, Percent difference: $percentdif")
end
