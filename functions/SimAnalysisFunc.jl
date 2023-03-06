# functions to analyze Molly simulations after running

############################################################################################################

"""
copy of molly's visualize fn (in src/makie.jl) for use w CairoMakie. used if isaac=true
"""
function myvisualize(coord_logger,
                    boundary,
                    out_filepath::AbstractString;
                    connections=Tuple{Int, Int}[],
                    connection_frames=[trues(length(connections)) for i in values(coord_logger)],
                    trails::Integer=0,
                    framerate::Integer=30,
                    color=:purple,
                    connection_color=:orange,
                    # markersize=0.001,
                    markersize=0.0005,
                    linewidth=2.0,
                    transparency=true,
                    show_boundary::Bool=true,
                    boundary_linewidth=2.0,
                    boundary_color=:black,
                    kwargs...)
    coords_start = first(values(coord_logger))
    dist_unit = unit(first(first(coords_start)))
    dims = n_dimensions(boundary)
    fig = Figure()

    if dims == 3
        PointType = Point3f
        # ax = Axis3(fig[1, 1], aspect=:data)
        ax = Axis3(fig[1, 1], aspect=(2.7,1.5,1))
        max_connection_dist = cbrt(box_volume(boundary)) / 2
    elseif dims == 2
        PointType = Point2f
        ax = Axis(fig[1, 1])
        ax.aspect = DataAspect()
        max_connection_dist = sqrt(box_volume(boundary)) / 2
    else
        throw(ArgumentError("Found $dims dimensions but can only visualize 2 or 3 dimensions"))
    end

    positions = Observable(PointType.(ustrip_vec.(coords_start)))
    scatter!(ax, positions; color=color, markersize=markersize, transparency=transparency,
                markerspace=:data, kwargs...)

    if show_boundary
        lines!(
            ax,
            bounding_box_lines(boundary, dist_unit)...;
            color=boundary_color,
            linewidth=boundary_linewidth,
        )
    end

    connection_nodes = []
    for (ci, (i, j)) in enumerate(connections)
        # Don't display connected atoms that are likely connected over the box edge
        if first(connection_frames)[ci] && norm(coords_start[i] - coords_start[j]) < max_connection_dist
            if dims == 3
                push!(connection_nodes, Observable(PointType.(
                        ustrip.([coords_start[i][1], coords_start[j][1]]),
                        ustrip.([coords_start[i][2], coords_start[j][2]]),
                        ustrip.([coords_start[i][3], coords_start[j][3]]))))
            elseif dims == 2
                push!(connection_nodes, Observable(PointType.(
                        ustrip.([coords_start[i][1], coords_start[j][1]]),
                        ustrip.([coords_start[i][2], coords_start[j][2]]))))
            end
        else
            if dims == 3
                push!(connection_nodes, Observable(PointType.([0.0, 0.0], [0.0, 0.0],
                                                        [0.0, 0.0])))
            elseif dims == 2
                push!(connection_nodes, Observable(PointType.([0.0, 0.0], [0.0, 0.0])))
            end
        end
    end
    for (ci, cn) in enumerate(connection_nodes)
        lines!(ax, cn;
                color=isa(connection_color, AbstractArray) ? connection_color[ci] : connection_color,
                linewidth=isa(linewidth, AbstractArray) ? linewidth[ci] : linewidth,
                transparency=transparency)
    end

    trail_positions = []
    for trail_i in 1:trails
        push!(trail_positions, Observable(PointType.(ustrip_vec.(coords_start))))
        col = parse.(Colorant, color)
        alpha = 1 - (trail_i / (trails + 1))
        alpha_col = RGBA.(red.(col), green.(col), blue.(col), alpha)
        scatter!(ax, trail_positions[end]; color=alpha_col,  markersize=markersize,
                    transparency=transparency, markerspace=:data, kwargs...)
    end

    boundary_conv = ustrip.(dist_unit, Molly.cubic_bounding_box(boundary))
    xlims!(ax, Molly.axis_limits(boundary_conv, coord_logger, 1))
    ylims!(ax, Molly.axis_limits(boundary_conv, coord_logger, 2))
    dims == 3 && zlims!(ax, Molly.axis_limits(boundary_conv, coord_logger, 3))

    record(fig, out_filepath, eachindex(values(coord_logger)); framerate=framerate) do frame_i
        coords = values(coord_logger)[frame_i]

        for (ci, (i, j)) in enumerate(connections)
            if connection_frames[frame_i][ci] && norm(coords[i] - coords[j]) < max_connection_dist
                if dims == 3
                    connection_nodes[ci][] = PointType.(
                                ustrip.([coords[i][1], coords[j][1]]),
                                ustrip.([coords[i][2], coords[j][2]]),
                                ustrip.([coords[i][3], coords[j][3]]))
                elseif dims == 2
                    connection_nodes[ci][] = PointType.(
                                ustrip.([coords[i][1], coords[j][1]]),
                                ustrip.([coords[i][2], coords[j][2]]))
                end
            else
                if dims == 3
                    connection_nodes[ci][] = PointType.([0.0, 0.0], [0.0, 0.0],
                                                        [0.0, 0.0])
                elseif dims == 2
                    connection_nodes[ci][] = PointType.([0.0, 0.0], [0.0, 0.0])
                end
            end
        end

        positions[] = PointType.(ustrip_vec.(coords))
        for (trail_i, trail_position) in enumerate(trail_positions)
            trail_position[] = PointType.(ustrip_vec.(values(coord_logger)[max(frame_i - trail_i, 1)]))
        end
    end
end

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
function outputsyscoords(sys,dt,i,path=".")
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
    write_xlsx(file,dt,i,stepslog,data)
end

############################################################################################################

"""
output atom coords w time to excel file at path. default path is current directory. used only in outputsysinfo but can be called independently

firstlastonly for speed. XLSX slow.
"""
function outputallsyscoords(sys,dt,path=".",firstlastonly=true)
    # coords from molly's loggers
    datasrc=sys.loggers.coords.history

    if firstlastonly
        # write excel file for first/last time step only.
        outputsyscoords(sys,dt,1,path)
        outputsyscoords(sys,dt,length(datasrc),path)
    else
        # write excel file for ea time step
        for i in eachindex(datasrc)
            outputsyscoords(sys,dt,i,path)
        end
    end
end

############################################################################################################

"""
output atom i z coords w time to excel file + make graph
"""
function outputallatomizcoords(sys,dt,atom_i,path=".")
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
        # visualize(sys.loggers.coords, sys.boundary, "$path/animation.mp4";show_boundary=false)
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
        # visualize(sys.loggers.coords, sys.boundary, "$path/animation.mp4";show_boundary=false,connections=connections,connection_color=connection_color,color=syscolors,)
    end
end

############################################################################################################

"""
output all system info. animation. logger quantities

run description needs to contain either "Au slab" or "NO/Au"
"""
function outputsysinfo(sys,dt,systype,path=".")
    # output animation of simulation in $path
    outputanimation(sys,path)

    # output quantities to excel file in separate folder
    outputallsyscoords(sys,dt,path)
    if systype==noaurundesc
        outputallatomizcoords(sys,dt,1,path)
    end
    outputsysE(sys,dt,systype,path)
    outputallsysforces(sys,dt,path)
    outputallsysvelocities(sys,dt,path)
end

############################################################################################################

"""
run MD on a system and output run info to results folder
"""
function runMDprintresults(sys,desc,simulator,steps)
    # time of run
    date=Dates.format(now(), "yyyy-mm-ddTHHMMSS")
    
    # make folder to store results of current run
    path=mkpath("results/$desc-$steps--$date")

    # time step in simulation
    dt=simulator.dt

    # run MD+give run time
    runtime=@elapsed simulate!(sys, simulator, steps)
    runtime*=u"s"

    # output all system data: animation, coords, last velocities/forces
    outputsysinfo(sys,dt,desc,path)

    # output summary of run
    outputsummary(sys,dt,desc,steps,runtime,path)
end
