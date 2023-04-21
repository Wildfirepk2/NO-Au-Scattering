# Define directory path
# dir_path = expanduser("~/scratch/NO-Au-results/run-4-13")
dir_path = expanduser("~/Downloads/test_combine_excel")
curdir=pwd()

############################################################################################################

using XLSX
using DataFrames
using CSV
using CairoMakie
using Dates
using StatsBase

############################################################################################################

"""
helper function: get path of au folder in directory. if not found, return nothing
"""
function getAuDirPath(T)
    aurundesc="Au_slab-T_$T"
    resultsfolder=readdir("$curdir/results";join=true)
    i_au=findfirst(contains.(resultsfolder,aurundesc))
    i_au===nothing ? error("Directory for T = $T not found") : resultsfolder[i_au]
end

"""
get equilibrated au coords from run at temp T
"""
function getEquilAuCoords(T)
    audir=getAuDirPath(T)
    if "syscoords.xlsx" in readdir(audir)
       coordsfile="$audir/syscoords.xlsx"
       xfcoord=XLSX.readxlsx(coordsfile)
       sheets=XLSX.sheetnames(xfcoord)
       sheetlastcoord=sheets[end]
       DataFrame(XLSX.readtable(coordsfile,sheetlastcoord))
    else
       readcoorddir=readdir("$audir/syscoords") # incorrectly sorts by number. eg 1,10,11,12,...,2,etc
       sortedcoorddir = sort(readcoorddir, by = x -> parse(Int, match(r"\d+", x).match)) # sorted by name. eg syscoords 1,2,3,4,etc
       lastcoords=sortedcoorddir[end]
       CSV.read("$audir/syscoords/$lastcoords",DataFrame)
    end
end

function getfilteredAu(T)
    df = getEquilAuCoords(T)
    # filter df for z col values greater than 7
    df_filtered = filter(row -> row.z > 6.5, df) # jump in z from ~5 to ~7
    # remove z col 
    df_filtered = select(df_filtered, Not(:z))
    convert.(Float64, df_filtered) # should have 132 rows (528 atoms/4 layers)
end

function filtertrajxy(df::DataFrame,temp,orient)
    e_filter=minimum(df.Ei)
    if "θorient" in names(df)
        # filter by orient and e_filter
        filter!(row -> row.θorient == orient && row.Ei == e_filter, df)
        # var="θ$orient"
    else
        # filter by temp and e_filter
        filter!(row -> row.T == temp && row.Ei == e_filter, df)
        # var="T$temp"
    end
    # delete cols not matching "xcom" or "ycom"
    select!(df, "xcom", "ycom")

    # rename cols to x and y
    rename!(df, "xcom"=>"x", "ycom"=>"y")

    # change col types
    df_new=convert.(Float64, df)

    # return df_new, var
    return df_new
end

function makeresultsfolder(desc::String)
    date=Dates.format(now(), "yyyy-mm-ddTHHMMSS-sss")
    path=mkpath("$desc--$date")
    println("Made folder: $path")
    return path
end

function outputgraph(dfvec,systype,temp,orient)
    # get au top layer xy coords
    df_au=getfilteredAu(temp)

    # plot au coords
    f = Figure(resolution = (800, 600))
    ax=Axis(f[1,1],xlabel="x (Å)",ylabel="y (Å)")
    xylim=[(-1.5, 34),(-0.5, 30.5)]#\fix. unknown error giving ax limits as (0,10). happens only if fig gen in loop
    limits!(ax,xylim...)
    scatter!(df_au.x,df_au.y,
            color=:gold,
            markersize=50,
            )

    # overlay heatmap
    nbin=30
    xyedge=Tuple(range(xy...,nbin) for xy in xylim)
    colors=[:Reds, :Blues]
    colorbarlabels=["Scattered", "Trapped"]
    transparancy_fac=count(!isempty, dfvec)==1 ? 0.7 : 0.45 # may \fix
    # figobjs=[]
    for (df,color,label) in zip(dfvec,colors,colorbarlabels)
        if !isempty(df)
            h = fit(Histogram, (df.x,df.y),xyedge)
            hm=heatmap!(h,
                    colormap=(color,transparancy_fac),
                    interpolate=false,
                    # label=label,
                    # transparency = true,
                    )
            Colorbar(f[:, end+1], hm, label=label,)
        end
    end

    # save fig
    var=contains(systype,"fix") ? "θ$orient" : "T$temp"
    save("$path/heatmap-$systype-$var.png",f)
end
 
xlims(ax::Axis=current_axis()) =
    (ax.finallimits[].origin[1], ax.finallimits[].widths[1])
ylims(ax::Axis=current_axis()) =
    (ax.finallimits[].origin[2], ax.finallimits[].widths[2])

############################################################################################################

t=@elapsed begin
    cd(dir_path)
    dir_folders=readdir()

    # Define unique folder names
    folder_names = ["NO-Au_sc-ISAAC-normalrun", "NO-Au_sc-ISAAC-fixorient", "O-Au_sc-ISAAC-10000"]
    systypes=["noau_norm", "noau_fix", "oau"]
    filtered_folders = [filter(x -> contains(x, name), dir_folders) for name in folder_names]

    # Define unique file names
    file_names = ["traj_scattered.xlsx", "traj_trapped.xlsx"]
    filetypes=["sc", "tr"]

    ts=[300, 400, 500]
    orients=[0, 180, 90]

    # dftraj=[[DataFrame() for _ in file_names] for _ in folder_names] # 1st vec is noau_norm/2nd vec is noau_fix/3rd vec is oau

    path=makeresultsfolder("heatmaps")
    # for (folderv,dfv) in zip(filtered_folders,dftraj), folder in folderv, (df,file) in zip(dfv,file_names)#\old
    # iterate through each sys type, then each T OR orient
    for (systype,folderv) in zip(systypes,filtered_folders), (t,orient) in zip(ts,orients)
        # container for systype-T/orient combo
        dfvec=[DataFrame() for _ in file_names]
        # iterate through each folder corresponding to systype, then each type of file: tr and sc
        for folder in folderv, (df,file) in zip(dfvec,file_names)
            # add to corresponding df in dfvec
            target_file = folder*"/"*file
            if isfile(target_file)
                data = DataFrame(XLSX.readtable(target_file, 1))
                filtered_data=filtertrajxy(data,t,orient)
                append!(df,filtered_data)
            end
        end
        try
            outputgraph(dfvec,systype,t,orient)
            # break
        catch e
            cd(curdir)
            error("Error occurred: $e")
        end
        # # save dfvec to file
        # for i in eachindex(dfvec)
        #     CSV.write("$path/$systype-$(filetypes[i])-T$t-θ$orient.csv",dfvec[i])
        # end
    end
    cd(curdir)
end

println("Elapsed time: $t seconds")