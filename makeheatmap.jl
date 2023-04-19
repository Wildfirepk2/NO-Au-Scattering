# Define directory path
# dir_path = expanduser("~/scratch/NO-Au-results/run-4-13")
dir_path = expanduser("~/Downloads/test_combine_excel")

############################################################################################################

using XLSX
using DataFrames
using CSV
using CairoMakie
using Unitful
using Dates
include("functions/MDUnits.jl")

############################################################################################################

"""
print excel files with detailed trajectory info of the multiple runs.
"""
function outputtrajinfo(namevec::Vector{String},dfvec::Vector{DataFrame},runpath::String=".")
    for (name,df) in zip(namevec,dfvec)
        outputtrajinfo(name,df,runpath)
    end
end

function outputtrajinfo(name::String,df::DataFrame,runpath::String=".")
    if !isempty(df)
        name="$runpath/$name.xlsx"
        write_xlsx(name,df)
    end
end

function write_xlsx(file::String, df::DataFrame)
    data = collect(eachcol(df))
    cols = names(df)
    XLSX.writetable(file, data, cols)
end

"""
helper function: get path of au folder in directory. if not found, return nothing
"""
function getAuDirPath(T)
    aurundesc="Au_slab-T_$T"
    resultsfolder=readdir("results";join=true)
    i_au=findfirst(contains.(resultsfolder,aurundesc))
    i_au isa Nothing ? i_au : resultsfolder[i_au]
end

"""
get equilibrated au coords from previous run
"""
function getEquilAuCoords()
    audir=getAuDirPath("results")
    if "syscoords.xlsx" in readdir(audir)
       coordsfile="$audir/syscoords.xlsx"
       xfcoord=XLSX.readxlsx(coordsfile)
       sheets=XLSX.sheetnames(xfcoord)
       sheetlastcoord=sheets[end]
       dfcoord=DataFrame(XLSX.readtable(coordsfile,sheetlastcoord))
       [SA[dfcoord[i,1],dfcoord[i,2],dfcoord[i,3]]u"d_MD" for i in 1:nrow(dfcoord)]
    else
       readcoorddir=readdir("$audir/syscoords") # incorrectly sorts by number. eg 1,10,11,12,...,2,etc
       sortedcoorddir = sort(readcoorddir, by = x -> parse(Int, match(r"\d+", x).match)) # sorted by name. eg syscoords 1,2,3,4,etc
       lastcoords=sortedcoorddir[end]
       dfcoord=CSV.read("$audir/syscoords/$lastcoords",DataFrame)
       [SA[dfcoord[i,1],dfcoord[i,2],dfcoord[i,3]]u"d_MD" for i in 1:nrow(dfcoord)]
    end
end

############################################################################################################

t=@elapsed begin
    curdir=pwd()
    cd(dir_path)
    dir_folders=readdir()

    # get au coords

    # filter au df for top layer

    # Define unique folder names
    folder_names = ["NO-Au_sc-ISAAC-normalrun", "NO-Au_sc-ISAAC-fixorient", "O-Au_sc-ISAAC-10000"]
    filtered_folders = [filter(x -> contains(x, name), dir_folders) for name in folder_names]

    # Define unique file names
    file_name = ["traj_scattered.xlsx", "traj_trapped.xlsx"]

    dfvec=[DataFrame() for _ in folder_names]

    for (df,type) in zip(dfvec,filtered_folders)
        df_type=DataFrame()
        for folder in type
            target_file = folder*"/"*file_name
            if isfile(target_file)
                data = DataFrame(XLSX.readtable(target_file, 1))
                append!(df_type, data)
            end
        end

        # if df_type still empty (no folders w 1 of the folder names), skip to next iteration
        if isempty(df_type); continue; end

        # find index of column name containing "avg"
        it = findfirst(x -> occursin(r"avg", x), names(df_type))

        # get headers for columns before "avg"
        valheaders = names(df_type)[1:it-1]

        # create iterator for unique values of each column
        iter = [unique(df_type[!,var]) for var in valheaders]

        # filter df for each combination of iter
        for values in Iterators.product(iter...)
            # create dictionary of column name-value pairs for current iteration
            dict = Dict(zip(valheaders, values))
            
            # filter df using dictionary values
            filtered_df = filter(row -> all(row[key] == value for (key, value) in dict), df_type)
            
            # do something with filtered_df
            maininfo=DataFrame(dict)
        
            n_scatter=sum(filtered_df.n_scatter)
            n_trap=sum(filtered_df.n_trap)
            n_total=n_scatter+n_trap
            frac_scatter=n_scatter/n_total
            frac_trap=n_trap/n_total
            otherinfo=DataFrame(n_scatter=n_scatter,
                                n_trap=n_trap,
                                n_total=n_total,
                                frac_scatter=frac_scatter,
                                frac_trap=frac_trap)
        
            # get average of columns that contain "avg" in header
            for col in names(filtered_df)
                if occursin(r"avg", col)
                    avg_col = filtered_df[!,col]
                    n_scatter_col = filtered_df[!,"n_scatter"]
                    avg_total = dot(avg_col,n_scatter_col) / n_scatter
                    dft=DataFrame(col=>avg_total)
                    maininfo=hcat(maininfo,dft)
                end
            end
            finaldf=hcat(maininfo,otherinfo)
            append!(df, finaldf)
        end
    end

    name_noau_n="summary_noau_normalrun"
    name_noau_fo="summary_noau_fixorient"
    name_oau="summary_oau"
    namevec=[name_noau_n,name_noau_fo,name_oau]
    outputtrajinfo(namevec,dfvec)
    cd(curdir)
end

println("Elapsed time: $t seconds")
