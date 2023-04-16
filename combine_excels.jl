# Define directory path
dir_path = expanduser("~/scratch/NO-Au-results/run-4-13")

############################################################################################################

# i have a folder (test_combine_excel) with folders, each of a particular naming pattern (type).
# they might start with one of the following: NO-Au_sc-ISAAC-fixorient, NO-Au_sc-ISAAC-normalrun, O-Au_sc-ISAAC-10000 
# for each folder in test_combine_excel, there might be a file in it named traj_summary.xlsx
# each xlsx file has the same headers, but different data

# write code in julia that will output a combined xlsx file for each type
# use XLSX and DataFrames pkgs.

using XLSX
using DataFrames
using LinearAlgebra

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

############################################################################################################

t=@elapsed begin
    curdir=pwd()
    cd(dir_path)
    dir_folders=readdir()

    # Define unique folder names
    folder_names = ["NO-Au_sc-ISAAC-normalrun", "NO-Au_sc-ISAAC-fixorient", "O-Au_sc-ISAAC-10000"]
    filtered_folders = [filter(x -> contains(x, name), dir_folders) for name in folder_names]

    # Define unique file names
    file_name = "traj_summary.xlsx"

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
