# Define directory path
dir_path = expanduser("~/scratch/NO-Au-results/runs--4-5")
# dir_path = expanduser("~/Downloads/test_combine_excel")

############################################################################################################

using XLSX
using DataFrames
using LinearAlgebra
using StatsBase

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
    summary_file_name = "traj_summary.xlsx"
    sc_file_name = "traj_scattered.xlsx"
    out_name = "traj_summary_with_SD"

    for folder in dir_folders
        summary_file = folder*"/"*summary_file_name
        try
            if isfile(summary_file)
                summarydata = DataFrame(XLSX.readtable(summary_file, 1))
    
                # find index of column name containing "avg"
                it = findfirst(x -> occursin(r"avg", x), names(summarydata))
    
                # get headers for columns before "avg"
                valheaders = names(summarydata)[1:it-1]
    
                # create iterator for unique values of each column
                iter = [unique(summarydata[!,var]) for var in valheaders]
    
                # erot std dev
                erot_SD_df = DataFrame()
                sc_file = folder*"/"*sc_file_name
                ecolname = contains(folder, "NO-Au") ? "Erot" : "Etransfer"
                if isfile(sc_file)
                    # read scattered file
                    scdata = DataFrame(XLSX.readtable(sc_file, 1))
    
                    # filter scdata for each combination of iter
                    for values in Iterators.product(iter...)
                        # create dictionary of column name-value pairs for current iteration
                        dict = Dict(zip(valheaders, values))
    
                        # filter df using dictionary values
                        filtered_df = filter(row -> all(row[key] == value for (key, value) in dict), scdata)
    
                        # calc std dev of $ecolname
                        erot_SD = std(Float64.(filtered_df[!,ecolname]))
                        ecolname_SD = ecolname*"_SD"
    
                        # maindf=DataFrame(dict)
                        # sddf=DataFrame(ecolname_SD=>erot_SD)
                        # finaldf=hcat(maindf,sddf)
                        # append!(erot_SD_df, finaldf)
    
                        sddf=DataFrame(ecolname_SD=>erot_SD)
                        append!(erot_SD_df, sddf)
                    end
                else
                    for values in Iterators.product(iter...)
                        # maindf=DataFrame(dict)
                        # sddf=DataFrame(ecolname_SD=>0)
                        # finaldf=hcat(maindf,sddf)
                        # append!(erot_SD_df, finaldf)
    
                        sddf=DataFrame(ecolname_SD=>0)
                        append!(erot_SD_df, sddf)
                    end
                end
    
                # trap SD
                trapt=[vcat(zeros(n_sc),ones(n_tr)) for (n_sc,n_tr) in zip(summarydata.n_scatter,summarydata.n_trap)]
                trap_SD=[std(vec) for vec in trapt]
                trap_SD_df=DataFrame(trap_SD=trap_SD)
    
                finaldf=hcat(summarydata,erot_SD_df,trap_SD_df)
                out_file = folder*"/"*out_name*".xlsx"
                if !isfile(out_file)
                    outputtrajinfo(out_name,finaldf,folder)
                end
                # break 
            end
        catch e
            cd(curdir)
            error("Error occurred: $e")
        end
    end
    cd(curdir)
end

println("Elapsed time: $t seconds")
