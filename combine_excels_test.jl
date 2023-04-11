using XLSX
using DataFrames

# A function to calculate weighted averages given two columns 
function weighted_avg(col, weight_col)
    dot(col, weight_col) / sum(weight_col)
end

# The main function that combines the xlsx files for each type
function combine_summaries(dir_path::AbstractString)
    # Define unique folder names
    folder_names = ["NO-Au_sc-ISAAC-normalrun", "NO-Au_sc-ISAAC-fixorient", "O-Au_sc-ISAAC-10000"]

    cd(dir_path)

    # Initialize a dataframe vector for each folder name
    dfvec = [DataFrame() for _ in folder_names]

    for (i, folder_name) in enumerate(folder_names)
        # Find all subdirectories matching the current folder name
        type_folders = filter(x -> occursin(folder_name, x), readdir())

        # Iterate over each subdirectory and read the excel file if it exists
        for type_folder in type_folders
            summary_path = joinpath(type_folder, "traj_summary.xlsx")
            if !isfile(summary_path)
                continue
            end

            summary_data = DataFrame(XLSX.readtable(summary_path, 1))

            # Keep only columns before the one containing "avg"
            avg_col_idx = findfirst(x -> occursin(r"avg", x), names(summary_data))
            val_headers = names(summary_data)[1:avg_col_idx-1]

            # Compute n_scatter and n_trap sums
            n_scatter, n_trap = sum(summary_data.n_scatter), sum(summary_data.n_trap)

            if n_scatter == 0
                continue # Skip empty dataframes
            end

            # Group the data by the remaining columns and compute the weighted average of the "avg" columns
            groupby_cols = [Symbol(val_header) for val_header in val_headers]
            summary_agg = combine(groupby(summary_data, groupby_cols), :n_scatter => sum, :n_trap => sum)

            for col in names(summary_data)
                if occursin(r"avg", col)
                    avg_val = weighted_avg(summary_data[!,col], summary_data.n_scatter)
                    push!(summary_agg, avg_val)
                end
            end

            # Compute additional columns and append to the respective dataframe vector
            summary_agg[:n_total] = n_scatter+n_trap
            summary_agg[:frac_scatter] = n_scatter / summary_agg.n_total
            summary_agg[:frac_trap] = n_trap / summary_agg.n_total
            summary_agg[:type] = folder_name # Add the folder name as a column
            append!(dfvec[i], summary_agg)
        end
    end

    cd("..") # Return to the original directory

    outputtrajinfo(dfvec) # Output the combined dataframes to separate sheets
end

# A function that outputs each dataframe to an xlsx file in the current directory
function outputtrajinfo(dfvec::Vector{DataFrame})
    runpath = "."
    for (i, df) in enumerate(dfvec)
        name = ["summary_", lowercase(replace(folder_names[i], r"[\W\d_]" => ""))]
        write_xlsx(joinpath(runpath, join(name)) * ".xlsx", df)
    end
end

# A function to write a dataframe to an xlsx file
function write_xlsx(file_path::AbstractString, df::DataFrame)
    XLSX.writetable(file_path, collect(eachcol(df)), names(df))
end