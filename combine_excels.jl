using XLSX
using DataFrames

# Define directory path
dir_path = "path/to/directory"

# Define unique folder names
folder_names = ["NO-Au_sc-ISAAC-fixorient", "NO-Au_sc-ISAAC-fixorient-run", "NO-Au_sc-ISAAC-normalrun"]

# Define unique file names
file_names = ["traj_summary.xlsx", "traj_trapped.xlsx", "traj_scattered.xlsx"]

# Loop through unique folder names and file names
for folder_name in folder_names
    for file_name in file_names
        # Define DataFrame to store all data from Excel files
        all_data = DataFrame()

        # Loop through all directories and files matching current folder and file names
        for path in readdir(joinpath(dir_path, folder_name))
            if isdir(joinpath(dir_path, folder_name, path))
                # Get path to current Excel file
                excel_file = joinpath(dir_path, folder_name, path, file_name)

                # Check if current Excel file exists
                if isfile(excel_file)
                    # Read current Excel file into DataFrame
                    data = DataFrame(XLSX.readdata(excel_file, "Sheet1")...)
                    
                    # Add directory name to DataFrame
                    data.directory = path

                    # Append current DataFrame to all_data DataFrame
                    append!(all_data, data)
                end
            end
        end

        # Check if any data was read in
        if size(all_data, 1) > 0
            # Write combined DataFrame to new Excel file
            XLSX.writetable(joinpath(dir_path, "$(folder_name)_$(file_name)"), "Sheet1", all_data)
        end
    end
end
