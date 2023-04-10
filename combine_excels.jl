# i have a directory (dir) with folders in it (shown below). 
# u can see that there are 3 unique naming patterns among the folders. 
# each folder can contain a max of 3 files: traj_summary.xlsx, traj_trapped.xlsx, or traj_scattered.xlsx. 
# each of those types of excel files have the same headers. 
# write a code in julia that will go through all the folders in dir and for each type of folder and for each type of excel file read them and output a combined version as excel into dir. 
# use XLSX and DataFrames pkgs.

# $ ls -l
# total 168
# drwxr-xr-x 2 btang7 tug2106 4096 Apr  6 13:18 NO-Au_sc-ISAAC-fixorient--2023-04-06T063035
# drwxr-xr-x 2 btang7 tug2106 4096 Apr  6 13:07 NO-Au_sc-ISAAC-fixorient--2023-04-06T063212
# drwxr-xr-x 2 btang7 tug2106 4096 Apr  6 13:15 NO-Au_sc-ISAAC-fixorient--2023-04-06T063857
# drwxr-xr-x 2 btang7 tug2106 4096 Apr  6 13:44 NO-Au_sc-ISAAC-fixorient--2023-04-06T064305
# drwxr-xr-x 2 btang7 tug2106 4096 Apr  6 13:59 NO-Au_sc-ISAAC-fixorient--2023-04-06T065624
# drwxr-xr-x 2 btang7 tug2106 4096 Apr  6 14:04 NO-Au_sc-ISAAC-fixorient--2023-04-06T065741
# drwxr-xr-x 2 btang7 tug2106 4096 Apr  6 13:42 NO-Au_sc-ISAAC-fixorient--2023-04-06T070301
# drwxr-xr-x 2 btang7 tug2106 4096 Apr  6 14:40 NO-Au_sc-ISAAC-fixorient--2023-04-06T071456
# drwxr-xr-x 2 btang7 tug2106 4096 Apr  6 14:51 NO-Au_sc-ISAAC-fixorient--2023-04-06T072145
# drwxr-xr-x 2 btang7 tug2106 4096 Apr  6 13:58 NO-Au_sc-ISAAC-fixorient--2023-04-06T072545
# drwxr-xr-x 2 btang7 tug2106 4096 Apr  6 14:58 NO-Au_sc-ISAAC-fixorient--2023-04-06T073117
# drwxr-xr-x 2 btang7 tug2106 4096 Apr  6 14:16 NO-Au_sc-ISAAC-fixorient--2023-04-06T073950
# drwxr-xr-x 2 btang7 tug2106 4096 Apr  6 14:23 NO-Au_sc-ISAAC-fixorient--2023-04-06T074439
# drwxr-xr-x 2 btang7 tug2106 4096 Apr  6 16:07 NO-Au_sc-ISAAC-fixorient--2023-04-06T080629
# drwxr-xr-x 2 btang7 tug2106 4096 Apr  6 06:30 NO-Au_sc-ISAAC-normalrun--2023-04-05T235135
# drwxr-xr-x 2 btang7 tug2106 4096 Apr  6 06:57 NO-Au_sc-ISAAC-normalrun--2023-04-05T235429
# drwxr-xr-x 2 btang7 tug2106 4096 Apr  6 06:43 NO-Au_sc-ISAAC-normalrun--2023-04-05T235643
# drwxr-xr-x 2 btang7 tug2106 4096 Apr  6 07:21 NO-Au_sc-ISAAC-normalrun--2023-04-06T000629
# drwxr-xr-x 2 btang7 tug2106 4096 Apr  6 07:14 NO-Au_sc-ISAAC-normalrun--2023-04-06T000651
# drwxr-xr-x 2 btang7 tug2106 4096 Apr  6 06:32 NO-Au_sc-ISAAC-normalrun--2023-04-06T000727
# drwxr-xr-x 2 btang7 tug2106 4096 Apr  6 06:38 NO-Au_sc-ISAAC-normalrun--2023-04-06T000731
# drwxr-xr-x 2 btang7 tug2106 4096 Apr  6 08:06 NO-Au_sc-ISAAC-normalrun--2023-04-06T000745
# drwxr-xr-x 2 btang7 tug2106 4096 Apr  6 06:56 NO-Au_sc-ISAAC-normalrun--2023-04-06T001025
# drwxr-xr-x 2 btang7 tug2106 4096 Apr  6 07:31 NO-Au_sc-ISAAC-normalrun--2023-04-06T001119
# drwxr-xr-x 2 btang7 tug2106 4096 Apr  6 07:03 NO-Au_sc-ISAAC-normalrun--2023-04-06T002952
# drwxr-xr-x 2 btang7 tug2106 4096 Apr  6 07:25 NO-Au_sc-ISAAC-normalrun--2023-04-06T005927
# drwxr-xr-x 2 btang7 tug2106 4096 Apr  6 07:39 NO-Au_sc-ISAAC-normalrun--2023-04-06T010726
# drwxr-xr-x 2 btang7 tug2106 4096 Apr  6 07:44 NO-Au_sc-ISAAC-normalrun--2023-04-06T010858
# drwxr-xr-x 2 btang7 tug2106 4096 Apr  6 17:37 O-Au_sc-ISAAC-10000--2023-04-06T130715
# drwxr-xr-x 2 btang7 tug2106 4096 Apr  6 17:47 O-Au_sc-ISAAC-10000--2023-04-06T131538
# drwxr-xr-x 2 btang7 tug2106 4096 Apr  6 18:07 O-Au_sc-ISAAC-10000--2023-04-06T131808
# drwxr-xr-x 2 btang7 tug2106 4096 Apr  6 18:15 O-Au_sc-ISAAC-10000--2023-04-06T134214
# drwxr-xr-x 2 btang7 tug2106 4096 Apr  6 18:37 O-Au_sc-ISAAC-10000--2023-04-06T134459
# drwxr-xr-x 2 btang7 tug2106 4096 Apr  6 18:29 O-Au_sc-ISAAC-10000--2023-04-06T135813
# drwxr-xr-x 2 btang7 tug2106 4096 Apr  6 19:11 O-Au_sc-ISAAC-10000--2023-04-06T135929
# drwxr-xr-x 2 btang7 tug2106 4096 Apr  6 19:03 O-Au_sc-ISAAC-10000--2023-04-06T140439
# drwxr-xr-x 2 btang7 tug2106 4096 Apr  6 18:50 O-Au_sc-ISAAC-10000--2023-04-06T141631
# drwxr-xr-x 2 btang7 tug2106 4096 Apr  6 19:07 O-Au_sc-ISAAC-10000--2023-04-06T142352
# drwxr-xr-x 2 btang7 tug2106 4096 Apr  6 19:51 O-Au_sc-ISAAC-10000--2023-04-06T144013
# drwxr-xr-x 2 btang7 tug2106 4096 Apr  6 20:08 O-Au_sc-ISAAC-10000--2023-04-06T145116
# drwxr-xr-x 2 btang7 tug2106 4096 Apr  6 20:06 O-Au_sc-ISAAC-10000--2023-04-06T145859
# drwxr-xr-x 2 btang7 tug2106 4096 Apr  6 21:51 O-Au_sc-ISAAC-10000--2023-04-06T160750

####################################################################################################

using XLSX
using DataFrames

# Define directory path
dir_path = "path/to/directory"

cd(dir_path)

dir_folders=readdir(dir_path,join=true)

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


# df1.+df2
# isfile