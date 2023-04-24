# Define directory path
dir_path = expanduser("~/OneDrive - University of Tennessee/graduate/roy research/no-au scattering/results/charge analysis")
# dir_path = expanduser("~/Downloads/charge analysis")
curdir=pwd()

############################################################################################################

using XLSX
using DataFrames
using GLMakie
using Dates

############################################################################################################

function getdffromfile(file::String)
    data = DataFrame(XLSX.readtable(file, 1))
    data.t*=0.01 # convert to ps
    convert.(Float64, data)
end

function makeresultsfolder(desc::String)
    date=Dates.format(now(), "yyyy-mm-ddTHHMMSS-sss")
    path=mkpath("$desc--$date")
    println("Made folder: $path")
    return path
end

############################################################################################################

t=@elapsed begin
    cd(dir_path)
    dir_files=readdir()

    # Define unique folder names
    systypes = ["-noau-", "-oau-"]
    filt_files_sys = [filter(x -> contains(x, name), dir_files) for name in systypes]

    # Define unique file names
    trajtypes=["sc", "tr"]
    filt_files_type = [[filter(x -> contains(x, name), ffs) for name in trajtypes] for ffs in filt_files_sys]

    ts=[300, 400, 500]
    legendlabels=["$t K" for t in ts]
    plot_kwargs = [:markersize => 5, :linewidth => 0.5]

    path=makeresultsfolder("charge_graphs")
    for (systype,filesv) in zip(systypes,filt_files_type), (trajtype,files) in zip(trajtypes,filesv)
        try
            f = Figure(resolution = (1200, 900))
            Axis(f[1,1],xlabel="t (ps)",ylabel="Charge (e⁻)")
            plotobjs=[]
            for file in files
                df=getdffromfile(file)
                obj=scatterlines!(df.t,df.Charge;plot_kwargs...)
                push!(plotobjs,obj)
            end
            Legend(f[1,2],plotobjs,legendlabels)
            save("$path/charge$systype$trajtype.png",f)
        catch e
            cd(curdir)
            error("Error occurred: $e")
        end
    end
    cd(curdir)
end

println("Elapsed time: $t seconds")