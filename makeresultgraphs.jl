# dir_path = expanduser("~/OneDrive - University of Tennessee/graduate/roy research/no-au scattering/results/final excels/")
dir_path = expanduser("~/Downloads/final_excels")

###########################################################################################

using XLSX
using DataFrames
using Unitful
using GLMakie
using Dates
include("functions/MDUnits.jl")

###########################################################################################

function initdffromxlsx(filename::String)
    df = DataFrame(XLSX.readtable(filename, 1))
    sort!(df, 1) # to connect the dots in the plot in order
    # \fix: put dict in fn

    # if path contains "_exp", use expdict. else use mydict
    unit_dict = contains(filename, "_exp") ? expunits : mdunits
    
    for col_name in names(df)
        # if col name is in mydict, convert to unit. else, yank the col
        if haskey(unit_dict, col_name)
            df[!, col_name] *= unit_dict[col_name]
            df[!, col_name] = df[!, col_name] .|> expunits[col_name] # convert units to expdict
        else
            select!(df, Not(col_name))
        end
    end

    return df
end

function outputgraph(mydf::DataFrame, expdf::DataFrame, xvar::String, comp::String, sd_comp::String, path::String=makeresultsfolder("roy plots"))
    systype=dfnamesdict[mydf]
    if "T" in names(mydf)
        # tr or Er: first filter mydf for rows where T == Texp
        Tcomp=unique(expdf[!, "T"])
        mydft=filter(row -> all(row["T"] == Tcomp[1]), mydf)
        plot_kwargs = []
    else
        # charge: as is
        mydft=copy(mydf)
        plot_kwargs = [:markersize => 5, :linewidth => 0.5]
    end

    dfvec=[mydft, expdf]
    legendlabels=["MD", "Literature"]
    xlabel=colnametolabel[xvar]
    ylabel=colnametolabel[comp]
    f = Figure(resolution = (900, 600)) # default 800x600
    Axis(f[1, 1], xlabel = xlabel, ylabel = ylabel)
    plotobjs=[]
    for df in dfvec
        xs=ustrip.(df[!,xvar])
        ys=ustrip.(df[!,comp])
        y_err=ustrip.(df[!,sd_comp])
        obj=scatterlines!(xs,ys;plot_kwargs...)
        errorbars!(xs, ys, y_err,
                    whiskerwidth = 10,
                    )
        push!(plotobjs,obj)
    end
    Legend(f[1, 2],plotobjs,legendlabels)
    save("$path/$systype-$comp v $xvar-roy comp.png", f)
    return f
end

function outputgraph(df::DataFrame, filt::String, xvar::String, yvar::String, path::String=makeresultsfolder("my plots"))
    systype=dfnamesdict[df]
    iters=unique(df[!, filt])
    xlabel=colnametolabel[xvar]
    ylabel=colnametolabel[yvar]
    legendlabels=string.(round.(Int64,unit(iters[1]),iters))

    f = Figure(resolution = (900, 600)) # default 800x600
    Axis(f[1, 1], xlabel = xlabel, ylabel = ylabel)
    plotobjs=[]
    for iter in iters
        # filter df for rows where vals in the $filter col equal $iter
        dft = filter(row -> all(row[filt] == iter), df)
        xs=ustrip.(dft[!,xvar])
        ys=ustrip.(dft[!,yvar])
        obj=scatterlines!(xs,ys)
        push!(plotobjs,obj)
    end
    Legend(f[1, 2],plotobjs,legendlabels)
    
    # save plot to $path
    save("$path/$systype-$yvar v $xvar-$filt series.png", f)
    return f
end

function makeresultsfolder(desc::String)
    date=Dates.format(now(), "yyyy-mm-ddTHHMMSS-sss")
    path=mkpath("$desc--$date")
    println("Made folder: $path")
    return path
end

###########################################################################################

# make dict for T, Ei, and θorient, avg_Erot_sc mapping to their respective units
mdunits = Dict(
    "T" => u"K",
    "Ei" => u"e_MD",
    "θorient" => u"°",
    "avg_Erot_sc" => u"e_MD",
    "avg_Etransfer_sc" => u"e_MD",
    "t" => u"t_MD",
    "Charge" => u"e⁻",
    "frac_trap" => unit(1), # hack to assign DimensionlessQuantity \fix
)

expunits = Dict(
    "T" => u"K",
    "Ei" => u"kJ/mol",
    "θorient" => u"°",
    "avg_Erot_sc" => u"kJ/mol",
    "avg_Etransfer_sc" => u"kJ/mol",
    "t" => u"ps",
    "Charge" => u"e⁻",
    "frac_trap" => unit(1), # hack to assign DimensionlessQuantity \fix
)

# colnametolabel = Dict(
#     "T" => L"T\ \(%$(expunits[String(:T)]))", 
#     "Ei" => L"E_i\ \(%$(expunits[String(:Ei)]))",
#     "θorient" => L"θ\ \(%$(expunits[String(:θorient)]))",
#     "avg_Erot_sc" => L"\langle E_{rot} \rangle\ (%$(expunits[String(:avg_Erot_sc)]))",#\fix use latex: L"x^{2}"
#     "avg_Etransfer_sc" => L"\langle E_{transfer} \rangle\ \(%$(expunits[String(:avg_Etransfer_sc)]))",
#     "t" => L"t\ \(%$(expunits[String(:t)]))",
#     "Charge" => L"Charge\ \(%$(expunits[String(:Charge)]))",
#     "frac_trap" => L"Trapping\ \Probability",
# )

colnametolabel = Dict(
    "T" => "T ($(expunits["T"]))", 
    "Ei" => "Eᵢ ($(expunits["Ei"]))",
    "θorient" => "θ ($(expunits["θorient"]))",
    "avg_Erot_sc" => "⟨Eᵣ⟩ ($(expunits["avg_Erot_sc"]))",#\fix use latex: L"x^{2}"
    "avg_Etransfer_sc" => "⟨Eₜ⟩ ($(expunits["avg_Etransfer_sc"]))",
    "t" => "t ($(expunits["t"]))",
    "Charge" => "Charge ($(expunits["Charge"]))",
    "frac_trap" => "Trapping Probability",
)

###########################################################################################

t=@elapsed begin
    curdir=pwd()
    cd(dir_path)
    # dir_folders=readdir()

    # read all excel files into dataframes
    noau_norm_df = initdffromxlsx("summary_noau_normalrun_final.xlsx")
    noau_fix_df = initdffromxlsx("summary_noau_fixorient_final.xlsx")
    oau_df = initdffromxlsx("summary_oau_final.xlsx")
    charge_df = initdffromxlsx("charge_final.xlsx")
    noau_norm_exp_df = initdffromxlsx("noau_normalrun_exp.xlsx")
    charge_exp_df = initdffromxlsx("charge_exp.xlsx")

    dfvec=[noau_norm_df, noau_fix_df, oau_df, charge_df, noau_norm_exp_df, charge_exp_df]
    dfvecnames=["noau_norm_df", "noau_fix_df", "oau_df", "charge_df", "noau_norm_exp_df", "charge_exp_df"]
    dfnamesdict=Dict(zip(dfvec,dfvecnames))

    # comp to roy
    dfvec_roy=[noau_norm_df,charge_df]
    dfveccomp_roy=[noau_norm_exp_df,charge_exp_df]
    compares=["avg_Erot_sc", "frac_trap", "Charge"]
    sd_compares=["Erot_SD", "trap_SD", ""]
    xvars=["Ei", "t"]
    path=makeresultsfolder("roy plots")
    for df in dfvec_roy, dfcomp in dfveccomp_roy, (compare,sd_comp) in zip(compares,sd_compares), xvar in xvars
        if all([xvar in names(df), compare in names(df), xvar in names(dfcomp), compare in names(dfcomp)])
            outputgraph(df, dfcomp, xvar, compare, sd_comp, path)
        end
        # break
    end

    # my studies
    dfvec_me=[noau_norm_df, noau_fix_df, oau_df]
    filters=["T", "θorient"]
    xvars=["Ei"]
    yvars=["frac_trap", "avg_Erot_sc", "avg_Etransfer_sc"]

    path=makeresultsfolder("my plots")
    for df in dfvec_me, filt in filters, xvar in xvars, yvar in yvars
        if all([filt in names(df), xvar in names(df), yvar in names(df)]) && !(df==noau_fix_df && filt=="T")
            outputgraph(df, filt, xvar, yvar, path)
        end
        # break
    end

    cd(curdir)
end

println("Elapsed time: $t seconds")
