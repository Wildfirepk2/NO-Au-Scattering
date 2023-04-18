
const base_path = expanduser("~/OneDrive - University of Tennessee/graduate/roy research/no-au scattering/results/final excels/")

using XLSX
using DataFrames
using LinearAlgebra
using Unitful
using GLMakie
include("functions/MDUnits.jl")

###########################################################################################

function initdffromxlsx(filename::String)
    path = base_path*filename
    df = DataFrame(XLSX.readtable(path, 1))
    # \fix: put dict in fn

    # if path contains "_exp", use expdict. else use mydict
    unit_dict = contains(filename, "_exp") ? expdict : mydict
    
    for col_name in names(df)
        # if col name is in mydict, convert to unit. else, yank the col
        if haskey(unit_dict, col_name)
            df[!, col_name] *= unit_dict[col_name]
            df[!, col_name] = df[!, col_name] .|> expdict[col_name] # convert units to expdict
        else
            select!(df, Not(col_name))
        end
    end

    return df
end

function outputgraph(df1::DataFrame, df2::DataFrame, var::String)

end

function outputgraph(df1::DataFrame, filter::String, xvar::String, yvar::String)
    # ustrip cols df1 according to expdict
    for col_name in names(df1)
        if haskey(expdict, col_name)
            df1[!, col_name] = ustrip.(df1[!, col_name], expdict[col_name])
        end
    end
    
    iters=unique(df1[!, filter])
    xunit=unit(df1[begin, xvar])
    yunit=unit(df1[begin, yvar])
    if !(xunit isa Unitful.DimensionlessUnits)
        xunit = "($xunit)"
    end
    if !(yunit isa Unitful.DimensionlessUnits)
        yunit = "($yunit)"
    end
    xlabel=colnametolabel[xvar]
    ylabel=colnametolabel[yvar]

    f = Figure(resolution = (900, 600)) # default 800x600
    Axis(f[1, 1], xlabel = "$xlabel $xunit", ylabel = "$ylabel $yunit")
    plotobjs=[]
    for iter in iters
        # filter df1 for rows where vals in filter col equal iter
        df = filter(row -> all(row[filter] == iter), df1)
        obj=scatterlines!(df[!,xvar],df[!,yvar],markersize=5,linewidth=0.5)
        push!(plotobjs,obj)
    end
    if length(plotobjs)>1;Legend(f[1, 2],plotobjs,iters);end
    
    # save plot to $path
    save("$filter-$xvar-$yvar.png", f)
    return f
end

###########################################################################################

# make dict for T, Ei, and θorient, avg_Erot_sc mapping to their respective units
mydict = Dict(
    "T" => u"K",
    "Ei" => u"e_MD",
    "θorient" => u"°",
    "avg_Erot_sc" => u"e_MD",
    "avg_Etransfer_sc" => u"e_MD",
    "t" => u"t_MD",
    "Charge" => u"e⁻",
    "frac_trap" => u"s/s", # hack to assign DimensionlessQuantity \fix
)

expdict = Dict(
    "T" => u"K",
    "Ei" => u"kJ/mol",
    "θorient" => u"°",
    "avg_Erot_sc" => u"kJ/mol",
    "avg_Etransfer_sc" => u"kJ/mol",
    "t" => u"ps",
    "Charge" => u"e⁻",
    "frac_trap" => u"s/s", # hack to assign DimensionlessQuantity \fix
)

colnametolabel = Dict(
    "T" => L"T", 
    "Ei" => L"E_i",
    "θorient" => L"θ",
    "avg_Erot_sc" => L"\langle E_{rot} \rangle",#\fix use latex: L"x^{2}"
    "avg_Etransfer_sc" => L"\langle E_{transfer} \rangle",
    "t" => L"t",
    "Charge" => L"Charge",
    "frac_trap" => L"Trapping Probability",
)

###########################################################################################

# read all excel files into dataframes
noau_norm_df = initdffromxlsx("summary_noau_normalrun_final.xlsx")
noau_fix_df = initdffromxlsx("summary_noau_fixorient_final.xlsx")
oau_df = initdffromxlsx("summary_oau_final.xlsx")
charge_df = initdffromxlsx("charge_final.xlsx")
noau_norm_exp_df = initdffromxlsx("noau_normalrun_exp.xlsx")
charge_exp_df = initdffromxlsx("charge_exp.xlsx")

# comp to roy

# my studies
dfvec=[noau_norm_df, noau_fix_df, oau_df]
filters=["T", "θorient", "T"]
xvars=["Ei", "Ei"]
yvars=["frac_trap", "avg_Erot_sc"]

for (df,filter) in zip(dfvec,filters), (xvar,yvar) in zip(xvars,yvars)
    outputgraph(df, filter, xvar, yvar)
    break
end