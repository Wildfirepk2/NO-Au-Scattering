# support functions for multiple runs

############################################################################################################

function getacttraj()
    ts=count(!ismissing,param.T)
    eis=count(!ismissing,no.Et_i)
    orients=count(!ismissing,no.θorient)

    basetrajs=param.Ntraj[1]
    fac=(runningnoau+runningoau)*ts*eis + eis*orients
    round(basetrajs/fac)
end

############################################################################################################

"""helper function: zip folders in dir"""
function zipfolders(dirpath::String)
    curdir=pwd()
    cd(dirpath)
    dircontents=readdir()
    archivename="traj_details.zip"
    run(`zip -qr $archivename $dircontents`)

    if isaac || !debug
        for i in eachindex(dircontents)
            rm(dircontents[i],recursive=true)
        end
    end
    cd(curdir)
end

############################################################################################################

function outputmultirunresults(trajscatter::DataFrame,trajtrap::DataFrame,t,runpath=".")
    counttraj=analyzetraj(trajscatter,trajtrap)
    name_sc="traj_scattered"
    name_tr="traj_trapped"
    name_an="traj_summary"
    namevec=[name_sc,name_tr,name_an]
    dfvec=[trajscatter,trajtrap,counttraj]
    outputtrajinfo(namevec,dfvec,runpath)
    outputmultirunsummary(t,trajscatter,trajtrap,counttraj,runpath)
end

############################################################################################################

"""
analyze scattered/trapped trajectories
"""
function analyzetraj(trajscatter::DataFrame,trajtrap::DataFrame)
    # handling of empty dataframes. give each same headers
    if isempty(trajtrap)
        trajtrap=copy(trajscatter)
        empty!(trajtrap)
    elseif isempty(trajscatter)
        trajscatter=copy(trajtrap)
        empty!(trajscatter)
    end

    # choose no/au or o/au
    oau = names(trajscatter)[end]=="Etransfer"

    counttraj=DataFrame()
    iter=[union(trajscatter.T,trajtrap.T),
        union(trajscatter.Ei,trajtrap.Ei)]
    if "θorient" in names(trajtrap)
        push!(iter,union(trajscatter.θorient,trajtrap.θorient))
    end

    # iterate through each combination of T,E,(orient)
    for values in Iterators.product(iter...) # values as tuple
        valheaders=names(trajtrap)[1:length(values)]
        valindf=collect(values)'
        maininfo=DataFrame(valindf,valheaders)

        # filter trajscatter and trajtrap by the current combination of values
        filtered_trajscatter = filter(row -> all(row[valheaders[i]] == valindf[1, i] for i in 1:length(valheaders)), trajscatter)
        filtered_trajtrap = filter(row -> all(row[valheaders[i]] == valindf[1, i] for i in 1:length(valheaders)), trajtrap)

        # data for this T, E, (orient)
        nscatter = nrow(filtered_trajscatter)
        ntrap = nrow(filtered_trajtrap)
        ntotal = nscatter + ntrap
        frac_scatter = nscatter / ntotal
        frac_trap = ntrap / ntotal

        if oau
            avg_Etrans_sc = isempty(filtered_trajscatter) ? NaN : mean(filtered_trajscatter.Etrans)
            avg_Etransfer_sc = isempty(filtered_trajscatter) ? NaN : mean(filtered_trajscatter.Etransfer)
            trackedE=DataFrame(avg_Etrans_sc=avg_Etrans_sc,
                                avg_Etransfer_sc=avg_Etransfer_sc)
        else
            avg_Erot_sc = isempty(filtered_trajscatter) ? NaN : mean(filtered_trajscatter.Erot)
            trackedE=DataFrame(avg_Erot_sc=avg_Erot_sc)
        end

        otherinfo=DataFrame(n_scatter=nscatter,
                            n_trap=ntrap,
                            n_total=ntotal,
                            frac_scatter=frac_scatter,
                            frac_trap=frac_trap)
        allinfo=hcat(maininfo,trackedE,otherinfo)
        append!(counttraj,allinfo)
    end

    return counttraj
end

############################################################################################################

"""
get sys info
"""
function getCPUinfo()
    file="sysinfo.txt"
    open(file,"w") do io
        getCPUinfo(io)
    end
end

function getCPUinfo(io::IO)
    ncpus = Sys.CPU_THREADS
    cpuname = Sys.CPU_NAME
    nmem = (10^-9)Sys.total_memory() 
    nworkers = Sys.Threads.nthreads()

    println(io,"CPU info")
    println(io)
    println(io,"Processor: $cpuname")
    println(io,"Number of cores: $ncpus")
    println(io,"Memory: $nmem GB")
    println(io,"Julia threads: $nworkers")
    println(io)
end

############################################################################################################

"""
print txt file with summary of the multiple runs.
"""
function outputmultirunsummary(t,trajscatter::DataFrame,trajtrap::DataFrame,counttraj::DataFrame=analyzetraj(trajscatter,trajtrap),runpath=".")
    # data for txt file
    daterun=Dates.format(now(), "yyyy-mm-dd HH:MM:SS")
    it=findfirst(x -> x == "ycom", names(trajscatter))
    vary=join(names(trajscatter)[1:it],", ")
    randomtraj=all(ismissing,no.xi) || all(ismissing,no.yi)
    ntraj=nrow(trajscatter)+nrow(trajtrap)
    tpertraj=t/ntraj

    file="$runpath/summary.txt"
    open(file,"w") do io
        println(io,"Summary of multiple runs")
        println(io)
        println(io,"Time of run: $daterun")
        println(io,"Total trajectories: $ntraj")
        println(io,"Variables varied: $vary")
        println(io,"Ran on ISAAC: $isaac")
        println(io,"Random trajectories?: $randomtraj")
        println(io,"Debugging?: $debug")
        println(io,"Short trajectory?: $shortrun")
        println(io,"Simple results?: $simplerun")
        println(io,"Simulation runtime: $t ($tpertraj/trajectory)")
        println(io)
        getCPUinfo(io)
        println(io,"Full trajectory details")
        println(io,counttraj)
    end
end

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

############################################################################################################

"""
run post analysis
"""
function runpostanalysis!(maininfo::DataFrame, sys::System,trajscatter::DataFrame,trajtrap::DataFrame)
    Efinfo=finalE_molec(sys)
    allinfo=hcat(maininfo,Efinfo)
    if checkscattering(sys)
        append!(trajscatter, allinfo)
    else
        append!(trajtrap, allinfo)
    end
end

############################################################################################################

"""
run multiple no/au trajectories and output run info to results folder
"""
function runMultiNOAuTrajectory(;fixorient::Bool=false,
                                    path::String=makeresultsfolder(noaurundesc,fixorient ? "fixorient" : "normalrun")
                                    )
    fixorient ? println("---Running fixed orientation NO/Au scattering---") : println("---Running NO/Au scattering---")
    println()
    trajtrap=DataFrame()
    trajscatter=DataFrame()
    
    randtraj=all(ismissing, no.xi) || all(ismissing, no.yi)
    if randtraj
        xs=[au.aPBCx[1]*rand() for _ in 1:acttraj]
        ys=[au.aPBCy[1]*rand() for _ in 1:acttraj]
    else # specified trajectories
        xs=skipmissing(no.xi)
        ys=skipmissing(no.yi)
    end

    t=@elapsed begin
        if fixorient
            tempt=Int64(ustrip(u"K",Torient))
            global aurundesc="Au_slab-T_$tempt"
            runAuSlabEquilibration(Torient)

            for ei in skipmissing(no.Et_i), orient in skipmissing(no.θorient), (x,y) in zip(xs,ys)
                eit=round(ustrip(u"e_MD",ei);digits=3)
                orientt=Int64(ustrip(u"°",orient))
                xt=round(ustrip(u"Å",x);digits=3)
                yt=round(ustrip(u"Å",y);digits=3)
                resultpath=makeresultsfolder("$path/T $tempt/Ei $eit/orient $orientt/x $xt y $yt")
                sys=runNOAuTrajectory(x,y,Torient,ei,resultpath;θorient=orient)
                df=DataFrame(T=tempt,Ei=eit,θorient=orientt,xcom=xt,ycom=yt)
                runpostanalysis!(df,sys,trajscatter,trajtrap)
                if debug; break; end
            end
        else # rand orient
            for temp in skipmissing(param.T), ei in skipmissing(no.Et_i), (x,y) in zip(xs,ys)
                tempt=Int64(ustrip(u"K",temp))
                global aurundesc="Au_slab-T_$tempt"
                runAuSlabEquilibration(temp)
                
                eit=round(ustrip(u"e_MD",ei);digits=3)
                xt=round(ustrip(u"Å",x);digits=3)
                yt=round(ustrip(u"Å",y);digits=3)
                resultpath=makeresultsfolder("$path/T $tempt/Ei $eit/x $xt y $yt")
                sys=runNOAuTrajectory(x,y,temp,ei,resultpath)
                df=DataFrame(T=tempt,Ei=eit,xcom=xt,ycom=yt)
                runpostanalysis!(df,sys,trajscatter,trajtrap)
                if debug; break; end
            end
        end
    end
    t*=u"s"

    # zip all folders
    zipfolders(path)
    outputmultirunresults(trajscatter,trajtrap,t,path)
    
    println("---NO/Au multiruns complete---")
    println("Total time taken: $t")
    println()

    # return sys # fails
end  

############################################################################################################

"""
run multiple o/au trajectories and output run info to results folder
"""
function runMultiOAuTrajectory(;path::String=makeresultsfolder(oaurundesc,steps_dyn_OAu))
    println("---Running O/Au scattering---")
    println()
    trajtrap=DataFrame()
    trajscatter=DataFrame()
    
    randtraj=all(ismissing, no.xi) || all(ismissing, no.yi)
    if randtraj
        xs=[au.aPBCx[1]*rand() for _ in 1:acttraj]
        ys=[au.aPBCy[1]*rand() for _ in 1:acttraj]
    else # specified trajectories
        xs=skipmissing(no.xi)
        ys=skipmissing(no.yi)
    end

    t=@elapsed begin
        for temp in skipmissing(param.T), ei in skipmissing(no.Et_i), (x,y) in zip(xs,ys)
            tempt=Int64(ustrip(u"K",temp))
            global aurundesc="Au_slab-T_$tempt"
            runAuSlabEquilibration(temp)
            
            eit=round(ustrip(u"e_MD",ei);digits=3)
            xt=round(ustrip(u"Å",x);digits=3)
            yt=round(ustrip(u"Å",y);digits=3)
            resultpath=makeresultsfolder("$path/T $tempt/Ei $eit/x $xt y $yt")
            sys=runOAuTrajectory(x,y,temp,ei,resultpath)
            df=DataFrame(T=tempt,Ei=eit,xcom=xt,ycom=yt)
            runpostanalysis!(df,sys,trajscatter,trajtrap)
            if debug; break; end
        end
    end
    t*=u"s"

    # zip all folders
    zipfolders(path)
    outputmultirunresults(trajscatter,trajtrap,t,path)

    println("---O/Au multiruns complete---")
    println("Total time taken: $t")
    println()

    # return sys # fails
end
