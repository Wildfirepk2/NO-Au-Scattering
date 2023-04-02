# support functions for multiple runs

############################################################################################################

function getacttraj()
    traj=param.Ntraj[1]
    lenT=count(!ismissing,param.T)
    lenEi=count(!ismissing,no.Et_i)
    # all(ismissing, no.xi) || all(ismissing, no.yi) ? lenxy=0 : lenxy=count(!ismissing,no.xi)
    all(ismissing, no.θorient) ? lenorient=0 : lenorient=count(!ismissing,no.θorient)

    fac=lenT*lenEi+lenorient*lenEi
    runningoau ? fac+=lenT*lenEi : fac+=0
    round(Int64,traj/fac)
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

function outputmultirunresults(trajscatter::DataFrame,trajtrap::DataFrame,runpath=".")
    outputmultirunsummary(trajscatter,trajtrap,runpath)
    outputtrajinfo(trajscatter,trajtrap,runpath)
end

############################################################################################################

"""
print txt file with summary of the multiple runs.
"""
function outputmultirunsummary(trajscatter::DataFrame,trajtrap::DataFrame,runpath=".")
    # time of run
    daterun=Dates.format(now(), "yyyy-mm-dd HH:MM:SS")

    randomtraj=all(ismissing,no.xi) || all(ismissing,no.yi)
    vary=""
    if count(!ismissing,param.T)>1; vary*="T, "; end
    if count(!ismissing,no.Et_i)>1; vary*="Ei, "; end
    vary*="xy"

    counttraj=DataFrame()
    if isempty(trajtrap)
        ts=unique(trajscatter.T)
        eis=unique(trajscatter.Ei)
    else
        ts=union(trajscatter.T,trajtrap.T)
        eis=union(trajscatter.Ei,trajtrap.Ei)
    end

    for T in ts, Ei in eis
        filtered_trajscatter = filter(row -> row.T == T && row.Ei == Ei, trajscatter)
        filtered_trajtrap = filter(row -> row.T == T && row.Ei == Ei, trajtrap)

        # data for this T and Ei
        avg_Erot_sc=mean(filtered_trajscatter.Erot)
        nscatter=nrow(filtered_trajscatter)
        ntrap=nrow(filtered_trajtrap)
        ntotal=nscatter+ntrap
        frac_scatter=nscatter/ntotal
        frac_trap=ntrap/ntotal

        df=DataFrame(T=T,
                    Ei=Ei,
                    avg_Erot_sc=avg_Erot_sc,
                    n_scatter=nscatter,
                    n_trap=ntrap,
                    n_total=ntotal,
                    frac_scatter=frac_scatter,
                    frac_trap=frac_trap)
        append!(counttraj,df)
    end
    
    file="$runpath/summary.txt"
    open(file,"w") do io
        println(io,"Summary of multiple runs")
        println(io)
        println(io,"Time of run: $daterun")
        println(io,"Variables varied: $vary")
        println(io,"Ran on ISAAC: $isaac")
        println(io,"Random Trajectories?: $randomtraj")
        println(io,"Debugging?: $debug")
        println(io)
        println(io,counttraj)
    end
end

############################################################################################################

"""
print excel files with detailed trajectory info of the multiple runs.
"""
function outputtrajinfo(trajscatter::DataFrame,trajtrap::DataFrame,runpath=".")
    if !isempty(trajscatter)
        name_sc="$runpath/traj_scattered.xlsx"
        write_xlsx(name_sc,trajscatter)
    end
    if !isempty(trajtrap)
        name_tr="$runpath/traj_trapped.xlsx"
        write_xlsx(name_tr,trajtrap)
    end
end

############################################################################################################

"""
run post analysis
"""
function runpostanalysis!(tempt, eit, xt, yt, sys::System,trajscatter::DataFrame,trajtrap::DataFrame)
    maininfo=DataFrame(T=tempt, Ei=eit, xcom=xt, ycom=yt)
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
    trajtrap=DataFrame()
    trajscatter=DataFrame()
    
    randtraj=all(ismissing, no.xi) || all(ismissing, no.yi)
    if randtraj
        xs=[au.aPBCx[1]*rand() for _ in Base.OneTo(acttraj)]
        ys=[au.aPBCy[1]*rand() for _ in Base.OneTo(acttraj)]
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
                eit=round(ustrip(u"e_MD",ei);digits=2)
                orientt=Int64(ustrip(u"°",orient))
                xt=round(ustrip(u"Å",x);digits=3)
                yt=round(ustrip(u"Å",y);digits=3)
                resultpath=makeresultsfolder("$path/T $tempt/Ei $eit/orient $orientt/x $xt y $yt")
                sys=runNOAuTrajectory(orient,x,y,Torient,ei,resultpath)
                runpostanalysis!(tempt, eit, xt, yt, sys,trajscatter,trajtrap)
                if debug; break; end
            end
        else # rand orient
            for temp in skipmissing(param.T), ei in skipmissing(no.Et_i), (x,y) in zip(xs,ys)
                tempt=Int64(ustrip(u"K",temp))
                global aurundesc="Au_slab-T_$tempt"
                runAuSlabEquilibration(temp)
                
                eit=round(ustrip(u"e_MD",ei);digits=2)
                xt=round(ustrip(u"Å",x);digits=3)
                yt=round(ustrip(u"Å",y);digits=3)
                resultpath=makeresultsfolder("$path/T $tempt/Ei $eit/x $xt y $yt")
                sys=runNOAuTrajectory(x,y,temp,ei,resultpath)
                runpostanalysis!(tempt, eit, xt, yt, sys,trajscatter,trajtrap)
                if debug; break; end
            end
        end

        # zip all folders
        zipfolders(path)
        outputmultirunresults(trajscatter,trajtrap,path)
    end
    t*=u"s"
    println("---NO/Au multiruns complete---")
    println("Total time taken: $t")
    println()
end  

############################################################################################################

"""
run multiple o/au trajectories and output run info to results folder
"""
function runMultiOAuTrajectory(ts, eis, xs, ys)
	for i in eachindex(ts)
        param.T[1]=ts[i]
        T=Int64(ustrip(u"K",param.T[1]))
        global aurundesc="Au_slab-T $T"
        runAuSlabEquilibration()
        tpath=mkpath("$path/T $T")

        for j in eachindex(eis)
            no.Et_i[1]=eis[j]
            ei=round(ustrip(u"e_MD",no.Et_i[1]);digits=2)
            eipath=mkpath("$tpath/Ei $ei")

            nscatter=0
            ntrap=0
            for k in eachindex(xs)
                xOi=xs[k]
                x=round(ustrip(u"Å",xOi);digits=3)
                yOi=ys[k]
                y=round(ustrip(u"Å",yOi);digits=3)
                xypath=makeresultsfolder("$eipath/x $x y $y")

                sys=runOAuTrajectory(xOi,yOi,xypath)

                finalOinfo=finalE_molec(sys)
                allinfo=vcat([T, ei], finalOinfo)
                if checkscattering(sys)
                    nscatter+=1
                    push!(trajscatter, allinfo)
                else
                    ntrap+=1
                    push!(trajtrap, allinfo)
                end
            end

            totaltraj=ntrap+nscatter
            fracscatter=nscatter/totaltraj
            fractrap=ntrap/totaltraj
            push!(counttraj, [T, ei, nscatter, ntrap, fracscatter, fractrap])
        end
    end

    # zip all folders
    zipfolders(path)
    outputmultirunsummary(vary,counttraj,path)
    outputtrajinfo(trajscatter,trajtrap,path)
end
