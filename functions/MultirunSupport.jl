# support functions for multiple runs

############################################################################################################

function inittrajcontainers()
    h=["T", "Ei", "xcom", "ycom", "xNf", "yNf", "zNf", "xOf", "yOf", "zOf", "vxNf", "vyNf", "vzNf", "vxOf", "vyOf", "vzOf", "KEtot", "Etrans", "Erot", "KEvib"]
    trajtrap=DataFrame([name => [] for name in h])
    trajscatter=DataFrame([name => [] for name in h])

    h2=["T", "Ei", "n_scatter", "n_trap", "frac_scatter", "frac_trap"]
    counttraj=DataFrame([name => [] for name in h2])

    return trajtrap, trajscatter, counttraj
end

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

"""
print txt file with summary of the multiple runs.
"""
function outputmultirunsummary(counttraj::DataFrame,runpath=".")
    randomtraj=all(ismissing,no.xi) || all(ismissing,no.yi)
    vary=""
    if count(!ismissing,param.T)>1; vary*="T, "; end
    if count(!ismissing,no.Et_i)>1; vary*="Ei, "; end
    vary*="xy"
    
    file="$runpath/summary.txt"
    open(file,"w") do io
        println(io,"Summary of multiple runs")
        println(io)
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
print txt file with summary of the multiple runs.
"""
function outputtrajinfo(trajscatter::DataFrame,trajtrap::DataFrame,runpath=".")
    name_sc="$runpath/traj_scattered.xlsx"
    name_tr="$runpath/traj_trapped.xlsx"
    write_xlsx(name_sc,trajscatter)
    write_xlsx(name_tr,trajtrap)
end

############################################################################################################

"""
run multiple no/au trajectories and output run info to results folder
"""
function runMultiNOAuTrajectory(path::String=makeresultsfolder(noaurundesc,"normalrun");fixorient=false)
	trajtrap, trajscatter, counttraj = inittrajcontainers()

    for temp in skipmissing(param.T)
        tempt=Int64(ustrip(u"K",temp))
        global aurundesc="Au_slab-T_$tempt"
        runAuSlabEquilibration(temp)
        tpath=mkpath("$path/T $tempt")

        for ei in skipmissing(no.Et_i)
            eit=round(ustrip(u"e_MD",ei);digits=2)
            eipath=mkpath("$tpath/Ei $eit")

            nscatter=0
            ntrap=0
            if fixorient && tempt==Torient
                for θ in skipmissing(no.θorient)
                    θt=Int64(ustrip(u"°",θ))
                    θpath=mkpath("$eipath/orient $θt")

                    if all(ismissing, no.xi) || all(ismissing, no.yi)
                        # run random trajectories
                        xs=[au.aPBCx[1]*rand() for _ in Base.OneTo(acttraj)]
                        ys=[au.aPBCy[1]*rand() for _ in Base.OneTo(acttraj)]
    
                        for (x,y) in zip(xs,ys)
                            xt=round(ustrip(u"Å",x);digits=3)
                            yt=round(ustrip(u"Å",y);digits=3)
                            xypath=makeresultsfolder("$θpath/x $xt y $yt")
            
                            sys=runNOAuTrajectory(x,y,temp,ei,xypath;fixorient=fixorient)
    
                            finalNOinfo=finalE_molec(sys)
                            allinfo=vcat([tempt, eit, xt, yt], finalNOinfo)
                            if checkscattering(sys)
                                nscatter+=1
                                push!(trajscatter, allinfo)
                            else
                                ntrap+=1
                                push!(trajtrap, allinfo)
                            end
                            if debug; break; end
                        end
                    else
                        # run specified trajectories
                        for (x,y) in zip(skipmissing(no.xi), skipmissing(no.yi))
                            xt=round(ustrip(u"Å",x);digits=3)
                            yt=round(ustrip(u"Å",y);digits=3)
                            xypath=makeresultsfolder("$θpath/x $xt y $yt")
            
                            sys=runNOAuTrajectory(x,y,temp,ei,xypath)
    
                            finalNOinfo=finalE_molec(sys)
                            allinfo=vcat([tempt, eit, xt, yt], finalNOinfo)
                            if checkscattering(sys)
                                nscatter+=1
                                push!(trajscatter, allinfo)
                            else
                                ntrap+=1
                                push!(trajtrap, allinfo)
                            end
                            if debug; break; end
                        end
                    end
                end
            else
                if all(ismissing, no.xi) || all(ismissing, no.yi)
                    # run random trajectories
                    xs=[au.aPBCx[1]*rand() for _ in Base.OneTo(acttraj)]
                    ys=[au.aPBCy[1]*rand() for _ in Base.OneTo(acttraj)]

                    for (x,y) in zip(xs,ys)
                        xt=round(ustrip(u"Å",x);digits=3)
                        yt=round(ustrip(u"Å",y);digits=3)
                        xypath=makeresultsfolder("$eipath/x $xt y $yt")
        
                        sys=runNOAuTrajectory(x,y,temp,ei,xypath;fixorient=fixorient)

                        finalNOinfo=finalE_molec(sys)
                        allinfo=vcat([tempt, eit, xt, yt], finalNOinfo)
                        if checkscattering(sys)
                            nscatter+=1
                            push!(trajscatter, allinfo)
                        else
                            ntrap+=1
                            push!(trajtrap, allinfo)
                        end
                        if debug; break; end
                    end
                else
                    # run specified trajectories
                    for (x,y) in zip(skipmissing(no.xi), skipmissing(no.yi))
                        xt=round(ustrip(u"Å",x);digits=3)
                        yt=round(ustrip(u"Å",y);digits=3)
                        xypath=makeresultsfolder("$eipath/x $xt y $yt")
        
                        sys=runNOAuTrajectory(x,y,temp,ei,xypath)

                        finalNOinfo=finalE_molec(sys)
                        allinfo=vcat([tempt, eit, xt, yt], finalNOinfo)
                        if checkscattering(sys)
                            nscatter+=1
                            push!(trajscatter, allinfo)
                        else
                            ntrap+=1
                            push!(trajtrap, allinfo)
                        end
                        if debug; break; end
                    end
                end
            end

            totaltraj=ntrap+nscatter
            fracscatter=nscatter/totaltraj
            fractrap=ntrap/totaltraj
            push!(counttraj, [tempt, eit, nscatter, ntrap, fracscatter, fractrap])
            if debug; break; end
        end
        if debug; break; end
    end

    # zip all folders
    zipfolders(path)
    outputmultirunsummary(counttraj,path)
    outputtrajinfo(trajscatter,trajtrap,path)
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
