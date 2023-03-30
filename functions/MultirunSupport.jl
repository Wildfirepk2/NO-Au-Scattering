# support functions for multiple runs

function initTExy()
    if debug
        ts=(300:50:350)u"K"
        eis=(20:20:40)u"kJ/mol"
        if randomtraj
            ntraj=2
            xs=[au.aPBCx[1]*rand() for _ in Base.OneTo(ntraj)]
            ys=[au.aPBCy[1]*rand() for _ in Base.OneTo(ntraj)]
        else
            xs=(1:2:3)u"Å"
            ys=xs
        end
    else
        ts=(300:50:300)u"K"
        eis=(erun:25:erun)u"kJ/mol" # \fix
        if randomtraj
            ntraj=400
            xs=[au.aPBCx[1]*rand() for _ in Base.OneTo(ntraj)]
            ys=[au.aPBCy[1]*rand() for _ in Base.OneTo(ntraj)]
        else
            xs=(3:2:21)u"Å"
            ys=xs
        end
    end

    return ts, eis, xs, ys
end

############################################################################################################

function inittrajcontainers()
    h=["T", "Ei", "xNOi", "yNOi", "xNf", "yNf", "zNf", "xOf", "yOf", "zOf", "vxNf", "vyNf", "vzNf", "vxOf", "vyOf", "vzOf", "KEtot", "Etrans", "Erot", "KEvib"]
    trajtrap=DataFrame([name => [] for name in h])
    trajscatter=DataFrame([name => [] for name in h])

    h2=["T", "Ei", "n_scatter", "n_trap", "frac_scatter", "frac_trap"]
    counttraj=DataFrame([name => [] for name in h2])

    return trajtrap, trajscatter, counttraj
end

############################################################################################################

"""
run multiple no/au trajectories and output run info to results folder
"""
function runMultiNOAuTrajectory(ts, eis, xs, ys)
	for i in eachindex(ts)
        param.T[1]=ts[i]
        T=Int64(ustrip(u"K",param.T[1]))
        global aurundesc="Au_slab-T_$T"
        runAuSlabEquilibration()
        tpath=mkpath("$outpath/T $T")

        for j in eachindex(eis)
            no.Et_i[1]=eis[j]
            ei=round(ustrip(u"e_MD",no.Et_i[1]);digits=2)
            eipath=mkpath("$tpath/Ei $ei")

            nscatter=0
            ntrap=0
            for k in eachindex(xs)
                xNOi=xs[k]
                x=round(ustrip(u"Å",xNOi);digits=3)
                yNOi=ys[k]
                y=round(ustrip(u"Å",yNOi);digits=3)
                xypath=makeresultsfolder("$eipath/x $x y $y")

                sys=runNOAuTrajectory(xNOi,yNOi,xypath)
                checkEconserved(sys)

                finalNOinfo=finalE_molec(sys)
                allinfo=vcat([T, ei, x, y], finalNOinfo)
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
    zipfolders(outpath)
    outputmultirunsummary(vary,counttraj,outpath)
    outputtrajinfo(trajscatter,trajtrap,outpath)
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
        tpath=mkpath("$outpath/T $T")

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
                checkEconserved(sys)

                finalOinfo=finalE_molec(sys)
                allinfo=vcat([T, ei, x, y], finalOinfo)
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
    zipfolders(outpath)
    outputmultirunsummary(vary,counttraj,outpath)
    outputtrajinfo(trajscatter,trajtrap,outpath)
end
