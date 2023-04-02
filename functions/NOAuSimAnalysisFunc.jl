# functions to analyze NO/Au Molly simulations after running

############################################################################################################

"""
helper function: get label of atom i. used in outputallatomizcoords
"""
function getatomilabel(sys::System{D, false, T, CU, A, AD, PI} where {D,T,CU,A,AD,PI<:Tuple{NOAuInteraction}},atom_i::Int64)
    if atom_i==1
        label_i="zN"
    elseif atom_i==2
        label_i="zO"
    else
        label_i="zAu$atom_i"
    end
    return label_i
end

############################################################################################################

"""
helper function: title for E graph based on sys
"""
function getEtitle(sys::System{D, false, T, CU, A, AD, PI} where {D,T,CU,A,AD,PI<:Tuple{NOAuInteraction}})
    "NO/Au Scattering: Checking energy conservation"
end

############################################################################################################

"""
helper func: get title for sys
"""
function getsystitle(s::System{D, false, T, CU, A, AD, PI} where {D,T,CU,A,AD,PI<:Tuple{NOAuInteraction}})
    "NO/Au(111) Scattering"
end

############################################################################################################

"""
helper function: get z graph desc based on sys
"""
function getzgraphdesc(s::System{D, false, T, CU, A, AD, PI} where {D,T,CU,A,AD,PI<:Tuple{NOAuInteraction}})
    "NO Height Above Surface with Time"
end

############################################################################################################

"""
print txt file with summary of run.
"""
function outputsummary(sys::System{D, false, T, CU, A, AD, PI} where {D,T,CU,A,AD,PI<:Tuple{NOAuInteraction}},dt,T,Ei,simsteps=NaN,runtime=NaN,runpath=".")
    # description of run
    title="NO/Au(111) Scattering"
    
    # time of run
    daterun=Dates.format(now(), "yyyy-mm-dd HH:MM:SS")

    # number of atoms in system
    natoms=length(sys.atoms)

    # time step of run (dt) in MD units. rounding due to floating point errors
    dtmd=round(u"t_MD",dt;digits=2)

    # total time of run in given/MD units. rounding due to floating point errors
    ttotal=simsteps*dt
    ttotalmd=round(u"t_MD",ttotal;digits=2)

    # mult runs
    eimd=round(u"e_MD",Ei;digits=2)

    # xy position of N/O
    xNi=sys.loggers.coords.history[1][1][1] |> u"Å"
    yNi=sys.loggers.coords.history[1][1][2] |> u"Å"
    zNi=sys.loggers.coords.history[1][1][3] |> u"Å"
    xOi=sys.loggers.coords.history[1][2][1] |> u"Å"
    yOi=sys.loggers.coords.history[1][2][2] |> u"Å"
    zOi=sys.loggers.coords.history[1][2][3] |> u"Å"

    # number of steps between logging points for quantities
    stepsbtwnlogsCoords=sys.loggers.coords.n_steps
    stepsbtwnlogsE=sys.loggers.et.n_steps
    stepsbtwnlogsV=sys.loggers.velocities.n_steps

    # time step between logging points for quantities in given/MD units. rounding due to floating point errors
    logdtCoords=stepsbtwnlogsCoords*dt
    logdtE=stepsbtwnlogsE*dt
    logdtV=stepsbtwnlogsV*dt
    logdtCoordsmd=round(u"t_MD",logdtCoords;digits=2)
    logdtEmd=round(u"t_MD",logdtE;digits=2)
    logdtVmd=round(u"t_MD",logdtV;digits=2)

    if !simplerun
        stepsbtwnlogsF=sys.loggers.forces.n_steps
        logdtF=stepsbtwnlogsF*dt
        logdtFmd=round(u"t_MD",logdtF;digits=2)
    else
        stepsbtwnlogsF=NaN
        logdtF=NaN
        logdtFmd=NaN
    end

    # final slab energy in given/MD units. rounding due to floating point errors
    finalE=sys.loggers.et.history[end]
    finalEmd=round(u"e_MD",finalE;digits=2)

    # write to txt
    file="$runpath/summary.txt"
    open(file,"w") do io
        println(io,title) 
        println(io)
        println(io,"Time of run: $daterun")
        println(io,"Number of atoms: $natoms")
        println(io,"Number of steps: $simsteps")
        println(io,"Time step: $dt ($dtmd)")
        println(io,"Total time: $ttotal ($ttotalmd)")
        println(io,"Simulation runtime: $runtime")
        println(io,"Ran on ISAAC: $isaac")
        println(io)
        println(io,"T: $T")
        println(io,"Incident energy of NO: $Ei ($eimd)")
        println(io,"Initial x position of N: $xNi")
        println(io,"Initial y position of N: $yNi")
        println(io,"Initial z position of N: $zNi")
        println(io,"Initial x position of O: $xOi")
        println(io,"Initial y position of O: $yOi")
        println(io,"Initial z position of O: $zOi")
        println(io)
        println(io,"PESs:")
        println(io,"    neutral_PES_active = $neutral_PES_active")
        println(io,"    ionic_PES_active = $ionic_PES_active")
        println(io,"    coupled_PES_active = $coupled_PES_active")
        println(io)
        println(io,"Steps between logging quantities")
        println(io,"    Coords: $stepsbtwnlogsCoords")
        println(io,"    Energies: $stepsbtwnlogsE")
        println(io,"    Forces: $stepsbtwnlogsF")
        println(io,"    Velocities: $stepsbtwnlogsV")
        println(io)
        println(io,"Time period between logging quantities")
        println(io,"    Coords: $logdtCoords ($logdtCoordsmd)")
        println(io,"    Energies: $logdtE ($logdtEmd)")
        println(io,"    Forces: $logdtF ($logdtFmd)")
        println(io,"    Velocities: $logdtV ($logdtVmd)")
        println(io)
        println(io,"Final System Total Energy: $finalE ($finalEmd)")
    end
end

############################################################################################################

"""
output animation of trajectory for system
"""
function outputanimation(sys::System{D, false, T, CU, A, AD, PI} where {D,T,CU,A,AD,PI<:Tuple{NOAuInteraction}},path=".")
    connections=[(1,2)]
    connection_color=[:purple]

    # colors of all atoms
    NOcolors=[:blue,:red]
    Aucolors=[:gold for _ in 1:au.N[1]]
    syscolors=vcat(NOcolors,Aucolors)

    # output animation of simulation in $path
    if isaac
        myvisualize(sys.loggers.coords, sys.boundary, "$path/animation.mp4";show_boundary=false,connections=connections,connection_color=connection_color,color=syscolors,)
    else
        visualize(sys.loggers.coords, sys.boundary, "$path/animation.mp4";show_boundary=false,connections=connections,connection_color=connection_color,color=syscolors,)
    end
end

############################################################################################################

"""
output all system info. animation. logger quantities

run description needs to contain either "Au slab" or "NO/Au"
"""
function outputsysinfo(sys::System{D, false, T, CU, A, AD, PI} where {D,T,CU,A,AD,PI<:Tuple{NOAuInteraction}},dt,path::String=".")
    # output quantities to excel file in separate folder
    if isaac
        outputallsyscoords(sys,dt,path,false,false)
    else
        outputallsyscoords(sys,dt,path)
    end
    outputallatomizcoords(sys,dt,1:2,path)
    if !simplerun
        outputanimation(sys,path)
        outputsysE(sys,dt,path)
        outputallsysforces(sys,dt,path)
        outputallsysvelocities(sys,dt,path)
    end
end

############################################################################################################

function checkscattering(s::System{D, false, T, CU, A, AD, PI} where {D,T,CU,A,AD,PI<:Tuple{NOAuInteraction}})
    mN=s.atoms[1].mass
    mO=s.atoms[2].mass
    rNf=s.loggers.coords.history[end][1]
    rOf=s.loggers.coords.history[end][2]

    cutoff=10u"Å"
    zcom_f=getzcom(mN,mO,rNf,rOf)

    zcom_f>=cutoff ? true : false
end

############################################################################################################

"""
get final NO coords/vel + energies. output as vec in MD units
"""
function finalE_molec(s::System{D, false, T, CU, A, AD, PI} where {D,T,CU,A,AD,PI<:Tuple{NOAuInteraction}})
    # NO masses
    mN=s.atoms[1].mass
    mO=s.atoms[2].mass
    mt=mN+mO
    μNO=mN*mO/mt

    # NO last coords
    finalrN=s.loggers.coords.history[end][1] .|> u"d_MD"
    finalrO=s.loggers.coords.history[end][2] .|> u"d_MD"
    rvec_NO=vector(finalrN,finalrO,s.boundary)
    ur_NO=normalize(rvec_NO)

    # NO last velocities
    finalvN=s.loggers.velocities.history[end][1] .|> u"v_MD"
    finalvO=s.loggers.velocities.history[end][2] .|> u"v_MD"
    vvec_NO=finalvO-finalvN

    # NO final translational energy
    vcom_sq=sum([getCOM(mN,mO,finalvN[i],finalvO[i])^2 for i in eachindex(finalvN)])
    finalEtrans=0.5*mt*vcom_sq*N_A |> u"e_MD"

    # NO final vibrational KE
    vvib_sq=dot(ur_NO,vvec_NO)^2
    finalKEvib=0.5*μNO*vvib_sq*N_A |> u"e_MD"

    # NO final total KE/rotational energy
    finalKEtot=0.5*sum(mN.*finalvN.^2+mO.*finalvO.^2)*N_A |> u"e_MD"
    finalErot=finalKEtot-finalEtrans-finalKEvib

    # output final NO coords, velocities, energies
    df=DataFrame(xNf=finalrN[1],
            yNf=finalrN[2],
            zNf=finalrN[3],
            xOf=finalrO[1],
            yOf=finalrO[2],
            zOf=finalrO[3],
            vxNf=finalvN[1],
            vyNf=finalvN[2],
            vzNf=finalvN[3],
            vxOf=finalvO[1],
            vyOf=finalvO[2],
            vzOf=finalvO[3],
            KEtot=finalKEtot,
            Etrans=finalEtrans,
            Erot=finalErot,
            KEvib=finalKEvib)
    ustrip.(df)
end
