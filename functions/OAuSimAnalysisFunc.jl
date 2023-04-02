# functions to analyze O/Au Molly simulations after running

############################################################################################################

function getatomilabel(sys::System{D, false, T, CU, A, AD, PI} where {D,T,CU,A,AD,PI<:Tuple{OAuInteraction}},atom_i::Int64)
    atom_i==1 ? "zO" : "zAu$atom_i"
end

############################################################################################################

function getEtitle(sys::System{D, false, T, CU, A, AD, PI} where {D,T,CU,A,AD,PI<:Tuple{OAuInteraction}})
    "O/Au Scattering: Checking energy conservation"
end

############################################################################################################

function getsystitle(s::System{D, false, T, CU, A, AD, PI} where {D,T,CU,A,AD,PI<:Tuple{OAuInteraction}})
    "O/Au(111) Scattering"
end

############################################################################################################

"""
helper function: get z graph desc based on sys
"""
function getzgraphdesc(s::System{D, false, T, CU, A, AD, PI} where {D,T,CU,A,AD,PI<:Tuple{OAuInteraction}})
    "O Height Above Surface with Time"
end

############################################################################################################

function outputsummary(sys::System{D, false, T, CU, A, AD, PI} where {D,T,CU,A,AD,PI<:Tuple{OAuInteraction}},dt,T,Ei,simsteps=NaN,runtime=NaN,runpath=".")
    # description of run
    title="O/Au(111) Scattering"
    
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

    # xy position of O
    xOi=sys.loggers.coords.history[1][1][1] |> u"Å"
    yOi=sys.loggers.coords.history[1][1][2] |> u"Å"
    zOi=sys.loggers.coords.history[1][1][3] |> u"Å"

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
        println(io,"Incident energy of O: $Ei ($eimd)")
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

function outputanimation(sys::System{D, false, T, CU, A, AD, PI} where {D,T,CU,A,AD,PI<:Tuple{OAuInteraction}},path=".")
    # colors of all atoms
    Ocolor=[:red]
    Aucolors=[:gold for _ in 1:au.N[1]]
    syscolors=vcat(Ocolor,Aucolors)
 
    # output animation of simulation in $path
    if isaac
        myvisualize(sys.loggers.coords, sys.boundary, "$path/animation.mp4";show_boundary=false,color=syscolors,)
    else
        visualize(sys.loggers.coords, sys.boundary, "$path/animation.mp4";show_boundary=false,color=syscolors,)
    end
end

############################################################################################################

function outputsysinfo(sys::System{D, false, T, CU, A, AD, PI} where {D,T,CU,A,AD,PI<:Tuple{OAuInteraction}},dt,path::String=".")
    # output quantities to excel file in separate folder
    if isaac
        outputallsyscoords(sys,dt,path,false,false)
    else
        outputallsyscoords(sys,dt,path)
    end
    outputallatomizcoords(sys,dt,1,path)
    if !simplerun
        outputanimation(sys,path)
        outputsysE(sys,dt,path)
        outputallsysforces(sys,dt,path)
        outputallsysvelocities(sys,dt,path)
    end
end

############################################################################################################

function checkscattering(s::System{D, false, T, CU, A, AD, PI} where {D,T,CU,A,AD,PI<:Tuple{OAuInteraction}})
    zOf=s.loggers.coords.history[end][1][3]-maximum(au.z)
    cutoff=10u"Å"

    zOf>=cutoff ? true : false
end

############################################################################################################

function finalE_molec(s::System{D, false, T, CU, A, AD, PI} where {D,T,CU,A,AD,PI<:Tuple{OAuInteraction}})
    # O mass
    mO=s.atoms[1].mass
 
    # O last coords
    finalrO=s.loggers.coords.history[end][1] .|> u"d_MD"
 
    # O last velocities
    finalvO=s.loggers.velocities.history[end][1] .|> u"v_MD"
 
    # O final translational energy
    finalEtrans=0.5*mO*sum(finalvO.^2)*N_A |> u"e_MD"
 
    # output final O coord, velocity, energy
    ustrip_vec(vcat(finalrO,finalvO,finalEtrans))
 end