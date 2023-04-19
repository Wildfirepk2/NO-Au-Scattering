# functions for setting up NO/Au system in Molly.

############################################################################################################

"""
get equilibrated au coords from previous run
"""
function getEquilAuCoords()
    audir=getAuDirPath("results")
    if "syscoords.xlsx" in readdir(audir)
       coordsfile="$audir/syscoords.xlsx"
       xfcoord=XLSX.readxlsx(coordsfile)
       sheets=XLSX.sheetnames(xfcoord)
       sheetlastcoord=sheets[end]
       dfcoord=DataFrame(XLSX.readtable(coordsfile,sheetlastcoord))
       [SA[dfcoord[i,1],dfcoord[i,2],dfcoord[i,3]]u"d_MD" for i in 1:nrow(dfcoord)]
    else
       readcoorddir=readdir("$audir/syscoords") # incorrectly sorts by number. eg 1,10,11,12,...,2,etc
       sortedcoorddir = sort(readcoorddir, by = x -> parse(Int, match(r"\d+", x).match)) # sorted by name. eg syscoords 1,2,3,4,etc
       lastcoords=sortedcoorddir[end]
       dfcoord=CSV.read("$audir/syscoords/$lastcoords",DataFrame)
       [SA[dfcoord[i,1],dfcoord[i,2],dfcoord[i,3]]u"d_MD" for i in 1:nrow(dfcoord)]
    end
end

############################################################################################################

"""
get last au velocities from previous run
"""
function getLastAuVelocities()
    audir=getAuDirPath("results")
    velfile="$audir/sysvelocities.xlsx"
    xfvel=XLSX.readxlsx(velfile)
    sheets=XLSX.sheetnames(xfvel)
    sheetlastcoord=sheets[end]
    dfvel=DataFrame(XLSX.readtable(velfile,sheetlastcoord))
    [SA[dfvel[i,1],dfvel[i,2],dfvel[i,3]]u"v_MD" for i in 1:nrow(dfvel)]
end

############################################################################################################

"""initial NO bond length/vib velocity"""
function getrvNO(fixvib::Bool=false)
    idx=fixvib ? no.vib_phase[1]+1 : rand(1:length(no.r))
    r=no.r[idx]
    v=no.v[idx]
    return r,v
end

############################################################################################################

"""random initial NO orientation. output unit vector from COM to N"""
function getNOorient()
    vec=2rand(3).-1
    normalize(vec)
end

"""
specified initial NO orientation

θ: angle between NO bond and surface normal (0°: N pointing down)
"""
function getNOorient(θ::Unitful.DimensionlessQuantity)
    [sin(θ),0,cos(θ)]
end

############################################################################################################

function initNOCoords(rNOi,u::Vector,xcom=au.aPBCx[1]*rand(),ycom=au.aPBCy[1]*rand(),zcom=12u"Å")
    # masses
    mN=no.mN[1]
    mO=no.mO[1]
    μ=mN*mO/(mN+mO)

    zcomact=zcom+maximum(au.z)
    com=SA[xcom,ycom,zcomact]

    rN=com-rNOi*μ/mN*u
    rO=com+rNOi*μ/mO*u
    [rN,rO]
end

############################################################################################################

function initNOAuCoords(rNOi,u::Vector,xcom=au.aPBCx[1]*rand(),ycom=au.aPBCy[1]*rand())
    nocoords=initNOCoords(rNOi,u,xcom,ycom)
    aucoords=getEquilAuCoords()
    vcat(nocoords,aucoords)
end

############################################################################################################

function initNOAuAtoms()
    noatoms=[Atom(index=1, mass=no.mN[1]), Atom(index=2, mass=no.mO[1])]
    auatoms=[Atom(index=i, mass=au.m[1]) for i in 3:au.N[1]+2]
    vcat(noatoms,auatoms)
end

############################################################################################################

function initNOVelocities(vrel,u::Vector,Ei=no.Et_i[1])
    # masses
    mN=no.mN[1]
    mO=no.mO[1]
    mt=mN+mO
    μ=mN*mO/mt

    # convert incident molecule energy to velocity
    θ=no.θi[1]
    mass=mt*N_A # now in kg/mol
    vmag=sqrt(2*Ei/mass)
    vcom=SA[vmag*sin(θ),0u"m/s",-vmag*cos(θ)] # negative: pointing down

    vN=vcom-μ/mN*vrel*u
    vO=vcom+μ/mO*vrel*u
    [vN,vO]
end

############################################################################################################

function initNOAuVelocities(vrel,u::Vector,Ei=no.Et_i[1])
    vno=initNOVelocities(vrel,u,Ei)
    auv=getLastAuVelocities()
    vcat(vno,auv)
end

############################################################################################################
#\fix. may remove later
"""
helper function: check if au folder exists in directory
"""
function checkAuDir(path::String="results")
    resultsinpath=readdir(path;join=true)
    i_au=findfirst(contains.(resultsinpath,aurundesc))
    !(i_au === nothing) # if nothing found, !true=false
end

############################################################################################################

"""
helper function: get path of au folder in directory. if not found, return nothing
"""
function getAuDirPath(path::String=".")
    resultsinpath=readdir(path;join=true)
    i_au=findfirst(contains.(resultsinpath,aurundesc))
    i_au isa Nothing ? i_au : resultsinpath[i_au]
end

############################################################################################################

"""
initiallize NO/Au system. NO placed at specified x,y position
"""
function initNOAuSys(xcom::Unitful.Length=au.aPBCx[1]*rand(),
                       ycom::Unitful.Length=au.aPBCy[1]*rand(),
                       Ei::EnergyPerMole=no.Et_i[1];
                       θorient=nothing,
                       fixvib::Bool=false,)
    # initial NO BL/vib vel
    rNOi,vrel=getrvNO(fixvib)

    # initial NO orientation
    u = θorient===nothing ? getNOorient() : getNOorient(θorient)

    # Molly sys params
    atoms=initNOAuAtoms()
    pairwise_inters=(NOAuInteraction(true),) # using neighbor list=true
    coords=initNOAuCoords(rNOi,u,xcom,ycom)
    velocities=initNOAuVelocities(vrel,u,Ei) # NO vel based on √(2Ei/m). Au slab set to last vels
    boundary=simboxdims # periodic in x,y
    neighbor_finder=NONeighborFinder(
                                    N_cutoff=PES_GS.AuNcutoff[1],
                                    O_cutoff=PES_GS.AuOcutoff[1],
                                    )
    if simplerun                                  
       loggers=(
                et=TotalEnergyLogger(steps_dyn), # checking energy conservation AT LAST TIME STEP ONLY
                pe=PotentialEnergyLogger(steps_dyn),
                ke=KineticEnergyLogger(steps_dyn),
                velocities=VelocityLogger(steps_dyn),
                charge=ChargeLogger(actsteplog),
                coords=CoordinateLogger(actsteplog),
                )
    else
       loggers=(
          et=TotalEnergyLogger(actsteplog), # checking energy conservation 
          pe=PotentialEnergyLogger(actsteplog),
          ke=KineticEnergyLogger(actsteplog),
          velocities=VelocityLogger(steps_dyn),
          forces=ForceLogger(steps_dyn), # E needs to be calculated before F
          charge=ChargeLogger(actsteplog),
          coords=CoordinateLogger(actsteplog),
       )
    end
    dt=param.dt[1]
    remove_CM_motion=false # dont remove center of mass motion to keep layer fixed. may revert.

    # defining system/simulator
    s = System(
             atoms=atoms,
             pairwise_inters=pairwise_inters,
             coords=coords,
             velocities=velocities,
             boundary=boundary,
             neighbor_finder=neighbor_finder,
             loggers=loggers,
             )
    simul = VelocityVerlet(
                             dt=dt,
                             remove_CM_motion=remove_CM_motion,
                          )

    return s, simul
end

############################################################################################################

"""
run no/au trajectory and output run info to results folder
"""
function runNOAuTrajectory(xcom::Unitful.Length=au.aPBCx[1]*rand(),
                             ycom::Unitful.Length=au.aPBCy[1]*rand(),
                             T::Unitful.Temperature=param.T[1],
                             Ei::EnergyPerMole=no.Et_i[1],
                             path::String=makeresultsfolder(noaurundesc,steps_dyn);
                             θorient=nothing,
                             fixvib=false,
                             )
    sys_NOAu, simulator_NOAu = initNOAuSys(xcom, ycom, Ei; θorient=θorient, fixvib=fixvib)
    println("Conditions: T=$T, Ei=$Ei, θorient=$θorient, xcom=$xcom, ycom=$ycom")

    t=@elapsed runMDprintresults(sys_NOAu, noaurundesc, simulator_NOAu, steps_dyn, path, T, Ei)
    checkEconserved(sys_NOAu)

    println("NO/Au trajectory is complete")
    println("Time to run: $t seconds")
    println()

    return sys_NOAu
end
