# support functions for custom Molly functions

############################################################################################################

@inline @inbounds function getVij_NOAu(i,j,distbtwn,cosθ,dz,dr,boundary)
    En=0u"e_MD"
    Ei=0u"e_MD"
    Ec=0u"e_MD"
    EAu=0u"e_MD"
    if i==1
        # NO interaction + 1 time V_image, WF, Ea calc
        if j==2
            if neutral_PES_active
                En=V00_NO(distbtwn)
            end
            if ionic_PES_active
                Ei=V11_NO(distbtwn)+V11_image(dz)+PES_ionic.ϕ[1]-PES_ionic.Ea[1]
            end
        # N-Au interaction
        else
            if neutral_PES_active
                En=V00_AuN(distbtwn)
            end
            if ionic_PES_active
                Ei=V11_AuN(distbtwn,cosθ)
            end
            if coupled_PES_active
                Ec=V01_AuN(distbtwn)
            end
        end
    elseif i==2
        # O-Au interactions. NOT double counting the NO interaction (j==1)
        if j>1
            if neutral_PES_active
                En=V00_AuO(distbtwn)
            end
            if ionic_PES_active
                Ei=V11_AuO(distbtwn)
            end
            if coupled_PES_active
                Ec=V01_AuO(distbtwn)
            end
        end
    else
        # Au-Au interactions
        it=i-2
        jt=j-2

        # find j in nn[i] and return index. needed to access items in variables based off nn
        idx_j=findfirst(isequal(jt), nn[it])

        # dist of nn pair away from equilibrium. static array
        rij=vector(r0ij[it][:,idx_j],dr,boundary)

        # V nn pair exerts on atom i. number. divide by 2 because of v definition
        EAu=dot(rij,Aijarray[it][idx_j],rij)/2
    end
    return En,Ei,Ec,EAu
end

############################################################################################################

@inline @inbounds function getFij_NOAu(i,j,vec_ij,distbtwn,rNO,uON,u,cosθ,dz,a,b,c,boundary)
    Fn=0u"N/mol"
    Fi=0u"N/mol"
    Fc=0u"N/mol"
    Fimg=0u"N/mol"
    F_AuNc=zeros(3)u"N/mol"
    F_AuAu=zeros(3)u"N/mol"

    if i==1
        # NO interaction
        if j==2
            if neutral_PES_active
                Fn=F00_NO(distbtwn)
            end
            if ionic_PES_active
                Fi=F11_NO(distbtwn)

                mN=no.mN[1]
                mO=no.mO[1]
                mt=mN+mO
                Fimg=F11_image(dz)*mN/mt
            end
        # N-Au interaction
        else
            if neutral_PES_active
                Fn=F00_AuN(distbtwn)
            end
            if ionic_PES_active
                Fi=F11_AuN(distbtwn,cosθ)
                F_AuNc=F11_AuNcutoff(distbtwn,cosθ,rNO,uON)
            end
            if coupled_PES_active
                Fc=F01_AuN(distbtwn)
            end
        end
    elseif i==2
        # ON interaction
        if j==1
            if neutral_PES_active
                Fn=F00_NO(distbtwn)
            end
            if ionic_PES_active
                Fi=F11_NO(distbtwn)

                mN=no.mN[1]
                mO=no.mO[1]
                mt=mN+mO
                Fimg=F11_image(dz)*mO/mt
                if !isempty(FNO_AuN)
                    F_AuNc=F11_AuNcutoff() 
                end
            end
        # O-Au interactions
        else
            if neutral_PES_active
                Fn=F00_AuO(distbtwn)
            end
            if ionic_PES_active
                Fi=F11_AuO(distbtwn)
            end
            if coupled_PES_active
                Fc=F01_AuO(distbtwn)
            end
        end
    else
        # Au-Au interactions
        it=i-2
        jt=j-2

        # find j in nn[i] and return index. needed to access items in variables based off nn
        idx_j=findfirst(isequal(jt), nn[it])

        # no forces on atoms in last layer
        if i<auatomcutoff # need ref to original i since auatomcutoff=399
            # normal operation

            # dist of nn pair away from equilibrium. static array
            rij=vector(r0ij[it][:,idx_j],vec_ij,boundary)

            # force nn pair exerts on atom i. static array
            F_AuAu=SVector{3}(-Aijarray[it][idx_j]*rij)
        else
            # no force on Au atoms
            F_AuAu=SVector{3}(zeros(3)u"N/mol")
        end
    end

    # case for neu+ion?
    if neutral_PES_active && ionic_PES_active && coupled_PES_active
        Fi_mod=Fi*u-[0u"N/mol",0u"N/mol",Fimg]+F_AuNc
        Fg=Fn*a*u+Fi_mod*b+Fc*c*u+F_AuAu
        return Fg
    elseif ionic_PES_active
        Fi_mod=Fi*u-[0u"N/mol",0u"N/mol",Fimg]+F_AuNc
        return Fi_mod+F_AuAu
    else
        return Fn*u+F_AuAu
    end
end

############################################################################################################

@inline @inbounds function getVij_OAu(i,j,distbtwn,dr,dz,boundary)
    En=0u"e_MD"
    Ei=0u"e_MD"
    Ec=0u"e_MD"
    EAu=0u"e_MD"

    if i==1
        if j==1
            # 1 time image, work function, and Ea calc
            Ei=V11_image(dz)+PES_ionic.ϕ[1]-PES_ionic.Ea[1]
        else
            # O-Au interactions
            if neutral_PES_active
                En=V00_AuO(distbtwn)
            end
            if ionic_PES_active
                Ei=V11_AuO(distbtwn)
            end
            if coupled_PES_active
                Ec=V01_AuO(distbtwn)
            end
        end
    else
        # Au-Au interactions
        it=i-1
        jt=j-1

        # find j in nn[i] and return index. needed to access items in variables based off nn
        idx_j=findfirst(isequal(jt), nn[it])

        # dist of nn pair away from equilibrium. static array
        rij=vector(r0ij[it][:,idx_j],dr,boundary)

        # V nn pair exerts on atom i. number. divide by 2 because of v definition
        EAu=dot(rij,Aijarray[it][idx_j],rij)/2
    end

    return En,Ei,Ec,EAu
end

############################################################################################################

@inline @inbounds function getFij_OAu(i,j,vec_ij,dz,boundary)
    Fn=0u"N/mol"
    Fi=0u"N/mol"
    Fc=0u"N/mol"
    Fimg=0u"N/mol"
    F_AuNc=zeros(3)u"N/mol"
    F_AuAu=zeros(3)u"N/mol"

    if i==1
        distbtwn=euclidean(vec_ij,zeros(3)u"Å")
        u=normalize(vec_ij)

        if neutral_PES_active && ionic_PES_active && coupled_PES_active
            # get eigenvalues
            λ1=storeEs[end,2]
            λ2=storeEs[end,3]
            a=λ1^2
            b=λ2^2
            c=2*λ1*λ2

            if j==1
                Fimg=F11_image(dz)
                return -[0u"N/mol",0u"N/mol",Fimg]*b
            else
                Fn=F00_AuO(distbtwn)
                Fi=F11_AuO(distbtwn)
                Fc=F01_AuO(distbtwn)
                Fg=(Fn*a+Fi*b+Fc*c)*u
                return Fg
            end
        elseif ionic_PES_active
            if j==1
                Fimg=F11_image(dz)
                return -[0u"N/mol",0u"N/mol",Fimg]
            else
                Fi=F11_AuO(distbtwn)
                return Fi*u
            end
        else
            Fn=F00_AuO(distbtwn)
            return Fn*u
        end
    else
        # Au-Au interactions
        it=i-1
        jt=j-1

        # find j in nn[i] and return index. needed to access items in variables based off nn
        idx_j=findfirst(isequal(jt), nn[it])

        # no forces on atoms in last layer
        if i<auatomcutoff # need ref to original i since auatomcutoff=398
            # normal operation

            # dist of nn pair away from equilibrium. static array
            rij=vector(r0ij[it][:,idx_j],vec_ij,boundary)

            # force nn pair exerts on atom i. static array
            F_AuAu=SVector{3}(-Aijarray[it][idx_j]*rij)
        else
            # no force on Au atoms
            F_AuAu=SVector{3}(zeros(3)u"N/mol")
        end

        return F_AuAu
    end
end

############################################################################################################

"""
copy of molly's visualize fn (in src/makie.jl) for use w CairoMakie. used if isaac=true
"""
function myvisualize(coord_logger,
                    boundary,
                    out_filepath::AbstractString;
                    connections=Tuple{Int, Int}[],
                    connection_frames=[trues(length(connections)) for i in values(coord_logger)],
                    trails::Integer=0,
                    framerate::Integer=30,
                    color=:purple,
                    connection_color=:orange,
                    # markersize=0.001,
                    markersize=0.0005,
                    linewidth=2.0,
                    transparency=true,
                    show_boundary::Bool=true,
                    boundary_linewidth=2.0,
                    boundary_color=:black,
                    kwargs...)
    coords_start = first(values(coord_logger))
    dist_unit = unit(first(first(coords_start)))
    dims = n_dimensions(boundary)
    fig = Figure()

    if dims == 3
        PointType = Point3f
        # ax = Axis3(fig[1, 1], aspect=:data)
        ax = Axis3(fig[1, 1], aspect=(2.7,1.5,1))
        max_connection_dist = cbrt(box_volume(boundary)) / 2
    elseif dims == 2
        PointType = Point2f
        ax = Axis(fig[1, 1])
        ax.aspect = DataAspect()
        max_connection_dist = sqrt(box_volume(boundary)) / 2
    else
        throw(ArgumentError("Found $dims dimensions but can only visualize 2 or 3 dimensions"))
    end

    positions = Observable(PointType.(ustrip_vec.(coords_start)))
    scatter!(ax, positions; color=color, markersize=markersize, transparency=transparency,
                markerspace=:data, kwargs...)

    if show_boundary
        lines!(
            ax,
            bounding_box_lines(boundary, dist_unit)...;
            color=boundary_color,
            linewidth=boundary_linewidth,
        )
    end

    connection_nodes = []
    for (ci, (i, j)) in enumerate(connections)
        # Don't display connected atoms that are likely connected over the box edge
        if first(connection_frames)[ci] && norm(coords_start[i] - coords_start[j]) < max_connection_dist
            if dims == 3
                push!(connection_nodes, Observable(PointType.(
                        ustrip.([coords_start[i][1], coords_start[j][1]]),
                        ustrip.([coords_start[i][2], coords_start[j][2]]),
                        ustrip.([coords_start[i][3], coords_start[j][3]]))))
            elseif dims == 2
                push!(connection_nodes, Observable(PointType.(
                        ustrip.([coords_start[i][1], coords_start[j][1]]),
                        ustrip.([coords_start[i][2], coords_start[j][2]]))))
            end
        else
            if dims == 3
                push!(connection_nodes, Observable(PointType.([0.0, 0.0], [0.0, 0.0],
                                                        [0.0, 0.0])))
            elseif dims == 2
                push!(connection_nodes, Observable(PointType.([0.0, 0.0], [0.0, 0.0])))
            end
        end
    end
    for (ci, cn) in enumerate(connection_nodes)
        lines!(ax, cn;
                color=isa(connection_color, AbstractArray) ? connection_color[ci] : connection_color,
                linewidth=isa(linewidth, AbstractArray) ? linewidth[ci] : linewidth,
                transparency=transparency)
    end

    trail_positions = []
    for trail_i in 1:trails
        push!(trail_positions, Observable(PointType.(ustrip_vec.(coords_start))))
        col = parse.(Colorant, color)
        alpha = 1 - (trail_i / (trails + 1))
        alpha_col = RGBA.(red.(col), green.(col), blue.(col), alpha)
        scatter!(ax, trail_positions[end]; color=alpha_col,  markersize=markersize,
                    transparency=transparency, markerspace=:data, kwargs...)
    end

    boundary_conv = ustrip.(dist_unit, Molly.cubic_bounding_box(boundary))
    xlims!(ax, Molly.axis_limits(boundary_conv, coord_logger, 1))
    ylims!(ax, Molly.axis_limits(boundary_conv, coord_logger, 2))
    dims == 3 && zlims!(ax, Molly.axis_limits(boundary_conv, coord_logger, 3))

    record(fig, out_filepath, eachindex(values(coord_logger)); framerate=framerate) do frame_i
        coords = values(coord_logger)[frame_i]

        for (ci, (i, j)) in enumerate(connections)
            if connection_frames[frame_i][ci] && norm(coords[i] - coords[j]) < max_connection_dist
                if dims == 3
                    connection_nodes[ci][] = PointType.(
                                ustrip.([coords[i][1], coords[j][1]]),
                                ustrip.([coords[i][2], coords[j][2]]),
                                ustrip.([coords[i][3], coords[j][3]]))
                elseif dims == 2
                    connection_nodes[ci][] = PointType.(
                                ustrip.([coords[i][1], coords[j][1]]),
                                ustrip.([coords[i][2], coords[j][2]]))
                end
            else
                if dims == 3
                    connection_nodes[ci][] = PointType.([0.0, 0.0], [0.0, 0.0],
                                                        [0.0, 0.0])
                elseif dims == 2
                    connection_nodes[ci][] = PointType.([0.0, 0.0], [0.0, 0.0])
                end
            end
        end

        positions[] = PointType.(ustrip_vec.(coords))
        for (trail_i, trail_position) in enumerate(trail_positions)
            trail_position[] = PointType.(ustrip_vec.(values(coord_logger)[max(frame_i - trail_i, 1)]))
        end
    end
end
