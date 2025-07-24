using DrWatson
@quickactivate "sst-mot-analysis"

using TMI


function LGM_MSL_filter(γ::Grid)
    return (γ.depth .>= 120)
end

function LGM_MSL_grid(γ::Grid)
    MSL_filter = LGM_MSL_filter(γ)
    return Grid((γ.lon,γ.lat,γ.depth[MSL_filter]),
    γ.wet[:, :, MSL_filter],
    γ.interior[:, :, MSL_filter],
    γ.wrap,
    γ.Δ)

end


function cellvolume_lgm(γ::Grid)
    MSL_filter = LGM_MSL_filter(γ)
    return cellvolume(γ).tracer[:, :, MSL_filter]
end


function layerthickness_lgm(γ_MSL_drop::Grid)
    z = γ_MSL_drop.depth
    
    N = length(z)
    f = Vector{typeof(z[1])}(undef, N+1)

    # interior faces
    for i in 1:N-1
        f[i+1] = (z[i] + z[i+1]) / 2
    end

    # extrapolate ends
    f[1]   = z[1] - (f[2] - z[1])
    f[end] = z[end] + (z[end] - f[end-1])

    # cell thicknesses
    dz = diff(f)  # returns a length-N vector
    return dz
end

function readfield_lgm(file,tracername,γ::Grid{A,N}) where {A,N} 
    MSL_filter = LGM_MSL_filter(γ)
    
    # The mode "r" stands for read-only. The mode "r" is the default mode and the parameter can be omitted.
    tracer, units, longname = TMI._read3d(file,tracername)
    T = eltype(tracer)
    tracer = tracer[:, :, MSL_filter]
    
    γ_MSL_drop = LGM_MSL_grid(γ)
    TMI.checkgrid!(tracer,γ_MSL_drop.wet)

    
    c = Field(tracer,γ_MSL_drop,TMI.tracerdict()[tracername],longname,units)
    return c
end

function volumefilled_lgm(TMIversion,Alu,γ)::BoundaryCondition
    v = cellvolume(γ)
    v.tracer[:, :, .!LGM_MSL_filter(γ)] .= 0.0
    area = cellarea(γ)
    vfilled = 0*area # answer will go here
    surf_idx =  findfirst(LGM_MSL_filter(γ))
    surf_wet = γ.wet[:, :, surf_idx]
    # effectively take inverse of transpose A matrix.
    #dVdd = zeros(γ) # pre-allocate array
    dVdd = Alu'\v #[γ.wet]

    # scale the sensitivity value by surface area so that converging meridians are taken into account.
    I = γ.I
    #volume = zeros(Float64,length(γ.lon),length(γ.lat))
    volume = zeros(γ.wet[:,:,1])
    # this step could use a function with γ.I argument

    for ii ∈ eachindex(I)
        if I[ii][3] == 1
            wetarea = area.tracer[I[ii][1],I[ii][2]] #* surf_wet[I[ii][1],I[ii][2]]
            volume[I[ii][1],I[ii][2]] = dVdd.tracer[I[ii][1],I[ii][2],1]./( wetarea) 
        end
    end
             
    volume = log10.(volume)
    volume[.!isfinite.(volume)] .= NaN
    # volume[.!surf_wet] .= NaN

    ∂V∂b  = BoundaryCondition(volume,(γ.lon,γ.lat),γ.depth[1],3,1,γ.wet[:,:,1],:V,"volume filled by surface gridcell","log₁₀(m³/m²)")
    
    return  ∂V∂b 
end