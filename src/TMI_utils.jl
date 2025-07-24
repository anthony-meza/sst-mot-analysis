using DrWatson
@quickactivate "sst-mot-analysis"
import DrWatson: datadir, srcdir, plotsdir

using TMI, Interpolations, Statistics, Distributions, LinearAlgebra, SparseArrays
import TMI: config, layerthickness

function config(path, TMIfilename; compute_lu = true)
    TMIfile = path * "/" * TMIfilename

    # println(TMIfile)

    println("Form water-mass matrix A")
    @time A = watermassmatrix(TMIfile)

    # LU factorization for efficient matrix inversions
    if compute_lu
        println("LU factorization of A")
        @time Alu = lu(A)
    else
        println("skip LU factorization of A")
        Alu = nothing
    end

    γ = Grid(TMIfile)
    
    # would be good to make this optional
    println("circulation matrix L=")
    @time L = circulationmatrix(TMIfile,A,γ)

    println("Boundary matrix B=")
    @time B = boundarymatrix(TMIfile,γ)
    
    return  A, Alu, γ, TMIfile, L, B
end

function layerthickness(γ::Grid)
    if γ.depth[1] == 0.0 
        zface= (γ.depth[1:end-1].+γ.depth[2:end])./2;
        dz = ([zface[1] ; diff(zface); 500]);
        return dz
    else
        zcenter = γ.depth; nz = length(zcenter)
        zf = ones(nz + 1); zf[1] = 0.0
        for i in 1:nz
            dzhalf = Float64(zcenter[i] - zf[i])
            zf[i+1] = zf[i] + (Float64(2.0)* dzhalf)
        end
        dz = Float64.((zf[2:end] .- zf[1:end-1]))
        return Float64.(dz)

    end
end