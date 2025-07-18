using DrWatson
@quickactivate "sst-mot-analysis"

using GH19, TMI
using ProgressMeter, Downloads

import DrWatson.datadir


"""
    explist -> experiment list
"""
explist() = ("EQ-0015","EQ-1750","OPT-0015")

    file = "https://www.ncei.noaa.gov/pub/data/paleo/gcmoutput/gebbie2019/Theta_anom_OPT-0015.nc"

rooturl() = "https://www.ncei.noaa.gov/pub/data/paleo/gcmoutput/gebbie2019"

rooturl(args...) = joinpath(rooturl(), args...)

function load_GH19_equillibrium(;expt = "OPT-0015", variable = "theta", time_idx = 1)
    A, Alu, γ, TMIfile, L, B = config(GH19.TMIversion()); 
    
    theta = GH19.readfield_snapshot(datadir("Theta_$expt.nc"), variable, time_idx, γ)
    thetap = GH19.readfield_snapshot(datadir("Theta_anom_$expt.nc"), variable, time_idx, γ)
    
    thetabar = theta - thetap

    thetabar_BC = TMI.getsurfaceboundary(thetabar)

    ## preformed phosphate
    thetabar_eq = TMI.steadyinversion(Alu,thetabar_BC,γ)
    
    return thetabar_eq
end


"""
    function download_exp(experiment::String;anomaly=false)

# Arguments
- `experiment::String`: name of experiment, use `explist()` to get possible names
- `anomaly::Bool`: true to load θ anomaly, false to load full θ
# Output
- `outputfile`: name of loaded file, found in the `datadir()` directory
Original Author: Geoffrey Gebbie. 
"""
function download_exp(experiment::String; anomaly=false, force = false)
    
    if anomaly
        filename = "Theta_anom_"*experiment*".nc"
    else
        filename = "Theta_"*experiment*".nc"
    end
    infile = rooturl(filename)
    outfile = datadir(filename)

    if !isfile(outfile) || force == true
        !isdir(datadir()) && mkpath(datadir()) 
        Downloads.download(infile,outfile,verbose=true)
    else
        println("$filename already downloaded; use `force=true` to re-download")
    end
    return outfile
end

function download_GH19_output(;force = false)
    println("Downloading OPT-0015 output\n takes about 10 minutes")
    _ = download_exp("OPT-0015";anomaly=false, force = force)
    _ = download_exp("OPT-0015";anomaly=true, force = force)
    return nothing
end

download_GH19_output()
