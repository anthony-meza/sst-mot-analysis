using DrWatson
@quickactivate "sst-mot-analysis"

using GH19, TMI

import DrWatson.datadir

function load_GH19_equillibrium(;expt = "OPT-0015", variable = "theta", time_idx = 1)
    A, Alu, γ, TMIfile, L, B = config(GH19.TMIversion()); 
    
    theta = GH19.readfield_snapshot(datadir("Theta_$expt.nc"), variable, time_idx, γ)
    thetap = GH19.readfield_snapshot(datadir("Theta_anom_$expt.nc"), variable, time_idx, γ)
    
    thetabar = theta - thetap
    return thetabar
end