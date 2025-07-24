using DrWatson
@quickactivate "sst-mot-analysis"


"""
    function download_exp(experiment::String;anomaly=false)

# Arguments
- `experiment::String`: name of experiment, use `explist()` to get possible names
- `anomaly::Bool`: true to load θ anomaly, false to load full θ
# Output
- `outputfile`: name of loaded file, found in the `datadir()` directory
"""
function download_exp(experiment::String,anomaly=false;force = false)
    
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
        println("File already downloaded; use `force=true` to re-download")
    end
    return outfile
end

function download_GH19_output()
    download_exp("OPT-0015",anomaly=false;force = false)
    download_exp("OPT-0015",anomaly=true;force = false)
end