using DrWatson
@quickactivate "sst-mot-analysis"
import DrWatson: plotsdir

using TMI, GH19, Plots, ProgressMeter
using Printf, Statistics
using Colors: colorant, hex
include("random_profiles.jl")
include("load_GH19.jl")

last_non_nan(v) = findlast(!isnan, v) !== nothing ? v[findlast(!isnan, v)] : NaN
last_non_nan(m::Matrix) = last_non_nan.(eachcol(m))
global_ocean_average(x::AbstractArray, γ::Grid) = sum(replace(x, NaN=>0.0) .* cellvolume(γ).tracer) / sum(cellvolume(γ).tracer)
global_surface_average(x::AbstractArray, γ::Grid) = sum(replace(x[:, :, 1], NaN=>0.0) .* cellvolume(γ).tracer[:, :, 1]) / sum(cellvolume(γ).tracer[:, :, 1])

function bootstrap_PI_lgm_differences(N_sample, Nboot; sampling_method=:uniform)
    # Configure LGM dataset
    TMIversion_lgm = "LGM_90x45x33_G14"
    A_lgm, Alu_lgm, γ_lgm, TMIfile_lgm, L_lgm, B_lgm = config(TMIversion_lgm)
    LGM_theta = readfield(TMIfile_lgm, "θ", γ_lgm).tracer
    
    # Configure PI dataset
    TMIversion_PI = "modern_180x90x33_GH11_GH12"
    A_PI, Alu_PI, γ_PI, TMIfile_PI, L_PI, B_PI = config(TMIversion_PI)
    PI_theta = load_GH19_equillibrium().tracer
    
    # Initialize bootstrap profiles dictionary
    bootstrapped_profiles = Dict(
        "PI_surface" => zeros(N_sample, Nboot), "PI_bottom" => zeros(N_sample, Nboot), 
        "LGM_surface" => zeros(N_sample, Nboot), "LGM_bottom" => zeros(N_sample, Nboot)
    )
    
    # iter = ProgressBar()
    desc = "Running $Nboot MC trials (n=$N_sample).."
    @showprogress dt=0.1 desc=desc for ib in 1:Nboot
        # Generate sampling locations
        locs = Vector{Tuple{Float64, Float64}}(undef, N_sample)
        [locs[i] = wetsurfacelocation(γ_lgm, γ_PI; sampling_method=sampling_method) for i in eachindex(locs)]

        # Sample LGM and PI temperature profiles
        y_lgm, _, _, _ = random_profiles(TMIversion_lgm, "θ", γ_lgm, N_sample; locs=locs)
        bootstrapped_profiles["LGM_surface"][:, ib] .= y_lgm[1, :] 
        bootstrapped_profiles["LGM_bottom"][:, ib] .= last_non_nan(y_lgm)

        y_PI, _, _, _ = random_profiles(PI_theta, γ_PI, N_sample; locs=locs)
        bootstrapped_profiles["PI_surface"][:, ib] .= y_PI[1, :] 
        bootstrapped_profiles["PI_bottom"][:, ib] .= last_non_nan(y_PI)
        # set_description(iter, string(@sprintf("Bootstrap Iter: %.2f", ib)))

    end
    results_dict = @strdict bootstrapped_profiles LGM_theta PI_theta γ_PI γ_lgm
    return results_dict
end

function generate_temperature_difference_plot(bootstrap_results,  
                                             sampling_method, N_sample, Nboot, output_filename)
    # Calculate average differences
    @unpack bootstrapped_profiles, LGM_theta, PI_theta, γ_PI, γ_lgm = bootstrap_results
        
    delta_sst = mean(bootstrapped_profiles["PI_surface"] .- bootstrapped_profiles["LGM_surface"], dims=1)[:]
    delta_mot = mean(bootstrapped_profiles["PI_bottom"] .- bootstrapped_profiles["LGM_bottom"], dims=1)[:]

    # Define consistent plot style for all plots
    plot_style = Dict(
        :xlabel => "\$\\overline{\\Delta{SST}}(^\\circ C)\$", 
        :ylabel => "\$\\overline{\\Delta{MOT}}(^\\circ C)\$",
        :title => "\$\\Delta\$Ocean Temp. (LGM vs. PI)", 
        :markerstrokewidth => 0.1, :markersize => 3,
        :legend => :bottomright, 
        :right_margin => 2Plots.mm, :left_margin => 2Plots.mm,
        :bottom_margin => 2Plots.mm,  # Add extra margin for the legend
        :size => (700, 700), :dpi => 1000,
        :titlefontsize=>18, :guidefontsize=>15,
        :tickfontsize=>13, :legendfontsize=>10,        
        :xlims => (-2, 5), :ylims => (-2, 5),
    )
    
    sampling_description = sampling_method == :uniform ? "uniform grid sampling" : "area-weighted spherical sampling"
    
    # Create the error bars from Seltzer 2024 (LGM - PI diff) 
    SeltΔMOT = 2.27; SeltΔMOT_1σerr = 0.46

    p = plot(collect(plot_style[:xlims]), fill(SeltΔMOT, length(collect(plot_style[:xlims]))),
            grid=false,ribbon=SeltΔMOT_1σerr,
            linewidth = 1.5,
            color = :deeppink,fillalpha=.2, 
            label = "Seltzer et. al., 2024")
    
    scatter!(delta_sst, delta_mot, alpha = 0.5,
        label="Bootstrap avgs, $Nboot dots total, each dot is avg of $N_sample samples\n($sampling_description)", color=:dodgerblue; plot_style...)
    plot!(collect(plot_style[:xlims]), collect(plot_style[:xlims]), label="1:1",
          linestyle=:dash, linewidth=2, color=:black)
    
    # Add mean and global average points
    highlight = Dict(:markerstrokewidth => 1.5, :markersize => 10, :order => 10, :alpha => 1.0)
    
    scatter!([mean(delta_sst)], [mean(delta_mot)], 
        label="Mean diff averaged over all samples", color=:dodgerblue, 
        marker = :rect; 
        highlight...)
    
    global_sst_diff = global_surface_average(PI_theta, γ_PI) - global_surface_average(LGM_theta, γ_lgm)
    global_mot_diff = global_ocean_average(PI_theta, γ_PI) - global_ocean_average(LGM_theta, γ_lgm)
    
    scatter!([global_sst_diff], [global_mot_diff], 
        label="TMI Volume Weighted Differences\n(GH19 PI - G14 LGM)", 
        color=:sienna4,   marker= :rect,; highlight...)
    
    savefig(plotsdir(output_filename))
    return p
end