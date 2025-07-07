using DrWatson
@quickactivate "sst-mot-analysis"
import DrWatson: plotsdir

using TMI, GH19, PythonPlot, ProgressMeter
using Printf, Statistics
using Colors, LaTeXStrings
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
                                              sampling_method,
                                              N_sample, Nboot,
                                              output_filename)
    # unpack inputs
    @unpack bootstrapped_profiles, LGM_theta, PI_theta, γ_PI, γ_lgm = bootstrap_results

    # compute deltas
    delta_sst = vec(mean(bootstrapped_profiles["PI_surface"] .-
                         bootstrapped_profiles["LGM_surface"], dims=1))
    delta_mot = vec(mean(bootstrapped_profiles["PI_bottom"] .-
                         bootstrapped_profiles["LGM_bottom"], dims=1))

    # sampling description
    sampling_description = sampling_method == :uniform ?
        "uniform grid sampling" : "area-weighted spherical sampling"

    # reference values
    SeltΔMOT = 2.27; SeltΔMOT_σ = 0.46

    # figure + axes
    fig, ax = subplots(figsize=(7,7), dpi=250)

    # x-range for reference lines
    x_min, x_max = 0, 5
    xs = collect(range(x_min, x_max, length=200))

    # Seltzer line + error ribbon
    ys = fill(SeltΔMOT, length(xs))
    ax.plot(xs, ys;
            color="hotpink", linewidth=1.5,
            label="Seltzer et al., 2024", zorder = 1)
    ax.fill_between(xs,
                    ys .- SeltΔMOT_σ, ys .+ SeltΔMOT_σ;
                    color="hotpink", alpha=0.2, zorder = 0)

    # bootstrap means scatter
    ax.scatter(delta_sst, delta_mot;
               s=9, alpha=0.4,
               color="dodgerblue",
               linewidths=0.2,
                label = "Pseudo Sample Average\n" * L"($n_s = $" * "$N_sample" * L"$)$",
               # label="Bootstrap avgs: $Nboot dots\n(mean of $N_sample samples, $sampling_description)", 
                zorder = 8)

    # 1:1 reference
    ax.plot(xs, xs, linestyle="--", linewidth=1.0, color="k")

    # bootstrapped mean point
    ax.scatter([mean(delta_sst)], [mean(delta_mot)],
               marker="s",
               s=100, edgecolor="k",
               facecolor="dodgerblue",
               linewidths=1.5,
               label = "All Pseudo Sample Averages\n" * L"($n_{t} = n_s $" * L"\times" * "$Nboot" * L"$)$",
                zorder = 10)

    # true global-average point
    global_sst_diff = global_surface_average(PI_theta, γ_PI) -
                      global_surface_average(LGM_theta, γ_lgm)
    global_mot_diff = global_ocean_average(PI_theta, γ_PI) -
                      global_ocean_average(LGM_theta, γ_lgm)
    
    ax.scatter([global_sst_diff], [global_mot_diff];
               marker="s",
               s=100, edgecolor="k", 
               facecolor="sienna",
               linewidths=1.5,
               label="TMI Volume Weighted Differences\n(GH19 PI – G14 LGM)", zorder = 10)

    # labels, title, limits, legend
    # ax.set_xlabel(raw"$\overline{\Delta \mathrm{SST}}\;(^\circ\mathrm{C})$")
    # ax.set_ylabel(raw"$\overline{\Delta \mathrm{MOT}}\;(^\circ\mathrm{C})$")
    ax.set_ylabel("Mean Ocean Temperature Change" * raw"$(^\circ\mathrm{C})$")
    ax.set_xlabel("Mean SST Change" * raw"$(^\circ\mathrm{C})$")
    # ax.set_title(raw"$\Delta\,$Ocean Temp. (LGM vs. PI)")
    ax.set_xlim(x_min, x_max)
    ax.set_ylim(x_min, x_max)
    ax.legend(loc="upper right", fontsize=8, ncols = 2)
    ax.spines["top"].set_visible(false)
    ax.spines["right"].set_visible(false)
    # layout + save
    # fig.tight_layout()
    fig.savefig(plotsdir(output_filename), dpi = 200, bbox_inches = "tight")
end
