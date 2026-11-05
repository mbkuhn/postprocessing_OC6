using DelimitedFiles
using Interpolations
using Statistics
using Printf
include("FFTAnalysis.jl")
include("get_evenly_spaced_time_series.jl")

lowfrequencies_range = [0.005,0.05]
wavefrequencies_range = [0.05519,0.1345]

# Model-scale center of gravity and full/model length ratio
x_CG = 0.
z_CG = -1.506e-1
rescaling = 50.

# MARIN's constant-filter smoothing
smoothing_type = 2
smoothing_parameter = 30

file_exp = "data/oc6.phase1a.exp0.lc53.txt"
file_exawind_pfx = "data/platform_disp_etc"
file_exawind_sfx = ".dat"
file_exawind_labels = ["0","1","2","3","4","5","6"]
output_dir = "outputs_abridged"

# Data from OC6 experiment (full scale)
EXP_data = Float64.(readdlm(file_exp))
EXP_time = EXP_data[:,1]
EXP_x = EXP_data[:,3]
EXP_z = EXP_data[:,5]
EXP_thy = EXP_data[:,7]

# Join ExaWind restarts, replacing overlapping samples with the newer data.
file_exawind = file_exawind_pfx * file_exawind_labels[1] * file_exawind_sfx
Nalu_data = Float64.(readdlm(file_exawind,skipstart=2))
for label in file_exawind_labels[2:end]
    file_exawind_tmp = file_exawind_pfx * label * file_exawind_sfx
    Nalu_data_tmp = Float64.(readdlm(file_exawind_tmp,skipstart=1))
    last_index = searchsortedfirst(Nalu_data[:,1],Nalu_data_tmp[1,1]) - 1
    global Nalu_data = vcat(Nalu_data[1:last_index,:],Nalu_data_tmp)
end

time_evenly_spaced, Nalu_EvenlySpaced_x =
    get_evenly_spaced_time_series(Nalu_data[:,1],Nalu_data[:,2])
_, Nalu_EvenlySpaced_z =
    get_evenly_spaced_time_series(Nalu_data[:,1],Nalu_data[:,4])

# Convert quaternions to pitch in degrees.
q0 = Nalu_data[:,5]
q1 = Nalu_data[:,6]
q2 = Nalu_data[:,7]
q3 = Nalu_data[:,8]
theta_y_Nalu_rad = -pi/2.0 .+ 2.0.*atan.(sqrt.(1.0 .+ 2.0.*(q0.*q2.-q1.*q3)),
    sqrt.(1.0 .- 2.0.*(q0.*q2.-q1.*q3)))
_, Nalu_EvenlySpaced_thy =
    get_evenly_spaced_time_series(Nalu_data[:,1],rad2deg.(theta_y_Nalu_rad))

NaluWave_time = time_evenly_spaced.*sqrt(rescaling)

# Retained interval in model-scale seconds
t_skip = 37.
t_final = t_skip + 540.
indices = findall(t -> t_skip*sqrt(rescaling) <= t <= t_final*sqrt(rescaling),
    NaluWave_time)
isempty(indices) && error("No ExaWind samples in the configured post-transient interval.")
NaluWave_time = NaluWave_time[indices]
NaluWave_x = (Nalu_EvenlySpaced_x[indices] .- x_CG).*rescaling
NaluWave_z = (Nalu_EvenlySpaced_z[indices] .- z_CG).*rescaling
NaluWave_thy = Nalu_EvenlySpaced_thy[indices]

comparison_start = max(first(EXP_time),first(NaluWave_time))
comparison_end = min(last(EXP_time),last(NaluWave_time))
comparison_start < comparison_end ||
    error("Experiment and ExaWind have no overlapping time interval.")
exp_indices = findall(t -> comparison_start <= t <= comparison_end,EXP_time)
exawind_indices = findall(t -> comparison_start <= t <= comparison_end,NaluWave_time)
length(exp_indices) >= 2 && length(exawind_indices) >= 2 ||
    error("At least two samples from each source are required in the comparison interval.")

exp_time = EXP_time[exp_indices]
exawind_time = NaluWave_time[exawind_indices]
motions = (
    (name="Surge",symbol="x",unit="m",experiment=EXP_x,exawind=NaluWave_x),
    (name="Heave",symbol="z",unit="m",experiment=EXP_z,exawind=NaluWave_z),
    (name="Pitch",symbol="thy",unit="deg",experiment=EXP_thy,exawind=NaluWave_thy),
)

mkpath(output_dir)
@printf("\nTime-series errors: ExaWind - Experiment, full scale [%.3f, %.3f] s\n",
    first(exawind_time),last(exawind_time))
println("Experiment linearly interpolated onto ExaWind times; no phase alignment or mean removal.")
@printf("%-10s %-6s %8s %14s %14s %14s %14s\n",
    "Motion","Unit","Samples","Bias","MAE","RMSE","Max abs. error")

time_plots = Plots.Plot[]
spectra = Plots.Plot[]
for motion in motions
    experiment = motion.experiment[exp_indices]
    exawind = motion.exawind[exawind_indices]
    experiment_itp = interpolate((EXP_time,),motion.experiment,Gridded(Linear()))
    differences = exawind .- experiment_itp.(exawind_time)
    @printf("%-10s %-6s %8d %14.6g %14.6g %14.6g %14.6g\n",
        motion.name,motion.unit,length(differences),mean(differences),
        mean(abs.(differences)),sqrt(mean(abs2,differences)),maximum(abs.(differences)))

    time_plot = plot(exp_time,experiment,label="Experiment",
        ylabel="$(motion.name) [$(motion.unit)]")
    plot!(time_plot,exawind_time,exawind,label="ExaWind",linestyle=:dash)
    push!(time_plots,time_plot)

    f_exp, _, psd_exp, _, _ = FFTAnalysis(exp_time,experiment,
        true,false,smoothing_type,smoothing_parameter)
    f_exawind, _, psd_exawind, _, _ = FFTAnalysis(exawind_time,exawind,
        true,false,smoothing_type,smoothing_parameter)
    spectrum = plot(f_exp,abs.(psd_exp),yscale=:log10,label="Experiment",
        xlims=(0,0.4),xlabel="Frequency [Hz]",
        ylabel="$(motion.name) PSD [$(motion.unit)^2/Hz]")
    plot!(spectrum,f_exawind,abs.(psd_exawind),label="ExaWind",linestyle=:dash)
    vspan!(spectrum,wavefrequencies_range,alpha=0.3,label="Wave-frequency range")
    vspan!(spectrum,lowfrequencies_range,alpha=0.3,label="Low-frequency range")
    savefig(spectrum,joinpath(output_dir,"$(motion.symbol)_PSD.pdf"))
    push!(spectra,spectrum)
end

displacement_time_plot = plot(time_plots...,layout=(3,1),link=:x,size=(1000,900),
    xlabel="Time [s] - Full Scale",legend=:outertopright)
savefig(displacement_time_plot,joinpath(output_dir,"displacement_time_series.pdf"))
