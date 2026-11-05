using DelimitedFiles
using Interpolations
using Statistics
using Printf
include("FFTAnalysis.jl")
include("get_evenly_spaced_time_series.jl")

# Model-scale center of gravity and full/model length ratio
x_CG = 0.
z_CG = -1.506e-1
rescaling = 50.

# Convert full-scale frequency bands to model scale.
lowfrequencies_range = [0.005,0.05].*sqrt(rescaling)
wavefrequencies_range = [0.05519,0.1345].*sqrt(rescaling)
frequency_limits = (0.,0.4*sqrt(rescaling))

# MARIN's constant-filter smoothing
smoothing_type = 2
smoothing_parameter = 30

file_new_run_placeholder = "data/oc6.phase1a.exp0.lc53.txt"
file_current_run_pfx = "data/platform_disp_etc"
file_current_run_sfx = ".dat"
file_current_run_labels = ["0","1","2","3","4","5","6"]
output_dir = "outputs_abridged"
current_run_label = "ExaWind - current run"
new_run_label = "ExaWind - new run (experiment placeholder)"

# Placeholder for the new run: convert experimental full-scale data to model scale.
# A real ExaWind replacement must supply model-scale displacements and pitch in degrees.
new_run_data = Float64.(readdlm(file_new_run_placeholder))
new_run_time = new_run_data[:,1]./sqrt(rescaling)
new_run_x = new_run_data[:,3]./rescaling
new_run_z = new_run_data[:,5]./rescaling
new_run_thy = new_run_data[:,7]

# Join current-run restarts, replacing overlapping samples with the newer data.
file_current_run = file_current_run_pfx * file_current_run_labels[1] * file_current_run_sfx
current_run_data = Float64.(readdlm(file_current_run,skipstart=2))
for label in file_current_run_labels[2:end]
    file_restart = file_current_run_pfx * label * file_current_run_sfx
    restart_data = Float64.(readdlm(file_restart,skipstart=1))
    last_index = searchsortedfirst(current_run_data[:,1],restart_data[1,1]) - 1
    global current_run_data = vcat(current_run_data[1:last_index,:],restart_data)
end

current_run_time, current_run_x =
    get_evenly_spaced_time_series(current_run_data[:,1],current_run_data[:,2] .- x_CG)
_, current_run_z =
    get_evenly_spaced_time_series(current_run_data[:,1],current_run_data[:,4] .- z_CG)

# Convert quaternions to pitch in degrees (unchanged by geometric scaling).
q0 = current_run_data[:,5]
q1 = current_run_data[:,6]
q2 = current_run_data[:,7]
q3 = current_run_data[:,8]
pitch_rad = -pi/2.0 .+ 2.0.*atan.(sqrt.(1.0 .+ 2.0.*(q0.*q2.-q1.*q3)),
    sqrt.(1.0 .- 2.0.*(q0.*q2.-q1.*q3)))
_, current_run_thy =
    get_evenly_spaced_time_series(current_run_data[:,1],rad2deg.(pitch_rad))

# Model-scale seconds: retain all shared post-transient data unless an end is specified.
t_skip = 37.
t_final = Inf
comparison_start = max(t_skip,first(new_run_time),first(current_run_time))
comparison_end = min(t_final,last(new_run_time),last(current_run_time))
comparison_start < comparison_end ||
    error("Current and new runs have no overlapping time interval within the configured limits.")
new_indices = findall(t -> comparison_start <= t <= comparison_end,new_run_time)
current_indices = findall(t -> comparison_start <= t <= comparison_end,current_run_time)
length(new_indices) >= 2 && length(current_indices) >= 2 ||
    error("At least two samples from each run are required in the retained interval.")

new_time = new_run_time[new_indices]
current_time = current_run_time[current_indices]
motions = (
    (name="Surge",symbol="x",unit="m",new_run=new_run_x,current_run=current_run_x),
    (name="Heave",symbol="z",unit="m",new_run=new_run_z,current_run=current_run_z),
    (name="Pitch",symbol="thy",unit="deg",new_run=new_run_thy,current_run=current_run_thy),
)

mkpath(output_dir)
@printf("\nRetained model-scale interval: [%.3f, %.3f] s\n",comparison_start,comparison_end)
println("New run: experimental placeholder converted to model scale.")
@printf("Time-series errors: new - current, model scale [%.3f, %.3f] s\n",
    first(new_time),last(new_time))
println("Current run linearly interpolated onto new-run times; no phase alignment or mean removal.")
@printf("%-10s %-6s %8s %14s %14s %14s %14s\n",
    "Motion","Unit","Samples","Bias","MAE","RMSE","Max abs. error")

time_plots = Plots.Plot[]
spectra = Plots.Plot[]
for motion in motions
    new_values = motion.new_run[new_indices]
    current_values = motion.current_run[current_indices]
    current_itp = interpolate((current_run_time,),motion.current_run,Gridded(Linear()))
    differences = new_values .- current_itp.(new_time)
    @printf("%-10s %-6s %8d %14.6g %14.6g %14.6g %14.6g\n",
        motion.name,motion.unit,length(differences),mean(differences),
        mean(abs.(differences)),sqrt(mean(abs2,differences)),maximum(abs.(differences)))

    time_plot = plot(current_time,current_values,label=current_run_label,
        ylabel="$(motion.name) [$(motion.unit)]")
    plot!(time_plot,new_time,new_values,label=new_run_label,linestyle=:dash)
    push!(time_plots,time_plot)

    f_current, _, psd_current, _, _ = FFTAnalysis(current_time,current_values,
        true,false,smoothing_type,smoothing_parameter)
    f_new, _, psd_new, _, _ = FFTAnalysis(new_time,new_values,
        true,false,smoothing_type,smoothing_parameter)
    spectrum = plot(f_current,abs.(psd_current),yscale=:log10,label=current_run_label,
        xlims=frequency_limits,xlabel="Frequency [Hz] - Model Scale",
        ylabel="$(motion.name) PSD [$(motion.unit)^2/Hz] - Model Scale")
    plot!(spectrum,f_new,abs.(psd_new),label=new_run_label,linestyle=:dash)
    vspan!(spectrum,wavefrequencies_range,alpha=0.3,label="Wave-frequency range")
    vspan!(spectrum,lowfrequencies_range,alpha=0.3,label="Low-frequency range")
    savefig(spectrum,joinpath(output_dir,"$(motion.symbol)_PSD.pdf"))
    push!(spectra,spectrum)
end

displacement_time_plot = plot(time_plots...,layout=(3,1),link=:x,size=(1200,900),
    xlabel="Time [s] - Model Scale",legend=:outertopright)
savefig(displacement_time_plot,joinpath(output_dir,"displacement_time_series.pdf"))
