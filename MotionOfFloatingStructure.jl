using DelimitedFiles
using Interpolations
using DataFrames
using LaTeXStrings
using StatsBase
using MAT
using Measurements
using Printf
include("FFTAnalysis.jl")
include("rangedTrapz.jl")
include("get_evenly_spaced_time_series.jl")

# Set relevant frequency ranges based on MARIN's set-up
lowfrequencies_range = [0.005 0.05]
wavefrequencies_range = [0.05519 0.1345]

# Set relevant Milestone's +- validation ranges
range_wave_frequency = 0.1
range_low_frequency = 0.2
plus_range_wave_frequency = 1. + range_wave_frequency
minus_range_wave_frequency = 1. - range_wave_frequency
plus_range_low_frequency = 1. + range_low_frequency
minus_range_low_frequency = 1. - range_low_frequency
wave_range_str = @sprintf("%2.0f",100*range_wave_frequency)
low_range_str = @sprintf("%2.0f",100*range_low_frequency)

# Model scale CG: to calculate displacements
x_CG = 0.
y_CG = 0.
z_CG = -1.506e-1

# Define smoothing type and relevant parameters
# Options are: 
# Ensemble Averaging -> smoothing_type = 1; smoothing_parameter = 10000
# Constant Filter (MARIN's Method) -> smoothing_type = 2; smoothing_parameter = 30
# Gaussian -> smoothing_type = 3; smoothing_parameter = 30
smoothing_type = 2
smoothing_parameter = -1
if smoothing_type==1
    smoothing_parameter = 10000
else
    smoothing_parameter = 30
end

file_exp = "/Users/mkuhn/testruns_data/OC6_PhaseIa_dataset_and_more/LC5_motion_moored_waves/oc6.phase1a.exp0.lc53.txt"
file_cfd_ref = "/Users/mkuhn/testruns_data/OC6_PhaseIa_dataset_and_more/LC5_motion_moored_waves/oc6.phase1a.nrel.lc53.txt"
file_exawind_pfx = "/Users/mkuhn/testruns_data/OC5_semisubmersible/HOS-NWT/floating_fixed_dt/platform_disp_etc"
file_exawind_sfx = ".dat"
file_exawind_labels = ["0","1","2","3"]


####################### FORCES ##########################################################################################################
# Data from OC6 Experiment
ZXP0_data = Float64.(readdlm(file_exp))
EXP_time = ZXP0_data[:,1]
EXP_wave = ZXP0_data[:,2]
EXP_x = ZXP0_data[:,3] 
EXP_z = ZXP0_data[:,5] 
EXP_thy = ZXP0_data[:,7]

# Data from StarCCM+ simulation
CFD_data = Float64.(readdlm(file_cfd_ref))
CFD_time = CFD_data[:,1]
CFD_wave = CFD_data[:,2]
CFD_x = CFD_data[:,3] 
CFD_z = CFD_data[:,5] 
CFD_thy = CFD_data[:,7]

# NaluWind DATA

# Conversion between full scale and model scale
# L_p = 50 L_m
# T_p = sqrt(50)T_m
# F_p = L^3 F_m
# M_y = L^4 M_y
# theta_p = theta_m
rescaling = 50. 

file_exawind = file_exawind_pfx * file_exawind_labels[1] * file_exawind_sfx
Nalu_data = readdlm(file_exawind,skipstart=2)
for label in file_exawind_labels[2:end]
    file_exawind_tmp = file_exawind_pfx * label * file_exawind_sfx
    Nalu_data_tmp = readdlm(file_exawind_tmp,skipstart=1)
    first_time = Nalu_data_tmp[1,1]
    last_index = length(Nalu_data[:,1])
    for ii=1:length(Nalu_data[:,1])
        last_time = Nalu_data[last_index,1]
        if (last_time >= first_time)
            last_index -= 1
        else
            break
        end
    end
    global Nalu_data = vcat(Nalu_data[1:last_index,:],Nalu_data_tmp)
end

NaluWAVE_data = Float64.(Nalu_data)
# Make unevenly space time series into an evenly spaced one
dt_array = zeros(length(NaluWAVE_data[:,1]))
for ii=1:length(dt_array)-1
    dt = NaluWAVE_data[ii+1,1].-NaluWAVE_data[ii,1]
    dt_array[ii] = dt
end
t_spacing = median(dt_array)
# Create evenly-spaced time array
time_evenly_spaced = collect(NaluWAVE_data[1,1]:t_spacing:NaluWAVE_data[end,1])

# Interpolate the unevely-spaced time series
# x
_ , Nalu_EvenlySpaced_x = get_evenly_spaced_time_series(NaluWAVE_data[:,1],NaluWAVE_data[:,2])
# z
_ , Nalu_EvenlySpaced_z = get_evenly_spaced_time_series(NaluWAVE_data[:,1],NaluWAVE_data[:,4])
# convert from quaternions to angle
q0 = NaluWAVE_data[:,5]
q1 = NaluWAVE_data[:,6]
q2 = NaluWAVE_data[:,7]
q3 = NaluWAVE_data[:,8]
theta_y_Nalu_rad = -pi/2.0.+ 2.0.*atan.(sqrt.(1.0.+2.0.*(q0.*q2.-q1.*q3)),sqrt.(1.0.-2.0.*(q0.*q2.-q1.*q3)))
theta_y_Nalu = theta_y_Nalu_rad.*(180.0./pi)
# convert to degrees
# theta y
_ , Nalu_EvenlySpaced_thy = get_evenly_spaced_time_series(NaluWAVE_data[:,1],theta_y_Nalu)

# Convert time to full scale
NaluWave_time = time_evenly_spaced.*sqrt(rescaling)

# In case one needs to discard any transient time data
#t_skip = NaluWave_time[1]./sqrt(rescaling) + 0.027*(NaluWave_time[end]./sqrt(rescaling)) # seconds (model scale) # if skipping some times is needed 
#t_final = t_skip + 545.0 # seconds (model scale)
t_skip = 39.1
t_final = t_skip + 545.
println("t_skip = ",t_skip)
indices = findall(NaluWave_time->(NaluWave_time>=(t_skip*sqrt(rescaling)) && NaluWave_time<=(t_final*sqrt(rescaling))),NaluWave_time)

NaluWave_time = NaluWave_time[indices]
println("Initial time (Model Scale) = ", NaluWave_time[1]/sqrt(rescaling))
println("Final time (Model Scale) = ", NaluWave_time[end]/sqrt(rescaling))

# Convert x, z, and theta y to full scale after skipping transient portion of array
# Also subtract center of gravity to get displacements, not position
NaluWave_x = (Nalu_EvenlySpaced_x[indices] .- x_CG).*(rescaling)
NaluWave_z = (Nalu_EvenlySpaced_z[indices] .- z_CG).*(rescaling)
NaluWave_thy = (Nalu_EvenlySpaced_thy[indices])

stepRange = collect(3000:108000) # We are skipping 5 minutes (300 seconds) of full scale time (roughly 42 seconds of model scale time), out of a 10800 seconds which are 180 minutes (3 hours) of full scale time (26 minutes of model scale)

stepRangeNalu = collect(1:length(NaluWave_time)) # include all of the Nalu time

# Mean Surge - Unrelated to PSD, just the mean value of the surge force across times
mean_Nalu_x = mean(NaluWave_x[stepRangeNalu])
mean_EXP_x = mean(EXP_x[stepRange])
mean_CFD_x = mean(CFD_x[stepRange])
println("Nalu Mean Surge = ", mean_Nalu_x)
println("Exp Mean Surge = ", mean_EXP_x)
println("CFD Mean Surge = ", mean_CFD_x)
println("Ratio = ", mean_Nalu_x/mean_EXP_x)
mean_surge_plot = bar((1,mean_EXP_x),label="Experiment",legend=:outertopright)
bar!((2,mean_CFD_x),label="StarCCM+")
bar!((3,mean_Nalu_x),label="ExaWind")
hline!(ones(length(mean_EXP_x))*(abs.(mean_EXP_x)*plus_range_wave_frequency),color=:red,linestyle=:dash,label="+"*wave_range_str*"% from Experiment")
hline!(ones(length(mean_EXP_x))*(abs.(mean_EXP_x)*minus_range_wave_frequency),color=:black,linestyle=:dash,label="-"*wave_range_str*"% from Experiment")
ylabel!("Mean Drift in Surge \$[m]\$")
savefig(mean_surge_plot,"Mean_Surge_Drift.pdf")

### **** PSD COMPUTATIONS AND PLOTS ******

# first one is for plotting purposes
# x EXP
f, P1, PSD, fdouble, P2 = FFTAnalysis(EXP_time[stepRange],EXP_x[stepRange],true,false,smoothing_type,smoothing_parameter)
x_plot = plot(f[1:end],abs.(PSD[1:end]),yscale=:log10,label="Experiment",xlims=(0,0.4)) #,ylims=(1e10,1e16))
f, P1, PSD, fdouble, P2 = FFTAnalysis(EXP_time[stepRange],EXP_x[stepRange],false,false,~,~)
PSDInt_x_EXP_low  = rangedTrapz(f,abs.(PSD),lowfrequencies_range)
PSDInt_x_EXP_wave = rangedTrapz(f,abs.(PSD),wavefrequencies_range)

# x CFD
f, P1, PSD, fdouble, P2 = FFTAnalysis(CFD_time[stepRange],CFD_x[stepRange],true,false,smoothing_type,smoothing_parameter)
plot!(f[1:end],abs.(PSD[1:end]),yscale=:log10,label="StarCCM+")
f, P1, PSD, fdouble, P2 = FFTAnalysis(CFD_time[stepRange],CFD_x[stepRange],false,false,~,~)
PSDInt_x_CFD_low  = rangedTrapz(f,abs.(PSD),lowfrequencies_range)
PSDInt_x_CFD_wave = rangedTrapz(f,abs.(PSD),wavefrequencies_range)

# NALU x
f, P1, PSD, fdouble, P2 = FFTAnalysis(NaluWave_time[stepRangeNalu],NaluWave_x[stepRangeNalu],true,false,smoothing_type,smoothing_parameter)
plot!(f[1:end],abs.(PSD[1:end]),yscale=:log10,label="ExaWind")
vspan!([0.05519,0.1345]; alpha = 0.3,label="Frequency Range")
vspan!([0.005,0.05]; alpha = 0.3,label="Low Frequency Range")
xlabel!("Frequency [Hz]")
ylabel!("Global \$x [m^2/Hz]\$")
savefig(x_plot,"x_PSD.pdf")
f, P1, PSD, fdouble, P2 = FFTAnalysis(NaluWave_time[stepRangeNalu],NaluWave_x[stepRangeNalu],false,false,~,~)
PSDInt_x_NaluWave_low  = rangedTrapz(f,abs.(PSD),lowfrequencies_range)
PSDInt_x_NaluWave_wave = rangedTrapz(f,abs.(PSD),wavefrequencies_range) # these are just to get the bar plots

# first one is for plotting purposes
# z EXP
f, P1, PSD, fdouble, P2 = FFTAnalysis(EXP_time[stepRange],EXP_z[stepRange],true,false,smoothing_type,smoothing_parameter)
z_plot = plot(f[1:end],abs.(PSD[1:end]),yscale=:log10,label="Experiment",xlims=(0,0.4))
f, P1, PSD, fdouble, P2 = FFTAnalysis(EXP_time[stepRange],EXP_z[stepRange],false,false,~,~)
PSDInt_z_EXP_low  = rangedTrapz(f,abs.(PSD),lowfrequencies_range)
PSDInt_z_EXP_wave = rangedTrapz(f,abs.(PSD),wavefrequencies_range)

# z CFD
f, P1, PSD, fdouble, P2 = FFTAnalysis(CFD_time[stepRange],CFD_z[stepRange],true,false,smoothing_type,smoothing_parameter)
plot!(f[1:end],abs.(PSD[1:end]),yscale=:log10,label="StarCCM+")
f, P1, PSD, fdouble, P2 = FFTAnalysis(CFD_time[stepRange],CFD_z[stepRange],false,false,~,~)
PSDInt_z_CFD_low  = rangedTrapz(f,abs.(PSD),lowfrequencies_range)
PSDInt_z_CFD_wave = rangedTrapz(f,abs.(PSD),wavefrequencies_range)

# NALU z
f, P1, PSD, fdouble, P2 = FFTAnalysis(NaluWave_time[stepRangeNalu],NaluWave_z[stepRangeNalu],true,false,smoothing_type,smoothing_parameter)
plot!(f[1:end],abs.(PSD[1:end]),yscale=:log10,label="ExaWind") # need to investigate about this truncation
vspan!([0.05519,0.1345]; alpha = 0.3,label="Frequency Range")
vspan!([0.005,0.05]; alpha = 0.3,label="Low Frequency Range")
xlabel!("Frequency [Hz]")
ylabel!("Global \$z [m^2/Hz]\$")
savefig(z_plot,"z_PSD.pdf")
f, P1, PSD, fdouble, P2 = FFTAnalysis(NaluWave_time[stepRangeNalu],NaluWave_z[stepRangeNalu],false,false,~,~)
PSDInt_z_NaluWave_low  = rangedTrapz(f,abs.(PSD),lowfrequencies_range)
PSDInt_z_NaluWave_wave = rangedTrapz(f,abs.(PSD),wavefrequencies_range) # these are just to get the bar plots

# first one is for plotting purposes
# theta y EXP
f, P1, PSD, fdouble, P2 = FFTAnalysis(EXP_time[stepRange],EXP_thy[stepRange],true,false,smoothing_type,smoothing_parameter)
thy_plot = plot(f[1:end],abs.(PSD[1:end]),yscale=:log10,label="Experiment",xlims=(0,0.4)) #,ylims=(1e13,1e18),legend=:topright) # in the paper, the xlim is [0,0.4], in the code is [0,0.14]
f, P1, PSD, fdouble, P2 = FFTAnalysis(EXP_time[stepRange],EXP_thy[stepRange],false,false,~,~)
PSDInt_thy_EXP_low  = rangedTrapz(f,abs.(PSD),lowfrequencies_range)
PSDInt_thy_EXP_wave = rangedTrapz(f,abs.(PSD),wavefrequencies_range)

# theta y CFD
f, P1, PSD, fdouble, P2 = FFTAnalysis(CFD_time[stepRange],CFD_thy[stepRange],true,false,smoothing_type,smoothing_parameter)
plot!(f[1:end],abs.(PSD[1:end]),yscale=:log10,label="StarCCM+")
f, P1, PSD, fdouble, P2 = FFTAnalysis(CFD_time[stepRange],CFD_thy[stepRange],false,false,~,~)
PSDInt_thy_CFD_low  = rangedTrapz(f,abs.(PSD),lowfrequencies_range)
PSDInt_thy_CFD_wave = rangedTrapz(f,abs.(PSD),wavefrequencies_range)

# theta y Nalu
f, P1, PSD, fdouble, P2 = FFTAnalysis(NaluWave_time[stepRangeNalu],NaluWave_thy[stepRangeNalu],true,false,smoothing_type,smoothing_parameter)
plot!(f[1:end],abs.(PSD[1:end]),yscale=:log10,label="ExaWind") # need to investigate about this truncation
vspan!([0.05519,0.1345]; alpha = 0.3,label="Frequency Range")
vspan!([0.005,0.05]; alpha = 0.3,label="Low Frequency Range")
xlabel!("Frequency [Hz]")
ylabel!("Global " * L"\theta_y [^o ^2]")
savefig(thy_plot,"thy_PSD.pdf")
f, P1, PSD, fdouble, P2 = FFTAnalysis(NaluWave_time[stepRangeNalu],NaluWave_thy[stepRangeNalu],false,false,~,~)
PSDInt_thy_NaluWave_low  = rangedTrapz(f,abs.(PSD),lowfrequencies_range)
PSDInt_thy_NaluWave_wave = rangedTrapz(f,abs.(PSD),wavefrequencies_range) # these are just to get the bar plots

####### PSD INTEGRAL PLOTS (Bar plots) ##################################################
bar_x = bar((1,abs.(PSDInt_x_EXP_low)),label="Experiment",legend=:outertopright)
bar!((2,abs.(PSDInt_x_CFD_low)),label="StarCCM+")
bar!((3,abs.(PSDInt_x_NaluWave_low)),label="ExaWind")
hline!(ones(length(PSDInt_x_EXP_low))*(abs.(PSDInt_x_EXP_low)*plus_range_low_frequency),color=:red,linestyle=:dash,label="+"*low_range_str*"% from Experiment")
hline!(ones(length(PSDInt_x_EXP_low))*(abs.(PSDInt_x_EXP_low)*minus_range_low_frequency),color=:black,linestyle=:dash,label="-"*low_range_str*"% from Experiment")
ylabel!("Low-Frequency PSD Integral of \$x [m^2]\$")
savefig(bar_x,"bar_x_PSD_low_frequency.pdf")

println("x_low")
println(PSDInt_x_EXP_low)
println(PSDInt_x_CFD_low)
println(PSDInt_x_NaluWave_low)

bar_x_wave = bar((1,abs.(PSDInt_x_EXP_wave)),label="Experiment",legend=:outertopright)
bar!((2,abs.(PSDInt_x_CFD_wave)),label="StarCCM+")
bar!((3,abs.(PSDInt_x_NaluWave_wave)),label="ExaWind")
hline!(ones(length(PSDInt_x_EXP_wave))*(abs.(PSDInt_x_EXP_wave)*plus_range_wave_frequency),color=:red,linestyle=:dash,label="+"*wave_range_str*"% from Experiment")
hline!(ones(length(PSDInt_x_EXP_wave))*(abs.(PSDInt_x_EXP_wave)*minus_range_wave_frequency),color=:black,linestyle=:dash,label="-"*wave_range_str*"% from Experiment")
ylabel!("Wave-frequency PSD Integral of \$x [m^2]\$")
savefig(bar_x_wave,"bar_x_PSD_wave_frequency.pdf")

println("x_wave")
println(PSDInt_x_EXP_wave)
println(PSDInt_x_CFD_wave)
println(PSDInt_x_NaluWave_wave)

bar_z = bar((1,abs.(PSDInt_z_EXP_low)),label="Experiment",legend=:outertopright)
bar!((2,abs.(PSDInt_z_CFD_low)),label="StarCCM+")
bar!((3,abs.(PSDInt_z_NaluWave_low)),label="ExaWind")
hline!(ones(length(PSDInt_z_EXP_low))*(abs.(PSDInt_z_EXP_low)*plus_range_low_frequency),color=:red,linestyle=:dash,label="+"*low_range_str*"% from Experiment")
hline!(ones(length(PSDInt_z_EXP_low))*(abs.(PSDInt_z_EXP_low)*minus_range_low_frequency),color=:black,linestyle=:dash,label="-"*low_range_str*"% from Experiment")
ylabel!("Low-Frequency PSD Integral of \$z [m^2]\$")
savefig(bar_z,"bar_z_PSD_low_frequency.pdf")

println("z_low")
println(PSDInt_z_EXP_low)
println(PSDInt_z_CFD_low)
println(PSDInt_z_NaluWave_low)

bar_z_wave = bar((1,abs.(PSDInt_z_EXP_wave)),label="Experiment",legend=:outertopright)
bar!((2,abs.(PSDInt_z_CFD_wave)),label="StarCCM+")
bar!((3,abs.(PSDInt_z_NaluWave_wave)),label="ExaWind")
hline!(ones(length(PSDInt_z_EXP_wave))*(abs.(PSDInt_z_EXP_wave)*plus_range_wave_frequency),color=:red,linestyle=:dash,label="+"*wave_range_str*"% from Experiment")
hline!(ones(length(PSDInt_z_EXP_wave))*(abs.(PSDInt_z_EXP_wave)*minus_range_wave_frequency),color=:black,linestyle=:dash,label="-"*wave_range_str*"% from Experiment")
ylabel!("Wave-frequency PSD Integral of \$z [m^2]\$")
savefig(bar_z_wave,"bar_z_PSD_wave_frequency.pdf")

println("z_wave")
println(PSDInt_z_EXP_wave)
println(PSDInt_z_CFD_wave)
println(PSDInt_z_NaluWave_wave)

bar_thy = bar((1,abs.(PSDInt_thy_EXP_low)),label="Experiment",legend=:outertopright)
bar!((2,abs.(PSDInt_thy_CFD_low)),label="StarCCM+")
bar!((3,abs.(PSDInt_thy_NaluWave_low)),label="ExaWind")
hline!(ones(length(PSDInt_thy_EXP_low))*(abs.(PSDInt_thy_EXP_low)*plus_range_low_frequency),color=:red,linestyle=:dash,label="+"*low_range_str*"% from Experiment")
hline!(ones(length(PSDInt_thy_EXP_low))*(abs.(PSDInt_thy_EXP_low)*minus_range_low_frequency),color=:black,linestyle=:dash,label="-"*low_range_str*"% from Experiment")
ylabel!("Low-frequency PSD Integral of " * L"\theta_y [^o^2]")
savefig(bar_thy,"bar_thy_PSD_low_frequency.pdf")

println("thy_low")
println(PSDInt_thy_EXP_low)
println(PSDInt_thy_CFD_low)
println(PSDInt_thy_NaluWave_low)

bar_thy_wave = bar((1,abs.(PSDInt_thy_EXP_wave)),label="Experiment",legend=:outertopright)
bar!((2,abs.(PSDInt_thy_CFD_wave)),label="StarCCM+")
bar!((3,abs.(PSDInt_thy_NaluWave_wave)),label="ExaWind")
hline!(ones(length(PSDInt_thy_EXP_wave))*(abs.(PSDInt_thy_EXP_wave)*plus_range_wave_frequency),color=:red,linestyle=:dash,label="+"*wave_range_str*"% from Experiment")
hline!(ones(length(PSDInt_thy_EXP_wave))*(abs.(PSDInt_thy_EXP_wave)*minus_range_wave_frequency),color=:black,linestyle=:dash,label="-"*wave_range_str*"% from Experiment")
ylabel!("Wave-frequency PSD Integral of " * L"\theta_y [^o^2]")
savefig(bar_thy_wave,"bar_thy_PSD_wave_frequency.pdf")

println("thy_wave")
println(PSDInt_thy_EXP_wave)
println(PSDInt_thy_CFD_wave)
println(PSDInt_thy_NaluWave_wave)