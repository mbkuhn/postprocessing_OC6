using DelimitedFiles
using Interpolations
using DataFrames
using LaTeXStrings
using StatsBase
using MAT
using Measurements
include("FFTAnalysis.jl")
include("rangedTrapz.jl")
include("get_evenly_spaced_time_series.jl")

# Set relevant frequency ranges based on MARIN's set-up
lowfrequencies_range = [0.005 0.05]
wavefrequencies_range = [0.05519 0.1345]

# Set relevant Milestone's +- validation ranges
plus_range_wave_frequency = 1.15
minus_range_wave_frequency = 0.85
plus_range_low_frequency = 1.25
minus_range_low_frequency = 0.75

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

file_exp = "/Users/mkuhn/testruns_data/OC6_PhaseIa_dataset_and_more/OC6_LC33_LC53_CFD_scripts_data/LC33_FixedStructure/ZXP0_LC33_full.txt"
file_cfd_ref = "/Users/mkuhn/testruns_data/OC6_PhaseIa_dataset_and_more/OC6_LC33_LC53_CFD_scripts_data/LC33_FixedStructure/NREL_LC33_CFD.txt"
file_ofast_exp = "/Users/mkuhn/testruns_data/OC6_PhaseIa_dataset_and_more/OC6_LC33_LC53_CFD_scripts_data/LC33_FixedStructure/NREL_LC33_FAST_EXPWave.txt"
file_ofast_cfd = "/Users/mkuhn/testruns_data/OC6_PhaseIa_dataset_and_more/OC6_LC33_LC53_CFD_scripts_data/LC33_FixedStructure/NREL_LC33_FAST_CFDWave.txt"
file_exawind_pfx = "/Users/mkuhn/testruns_data/OC5_semisubmersible/HOS-NWT/pointwise_mesh/floater_forces"
file_exawind_sfx = ".dat"
file_exawind_labels = ["1","2","3","4","5"]


file_amrwind = ""
file_hosnwt = ""

if (!isempty(file_amrwind))
    rows_to_skip = 49
    # Data from probes for JONSWAP spectrum from HOS-NWT simulation (only needed if the comparison with wave elevation is needed)
    HOS_NWT_Case4_data = readdlm(file_hosnwt)
    HOS_NWT_Case4_data = HOS_NWT_Case4_data[rows_to_skip:end,1:3]
    HOS_NWT_Case4_data = Float64.(HOS_NWT_Case4_data)
    HOS_NWT_Case4_time = HOS_NWT_Case4_data[:,1]
    HOS_NWT_Case4_wave = HOS_NWT_Case4_data[:,3]

    # AMR-Wind Standalone (only needed if the comparison with wave elevation is needed)
    amrwind_data = readdlm(file_amrwind)
    amrwind_data = Float64.(amrwind_data)
    amrwind_time = amrwind_data[:,1]
    amrwind_wave = amrwind_data[:,2]
end

####################### FORCES ##########################################################################################################
# Data from OC6 Experiment
ZXP0_data = readdlm(file_exp)
ZXP0_data = ZXP0_data[2001:end,:]
ZXP0_data = Float64.(ZXP0_data)
ZXP0_data[:,1] = ZXP0_data[:,1].-200 # "." dot is to subtract elementwise for the column
EXP_time = ZXP0_data[:,1]
EXP_wave = ZXP0_data[:,2]
EXP_Fx = ZXP0_data[:,3] 
EXP_Fz = ZXP0_data[:,5] 
EXP_My = ZXP0_data[:,7]

# Data from StarCCM+ simulation
CFD_data = readdlm(file_cfd_ref)
CFD_data = Float64.(CFD_data[2:end,:])
CFD_time = CFD_data[:,1]
CFD_wave = CFD_data[:,2]
CFD_Fx = CFD_data[:,3] 
CFD_Fz = CFD_data[:,5] 
CFD_My = CFD_data[:,7]

# Data from OpenFAST Simulation using Experimental Waves as input
EXPWAVE_data = readdlm(file_ofast_exp)
EXPWAVE_data = EXPWAVE_data[2001:end,:] 
EXPWAVE_data = Float64.(EXPWAVE_data)
EXPWAVE_data[:,1] = EXPWAVE_data[:,1].-200;
FAST_EXPWave_time = EXPWAVE_data[:,1] 
FAST_EXPWave_wave = EXPWAVE_data[:,2]
FAST_EXPWave_Fx = EXPWAVE_data[:,3] 
FAST_EXPWave_Fz = EXPWAVE_data[:,5] 
FAST_EXPWave_My = EXPWAVE_data[:,7]

# Data from OpenFAST Simulation using StarCCM+ Waves as input
CFDWAVE_data = readdlm(file_ofast_cfd)
CFDWAVE_data = Float64.(CFDWAVE_data[2:end,:])
FASTCFDWave_time = CFDWAVE_data[:,1]
FASTCFDWave_wave = CFDWAVE_data[:,2]
FASTCFDWave_Fx = CFDWAVE_data[:,3] 
FASTCFDWave_Fz = CFDWAVE_data[:,5] 
FASTCFDWave_My = CFDWAVE_data[:,7]

# NaluWind DATA

# Conversion between full scale and model scale
# L_p = 50 L_m
# T_p = sqrt(50)T_m
# F_p = L^3 F_m
# M_y = L^4 M_y
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
    #global Nalu_data = Nalu_data[1:last_index,:]
    global Nalu_data = vcat(Nalu_data[1:last_index,:],Nalu_data_tmp)
end
# converting to full (prototype) scale
NaluWAVE_data = Float64.(Nalu_data[2:end,:])
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
# Fx
Fx_Nalu = NaluWAVE_data[:,2].+NaluWAVE_data[:,5]
_ , Nalu_EvenlySpaced_Fx = get_evenly_spaced_time_series(NaluWAVE_data[:,1],Fx_Nalu)
# Fz
Fz_Nalu = NaluWAVE_data[:,4].+NaluWAVE_data[:,7]
_ , Nalu_EvenlySpaced_Fz = get_evenly_spaced_time_series(NaluWAVE_data[:,1],Fz_Nalu)
# My
_ , Nalu_EvenlySpaced_My = get_evenly_spaced_time_series(NaluWAVE_data[:,1],NaluWAVE_data[:,9])

# Convert time to full scale
NaluWave_time = time_evenly_spaced.*sqrt(rescaling)

# In case one needs to discard any transient time data
#t_final = NaluWave_time[1]./sqrt(rescaling) + 0.027*(NaluWave_time[end]./sqrt(rescaling)) + 545.0 # seconds (model scale)
#t_skip = NaluWave_time[1]./sqrt(rescaling) + 0.027*(NaluWave_time[end]./sqrt(rescaling)) # seconds (model scale) # if skipping some times is needed 
t_skip = 20.
t_final = t_skip + 545.
println("t_skip = ",t_skip)
indices = findall(NaluWave_time->(NaluWave_time>=(t_skip*sqrt(rescaling)) && NaluWave_time<=(t_final*sqrt(rescaling))),NaluWave_time)

NaluWave_time = NaluWave_time[indices]
println("Initial time (Model Scale) = ", NaluWave_time[1]/sqrt(rescaling))
println("Final time (Model Scale) = ", NaluWave_time[end]/sqrt(rescaling))

# Convert Fx, Fz, and My to full scale after skipping transient portion of array
NaluWave_Fx = (Nalu_EvenlySpaced_Fx[indices]).*(rescaling^3)
NaluWave_Fz = (Nalu_EvenlySpaced_Fz[indices]).*(rescaling^3)
NaluWave_My = (Nalu_EvenlySpaced_My[indices]).*(rescaling^4)

# this was to check the difference coming from the density of water
#NaluWave_Fx *= 1.025 
#NaluWave_Fz *= 1.025 
#NaluWave_My *= 1.025 

stepRange = collect(3000:108000) # We are skipping 5 minutes (300 seconds) of full scale time (roughly 42 seconds of model scale time), out of a 10800 seconds which are 180 minutes (3 hours) of full scale time (26 minutes of model scale)

stepRangeNalu = collect(1:length(NaluWave_time)) # include all of the Nalu time

# Plot Waves Time series (Model Scale)
if (!isempty(file_amrwind))
    Time_series_plot = plot(EXP_time./sqrt(rescaling),EXP_wave./(rescaling),label="Experiment",xlims=(0,t_final),legend=:outertop)
    plot!(HOS_NWT_Case4_time,HOS_NWT_Case4_wave,label="HOS-NWT",linestyle=:dashdotdot)
    plot!(amrwind_time,amrwind_wave,label="AMR-Wind Standalone",linestyle=:dash)
    xlabel!("Time \$[sec]\$ - Model Scale")
    ylabel!("Wave Elevation  \$[m]\$ - Model Scale")
    savefig(Time_series_plot,"wave_elevation_time_series.pdf")
end

# Mean Surge - Unrelated to PSD, just the mean value of the surge force across times
mean_NaluFx = mean(NaluWave_Fx[stepRangeNalu])
mean_EXP_Fx = mean(EXP_Fx[stepRange])
mean_CFD_Fx = mean(CFD_Fx[stepRange])
mean_FastEXP_Fx = mean(FAST_EXPWave_Fx[stepRange])
mean_FastCFD_Fx = mean(FASTCFDWave_Fx[stepRange])
println("Nalu Mean Surge Force = ", mean_NaluFx)
println("Exp Mean Surge Force = ", mean_EXP_Fx)
println("CFD Mean Surge Force = ", mean_CFD_Fx)
println("OPEN_CFD Mean Surge Force = ", mean_FastCFD_Fx)
println("OPEN_EXP Mean Surge Force = ", mean_FastEXP_Fx)
println("Ratio = ", mean_NaluFx/mean_EXP_Fx)
mean_surge_plot = bar((1,mean_EXP_Fx),label="Experiment",legend=:outertopright)
bar!((2,mean_CFD_Fx),label="StarCCM+")
bar!((3,mean_FastEXP_Fx),label="OpenFAST (Experiment Input Waves)")
bar!((4,mean_FastCFD_Fx),label="OpenFAST (StarCCM+ Input Waves)")
bar!((5,mean_NaluFx),label="ExaWind")
hline!(ones(length(mean_EXP_Fx))*(abs.(mean_EXP_Fx)*plus_range_low_frequency),color=:red,linestyle=:dash,label="+25% from Experiment")
hline!(ones(length(mean_EXP_Fx))*(abs.(mean_EXP_Fx)*minus_range_low_frequency),color=:black,linestyle=:dash,label="-25% from Experiment")
ylabel!("Mean Surge Force \$[N]\$")
savefig(mean_surge_plot,"Mean_Surge.pdf")

### **** PSD COMPUTATIONS AND PLOTS ******

# first one is for plotting purposes
# wave EXP
f, P1, PSD, fdouble, P2 = FFTAnalysis(EXP_time[stepRange],EXP_wave[stepRange],true,false,smoothing_type,smoothing_parameter)
wave_plot = plot(f[1:end],abs.(PSD[1:end]),yscale=:log10,label="Experiment",xlims=(0,0.4),legend=:bottomright)
f, P1, PSD, fdouble, P2 = FFTAnalysis(EXP_time[stepRange],EXP_wave[stepRange],false,false,~,~)
PSDInt_wave_EXP_low  = rangedTrapz(f,abs.(PSD),lowfrequencies_range)
PSDInt_wave_EXP_wave = rangedTrapz(f,abs.(PSD),wavefrequencies_range)

# wave CFD
f, P1, PSD, fdouble, P2 = FFTAnalysis(CFD_time[stepRange],CFD_wave[stepRange],true,false,smoothing_type,smoothing_parameter)
plot!(f[1:end],abs.(PSD[1:end]),yscale=:log10,label="StarCCM+")
f, P1, PSD, fdouble, P2 = FFTAnalysis(CFD_time[stepRange],CFD_wave[stepRange],false,false,~,~)
PSDInt_wave_CFD_low  = rangedTrapz(f,abs.(PSD),lowfrequencies_range)
PSDInt_wave_CFD_wave = rangedTrapz(f,abs.(PSD),wavefrequencies_range)

# wave FAST_EXPWave
f, P1, PSD, fdouble, P2 = FFTAnalysis(FAST_EXPWave_time[stepRange],FAST_EXPWave_wave[stepRange],true,false,smoothing_type,smoothing_parameter)
f, P1, PSD, fdouble, P2 = FFTAnalysis(FAST_EXPWave_time[stepRange],FAST_EXPWave_wave[stepRange],false,false,~,~)
PSDInt_wave_FAST_EXPWAVE_low  = rangedTrapz(f,abs.(PSD),lowfrequencies_range)
PSDInt_wave_FAST_EXPWAVE_wave = rangedTrapz(f,abs.(PSD),wavefrequencies_range)

# wave FASTCFDWave
f, P1, PSD, fdouble, P2 = FFTAnalysis(FASTCFDWave_time[stepRange],FASTCFDWave_wave[stepRange],true,false,smoothing_type,smoothing_parameter)
f, P1, PSD, fdouble, P2 = FFTAnalysis(FASTCFDWave_time[stepRange],FASTCFDWave_wave[stepRange],false,false,~,~)
PSDInt_wave_FAST_CFDWAVE_low  = rangedTrapz(f,abs.(PSD),lowfrequencies_range)
PSDInt_wave_FAST_CFDWAVE_wave = rangedTrapz(f,abs.(PSD),wavefrequencies_range)

# wave AMR-Wind # Not part of the Milestone: only to validate AMR-Wind codebase
if (!isempty(file_amrwind))
    AMRWind_WAVE_time_interp = Interpolations.deduplicate_knots!(amrwind_time)
    itp_waves_amrwind =  interpolate((AMRWind_WAVE_time_interp,), amrwind_wave, Gridded(Linear()))
    amrwind_wave_evenly_spaced = itp_waves_amrwind.(time_evenly_spaced)
    AMRWind_time = time_evenly_spaced.*sqrt(rescaling)
    AMRWind_Wave = amrwind_wave_evenly_spaced.*rescaling

    # adding also the curve with just 640 seconds
    indices_amrwind = findall(AMRWind_time->(AMRWind_time>=(t_skip.*sqrt(rescaling)) && AMRWind_time<=(t_final.*sqrt(rescaling))),AMRWind_time)

    #f, P1, PSD, fdouble, P2 = FFTAnalysis(amrwind_time.*sqrt(rescaling),amrwind_wave.*rescaling,true,false,smoothing_type,smoothing_parameter)
    f, P1, PSD, fdouble, P2 = FFTAnalysis(AMRWind_time,AMRWind_Wave,true,false,smoothing_type,smoothing_parameter)
    f, P1, PSD, fdouble, P2 = FFTAnalysis(AMRWind_time[indices_amrwind],AMRWind_Wave[indices_amrwind],true,false,smoothing_type,smoothing_parameter)
    plot!(f[1:end],abs.(PSD[1:end]),yscale=:log10,label="AMR-Wind Standalone",linestyle=:dash)
    vspan!([0.05519,0.1345]; alpha = 0.3,label="Frequency Range")
    vspan!([0.005,0.05]; alpha = 0.3,label="Low Frequency Range")
    xlabel!("Frequency [Hz]")
    ylabel!("PSD of Incident Wave \$[m^2/Hz]\$")
    savefig(wave_plot,"wave_PSD.pdf")
    f, P1, PSD, fdouble, P2 = FFTAnalysis(AMRWind_time,AMRWind_Wave,false,false,~,~)
    PSDInt_wave_amrwind_low  = rangedTrapz(f,abs.(PSD),lowfrequencies_range)
    PSDInt_wave_amrwind_wave = rangedTrapz(f,abs.(PSD),wavefrequencies_range)
end

# first one is for plotting purposes
# Fx EXP
f, P1, PSD, fdouble, P2 = FFTAnalysis(EXP_time[stepRange],EXP_Fx[stepRange],true,false,smoothing_type,smoothing_parameter)
Fx_plot = plot(f[1:end],abs.(PSD[1:end]),yscale=:log10,label="Experiment",xlims=(0,0.4),ylims=(1e10,1e16))
f, P1, PSD, fdouble, P2 = FFTAnalysis(EXP_time[stepRange],EXP_Fx[stepRange],false,false,~,~)
PSDInt_Fx_EXP_low  = rangedTrapz(f,abs.(PSD),lowfrequencies_range)
PSDInt_Fx_EXP_wave = rangedTrapz(f,abs.(PSD),wavefrequencies_range)

# Fx CFD
f, P1, PSD, fdouble, P2 = FFTAnalysis(CFD_time[stepRange],CFD_Fx[stepRange],true,false,smoothing_type,smoothing_parameter)
plot!(f[1:end],abs.(PSD[1:end]),yscale=:log10,label="StarCCM+")
f, P1, PSD, fdouble, P2 = FFTAnalysis(CFD_time[stepRange],CFD_Fx[stepRange],false,false,~,~)
PSDInt_Fx_CFD_low  = rangedTrapz(f,abs.(PSD),lowfrequencies_range)
PSDInt_Fx_CFD_wave = rangedTrapz(f,abs.(PSD),wavefrequencies_range)

# Fx FAST_EXPWave
f, P1, PSD, fdouble, P2 = FFTAnalysis(FAST_EXPWave_time[stepRange],FAST_EXPWave_Fx[stepRange],true,false,smoothing_type,smoothing_parameter)
plot!(f[1:end],abs.(PSD[1:end]),yscale=:log10,label="OpenFAST (Experiment Input Waves)")
f, P1, PSD, fdouble, P2 = FFTAnalysis(FAST_EXPWave_time[stepRange],FAST_EXPWave_Fx[stepRange],false,false,~,~)
PSDInt_Fx_FAST_EXPWAVE_low  = rangedTrapz(f,abs.(PSD),lowfrequencies_range)
PSDInt_Fx_FAST_EXPWAVE_wave = rangedTrapz(f,abs.(PSD),wavefrequencies_range)

# Fx FASTCFDWave
f, P1, PSD, fdouble, P2 = FFTAnalysis(FASTCFDWave_time[stepRange],FASTCFDWave_Fx[stepRange],true,false,smoothing_type,smoothing_parameter)
plot!(f[1:end],abs.(PSD[1:end]),yscale=:log10,label="OpenFAST (StarCCM+ Input Waves)")
f, P1, PSD, fdouble, P2 = FFTAnalysis(FASTCFDWave_time[stepRange],FASTCFDWave_Fx[stepRange],false,false,~,~)
PSDInt_Fx_FAST_CFDWAVE_low  = rangedTrapz(f,abs.(PSD),lowfrequencies_range)
PSDInt_Fx_FAST_CFDWAVE_wave = rangedTrapz(f,abs.(PSD),wavefrequencies_range)

# NALU Fx
f, P1, PSD, fdouble, P2 = FFTAnalysis(NaluWave_time[stepRangeNalu],NaluWave_Fx[stepRangeNalu],true,false,smoothing_type,smoothing_parameter)
plot!(f[1:end],abs.(PSD[1:end]),yscale=:log10,label="ExaWind")
vspan!([0.05519,0.1345]; alpha = 0.3,label="Frequency Range")
vspan!([0.005,0.05]; alpha = 0.3,label="Low Frequency Range")
xlabel!("Frequency [Hz]")
ylabel!("Global \$F_x [N^2/Hz]\$")
savefig(Fx_plot,"Fx_PSD.pdf")
f, P1, PSD, fdouble, P2 = FFTAnalysis(NaluWave_time[stepRangeNalu],NaluWave_Fx[stepRangeNalu],false,false,~,~)
PSDInt_Fx_NaluWave_low  = rangedTrapz(f,abs.(PSD),lowfrequencies_range)
PSDInt_Fx_NaluWave_wave = rangedTrapz(f,abs.(PSD),wavefrequencies_range) # these are just to get the bar plots

# first one is for plotting purposes
# Fz EXP
f, P1, PSD, fdouble, P2 = FFTAnalysis(EXP_time[stepRange],EXP_Fz[stepRange],true,false,smoothing_type,smoothing_parameter)
Fz_plot = plot(f[1:end],abs.(PSD[1:end]),yscale=:log10,label="Experiment",xlims=(0,0.4))
f, P1, PSD, fdouble, P2 = FFTAnalysis(EXP_time[stepRange],EXP_Fz[stepRange],false,false,~,~)
PSDInt_Fz_EXP_low  = rangedTrapz(f,abs.(PSD),lowfrequencies_range)
PSDInt_Fz_EXP_wave = rangedTrapz(f,abs.(PSD),wavefrequencies_range)

# Fz CFD
f, P1, PSD, fdouble, P2 = FFTAnalysis(CFD_time[stepRange],CFD_Fz[stepRange],true,false,smoothing_type,smoothing_parameter)
plot!(f[1:end],abs.(PSD[1:end]),yscale=:log10,label="StarCCM+")
f, P1, PSD, fdouble, P2 = FFTAnalysis(CFD_time[stepRange],CFD_Fz[stepRange],false,false,~,~)
PSDInt_Fz_CFD_low  = rangedTrapz(f,abs.(PSD),lowfrequencies_range)
PSDInt_Fz_CFD_wave = rangedTrapz(f,abs.(PSD),wavefrequencies_range)

# Fz FAST_EXPWave
f, P1, PSD, fdouble, P2 = FFTAnalysis(FAST_EXPWave_time[stepRange],FAST_EXPWave_Fz[stepRange],true,false,smoothing_type,smoothing_parameter)
plot!(f[1:end],abs.(PSD[1:end]),yscale=:log10,label="OpenFAST (Experiment Input Waves)")
f, P1, PSD, fdouble, P2 = FFTAnalysis(FAST_EXPWave_time[stepRange],FAST_EXPWave_Fz[stepRange],false,false,~,~)
PSDInt_Fz_FAST_EXPWAVE_low  = rangedTrapz(f,abs.(PSD),lowfrequencies_range)
PSDInt_Fz_FAST_EXPWAVE_wave = rangedTrapz(f,abs.(PSD),wavefrequencies_range)

# Fz FASTCFDWave
f, P1, PSD, fdouble, P2 = FFTAnalysis(FASTCFDWave_time[stepRange],FASTCFDWave_Fz[stepRange],true,false,smoothing_type,smoothing_parameter)
plot!(f[1:end],abs.(PSD[1:end]),yscale=:log10,label="OpenFAST (StarCCM+ Input Waves)")
f, P1, PSD, fdouble, P2 = FFTAnalysis(FASTCFDWave_time[stepRange],FASTCFDWave_Fz[stepRange],false,false,~,~)
PSDInt_Fz_FAST_CFDWAVE_low  = rangedTrapz(f,abs.(PSD),lowfrequencies_range)
PSDInt_Fz_FAST_CFDWAVE_wave = rangedTrapz(f,abs.(PSD),wavefrequencies_range)

# NALU Fz
f, P1, PSD, fdouble, P2 = FFTAnalysis(NaluWave_time[stepRangeNalu],NaluWave_Fz[stepRangeNalu],true,false,smoothing_type,smoothing_parameter)
plot!(f[1:end],abs.(PSD[1:end]),yscale=:log10,label="ExaWind") # need to investigate about this truncation
vspan!([0.05519,0.1345]; alpha = 0.3,label="Frequency Range")
vspan!([0.005,0.05]; alpha = 0.3,label="Low Frequency Range")
xlabel!("Frequency [Hz]")
ylabel!("Global \$F_z [N^2/Hz]\$")
savefig(Fz_plot,"Fz_PSD.pdf")
f, P1, PSD, fdouble, P2 = FFTAnalysis(NaluWave_time[stepRangeNalu],NaluWave_Fz[stepRangeNalu],false,false,~,~)
PSDInt_Fz_NaluWave_low  = rangedTrapz(f,abs.(PSD),lowfrequencies_range)
PSDInt_Fz_NaluWave_wave = rangedTrapz(f,abs.(PSD),wavefrequencies_range) # these are just to get the bar plots

# first one is for plotting purposes
# My EXP
f, P1, PSD, fdouble, P2 = FFTAnalysis(EXP_time[stepRange],EXP_My[stepRange],true,false,smoothing_type,smoothing_parameter)
My_plot = plot(f[1:end],abs.(PSD[1:end]),yscale=:log10,label="Experiment",xlims=(0,0.4),ylims=(1e13,1e18),legend=:topright) # in the paper, the xlim is [0,0.4], in the code is [0,0.14]
f, P1, PSD, fdouble, P2 = FFTAnalysis(EXP_time[stepRange],EXP_My[stepRange],true,false,smoothing_type,smoothing_parameter)
PSDInt_My_EXP_low  = rangedTrapz(f,abs.(PSD),lowfrequencies_range)
PSDInt_My_EXP_wave = rangedTrapz(f,abs.(PSD),wavefrequencies_range)

# My CFD
f, P1, PSD, fdouble, P2 = FFTAnalysis(CFD_time[stepRange],CFD_My[stepRange],true,false,smoothing_type,smoothing_parameter)
plot!(f[1:end],abs.(PSD[1:end]),yscale=:log10,label="StarCCM+")
f, P1, PSD, fdouble, P2 = FFTAnalysis(CFD_time[stepRange],CFD_My[stepRange],true,false,smoothing_type,smoothing_parameter)
PSDInt_My_CFD_low  = rangedTrapz(f,abs.(PSD),lowfrequencies_range)
PSDInt_My_CFD_wave = rangedTrapz(f,abs.(PSD),wavefrequencies_range)

# My FAST_EXPWave
f, P1, PSD, fdouble, P2 = FFTAnalysis(FAST_EXPWave_time[stepRange],FAST_EXPWave_My[stepRange],true,false,smoothing_type,smoothing_parameter)
plot!(f[1:end],abs.(PSD[1:end]),yscale=:log10,label="OpenFAST (Experiment Input Waves)")
f, P1, PSD, fdouble, P2 = FFTAnalysis(FAST_EXPWave_time[stepRange],FAST_EXPWave_My[stepRange],true,false,smoothing_type,smoothing_parameter)
PSDInt_My_FAST_EXPWAVE_low  = rangedTrapz(f,abs.(PSD),lowfrequencies_range)
PSDInt_My_FAST_EXPWAVE_wave = rangedTrapz(f,abs.(PSD),wavefrequencies_range)

# My FASTCFDWave
f, P1, PSD, fdouble, P2 = FFTAnalysis(FASTCFDWave_time[stepRange],FASTCFDWave_My[stepRange],true,false,smoothing_type,smoothing_parameter)
plot!(f[1:end],abs.(PSD[1:end]),yscale=:log10,label="OpenFAST (StarCCM+ Input Waves)")
f, P1, PSD, fdouble, P2 = FFTAnalysis(FASTCFDWave_time[stepRange],FASTCFDWave_My[stepRange],true,false,smoothing_type,smoothing_parameter)
PSDInt_My_FAST_CFDWAVE_low  = rangedTrapz(f,abs.(PSD),lowfrequencies_range)
PSDInt_My_FAST_CFDWAVE_wave = rangedTrapz(f,abs.(PSD),wavefrequencies_range)

# NALU
f, P1, PSD, fdouble, P2 = FFTAnalysis(NaluWave_time[stepRangeNalu],NaluWave_My[stepRangeNalu],true,false,smoothing_type,smoothing_parameter)
plot!(f[1:end],abs.(PSD[1:end]),yscale=:log10,label="ExaWind") # need to investigate about this truncation
vspan!([0.05519,0.1345]; alpha = 0.3,label="Frequency Range")
vspan!([0.005,0.05]; alpha = 0.3,label="Low Frequency Range")
xlabel!("Frequency [Hz]")
ylabel!("Global \$M_y [(Nm)^2/Hz]\$")
savefig(My_plot,"My_PSD.pdf")
f, P1, PSD, fdouble, P2 = FFTAnalysis(NaluWave_time[stepRangeNalu],NaluWave_My[stepRangeNalu],true,false,smoothing_type,smoothing_parameter) # HERE
PSDInt_My_NaluWave_low  = rangedTrapz(f,abs.(PSD),lowfrequencies_range)
PSDInt_My_NaluWave_wave = rangedTrapz(f,abs.(PSD),wavefrequencies_range) # these are just to get the bar plots

####### PSD INTEGRAL PLOTS (Bar plots) ##################################################
bar_Fx = bar((1,abs.(PSDInt_Fx_EXP_low)),label="Experiment",legend=:outertopright)
bar!((2,abs.(PSDInt_Fx_CFD_low)),label="StarCCM+")
bar!((3,abs.(PSDInt_Fx_FAST_EXPWAVE_low)),label="OpenFAST (Experiment Input Waves)")
bar!((4,abs.(PSDInt_Fx_FAST_CFDWAVE_low)),label="OpenFAST (StarCCM+ Input Waves)")
bar!((5,abs.(PSDInt_Fx_NaluWave_low)),label="ExaWind")
hline!(ones(length(PSDInt_Fx_EXP_low))*(abs.(PSDInt_Fx_EXP_low)*plus_range_low_frequency),color=:red,linestyle=:dash,label="+25% from Experiment")
hline!(ones(length(PSDInt_Fx_EXP_low))*(abs.(PSDInt_Fx_EXP_low)*minus_range_low_frequency),color=:black,linestyle=:dash,label="-25% from Experiment")
ylabel!("Low-Frequency PSD Integral of \$F_x [N^2]\$")
savefig(bar_Fx,"bar_Fx_PSD_low_frequency.pdf")

println("Fx_low")
println(PSDInt_Fx_EXP_low)
println(PSDInt_Fx_CFD_low)
println(PSDInt_Fx_FAST_EXPWAVE_low)
println(PSDInt_Fx_FAST_CFDWAVE_low)
println(PSDInt_Fx_NaluWave_low)

bar_Fx_wave = bar((1,abs.(PSDInt_Fx_EXP_wave)),label="Experiment",legend=:outertopright)
bar!((2,abs.(PSDInt_Fx_CFD_wave)),label="StarCCM+")
bar!((3,abs.(PSDInt_Fx_FAST_EXPWAVE_wave)),label="OpenFAST (Experiment Input Waves)")
bar!((4,abs.(PSDInt_Fx_FAST_CFDWAVE_wave)),label="OpenFAST (StarCCM+ Input Waves)")
bar!((5,abs.(PSDInt_Fx_NaluWave_wave)),label="ExaWind")
hline!(ones(length(PSDInt_Fx_EXP_wave))*(abs.(PSDInt_Fx_EXP_wave)*plus_range_wave_frequency),color=:red,linestyle=:dash,label="+15% from Experiment")
hline!(ones(length(PSDInt_Fx_EXP_wave))*(abs.(PSDInt_Fx_EXP_wave)*minus_range_wave_frequency),color=:black,linestyle=:dash,label="-15% from Experiment")
ylabel!("Wave-frequency PSD Integral of \$F_x [N^2]\$")
savefig(bar_Fx_wave,"bar_Fx_PSD_wave_frequency.pdf")

println("Fx_wave")
println(PSDInt_Fx_EXP_wave)
println(PSDInt_Fx_CFD_wave)
println(PSDInt_Fx_FAST_EXPWAVE_wave)
println(PSDInt_Fx_FAST_CFDWAVE_wave)
println((PSDInt_Fx_NaluWave_wave))

bar_Fz = bar((1,abs.(PSDInt_Fz_EXP_low)),label="Experiment",legend=:outertopright)
bar!((2,abs.(PSDInt_Fz_CFD_low)),label="StarCCM+")
bar!((3,abs.(PSDInt_Fz_FAST_EXPWAVE_low)),label="OpenFAST (Experiment Input Waves)")
bar!((4,abs.(PSDInt_Fz_FAST_CFDWAVE_low)),label="OpenFAST (StarCCM+ Input Waves)")
bar!((5,abs.(PSDInt_Fz_NaluWave_low)),label="ExaWind")
hline!(ones(length(PSDInt_Fz_EXP_low))*(abs.(PSDInt_Fz_EXP_low)*plus_range_low_frequency),color=:red,linestyle=:dash,label="+25% from Experiment")
hline!(ones(length(PSDInt_Fz_EXP_low))*(abs.(PSDInt_Fz_EXP_low)*minus_range_low_frequency),color=:black,linestyle=:dash,label="-25% from Experiment")
ylabel!("Low-Frequency PSD Integral of \$F_z [N^2]\$")
savefig(bar_Fz,"bar_Fz_PSD_low_frequency.pdf")

println("Fz_low")
println(PSDInt_Fz_EXP_low)
println(PSDInt_Fz_CFD_low)
println(PSDInt_Fz_FAST_EXPWAVE_low)
println(PSDInt_Fz_FAST_CFDWAVE_low)

bar_Fz_wave = bar((1,abs.(PSDInt_Fz_EXP_wave)),label="Experiment",legend=:outertopright)
bar!((2,abs.(PSDInt_Fz_CFD_wave)),label="StarCCM+")
bar!((3,abs.(PSDInt_Fz_FAST_EXPWAVE_wave)),label="OpenFAST (Experiment Input Waves)")
bar!((4,abs.(PSDInt_Fz_FAST_CFDWAVE_wave)),label="OpenFAST (StarCCM+ Input Waves)")
bar!((5,abs.(PSDInt_Fz_NaluWave_wave)),label="ExaWind")
hline!(ones(length(PSDInt_Fz_EXP_wave))*(abs.(PSDInt_Fz_EXP_wave)*plus_range_wave_frequency),color=:red,linestyle=:dash,label="+15% from Experiment")
hline!(ones(length(PSDInt_Fz_EXP_wave))*(abs.(PSDInt_Fz_EXP_wave)*minus_range_wave_frequency),color=:black,linestyle=:dash,label="-15% from Experiment")
ylabel!("Wave-frequency PSD Integral of \$F_z [N^2]\$")
savefig(bar_Fz_wave,"bar_Fz_PSD_wave_frequency.pdf")

println("Fz_wave")
println(PSDInt_Fz_EXP_wave)
println(PSDInt_Fz_CFD_wave)
println(PSDInt_Fz_FAST_EXPWAVE_wave)
println(PSDInt_Fz_FAST_CFDWAVE_wave)

bar_My = bar((1,abs.(PSDInt_My_EXP_low)),label="Experiment",legend=:outertopright)
bar!((2,abs.(PSDInt_My_CFD_low)),label="StarCCM+")
bar!((3,abs.(PSDInt_My_FAST_EXPWAVE_low)),label="OpenFAST (Experiment Input Waves)")
bar!((4,abs.(PSDInt_My_FAST_CFDWAVE_low)),label="OpenFAST (StarCCM+ Input Waves)")
bar!((5,abs.(PSDInt_My_NaluWave_low)),label="ExaWind")
hline!(ones(length(PSDInt_My_EXP_low))*(abs.(PSDInt_My_EXP_low)*plus_range_low_frequency),color=:red,linestyle=:dash,label="+25% from Experiment")
hline!(ones(length(PSDInt_My_EXP_low))*(abs.(PSDInt_My_EXP_low)*minus_range_low_frequency),color=:black,linestyle=:dash,label="-25% from Experiment")
ylabel!("Low-Frequency PSD Integral of \$M_y [(Nm)^2]\$")
savefig(bar_My,"bar_My_PSD_low_frequency.pdf")

println("My_low")
println(PSDInt_My_EXP_low)
println(PSDInt_My_CFD_low)
println(PSDInt_My_FAST_EXPWAVE_low)
println(PSDInt_My_FAST_CFDWAVE_low)
println(PSDInt_My_NaluWave_low)

bar_My_wave = bar((1,abs.(PSDInt_My_EXP_wave)),label="Experiment",legend=:outertopright)
bar!((2,abs.(PSDInt_My_CFD_wave)),label="StarCCM+")
bar!((3,abs.(PSDInt_My_FAST_EXPWAVE_wave)),label="OpenFAST (Experiment Input Waves)")
bar!((4,abs.(PSDInt_My_FAST_CFDWAVE_wave)),label="OpenFAST (StarCCM+ Input Waves)")
bar!((5,abs.(PSDInt_My_NaluWave_wave)),label="ExaWind")
hline!(ones(length(PSDInt_My_EXP_wave))*(abs.(PSDInt_My_EXP_wave)*plus_range_wave_frequency),color=:red,linestyle=:dash,label="+15% from Experiment")
hline!(ones(length(PSDInt_My_EXP_wave))*(abs.(PSDInt_My_EXP_wave)*minus_range_wave_frequency),color=:black,linestyle=:dash,label="-15% from Experiment")
ylabel!("Wave-frequency PSD Integral of \$M_y [(Nm)^2]\$")
savefig(bar_My_wave,"bar_My_PSD_wave_frequency.pdf")

println("My_wave")
println(PSDInt_My_EXP_wave)
println(PSDInt_My_CFD_wave)
println(PSDInt_My_FAST_EXPWAVE_wave)
println(PSDInt_My_FAST_CFDWAVE_wave)
println(PSDInt_My_NaluWave_wave)