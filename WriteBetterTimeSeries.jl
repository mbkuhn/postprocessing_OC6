using DelimitedFiles
using Interpolations
using DataFrames
using LaTeXStrings
using StatsBase
using Measurements
using Printf
include("get_evenly_spaced_time_series.jl")

# Model scale CG: to calculate displacements
x_CG = 0.
y_CG = 0.
z_CG = -1.506e-1

dt = 0.002

file_amrwind = "/Users/mkuhn/testruns_data/OC5_semisubmersible/HOS-NWT/amr-wind-standalone/same_res/output_same_res.txt"
file_exawind_pfx = "/Users/mkuhn/testruns_data/OC5_semisubmersible/HOS-NWT/floating_fixed_dt/platform_disp_etc"
file_exawind_sfx = ".dat"
file_exawind_labels = ["0","1","2","3","4","5","6"]

# NaluWind DATA

# Conversion between full scale and model scale
# L_p = 50 L_m
# T_p = sqrt(50)T_m
# F_p = L^3 F_m
# M_y = L^4 M_y
# theta_p = theta_m
rescaling = 50. 

file_exawind = file_exawind_pfx * file_exawind_labels[1] * file_exawind_sfx
Nalu_data = readdlm(file_exawind)
for label in file_exawind_labels
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

# convert from quaternions to angle
q0 = NaluWAVE_data[:,5]
q1 = NaluWAVE_data[:,6]
q2 = NaluWAVE_data[:,7]
q3 = NaluWAVE_data[:,8]
theta_y_Nalu_rad = -pi/2.0.+ 2.0.*atan.(sqrt.(1.0.+2.0.*(q0.*q2.-q1.*q3)),sqrt.(1.0.-2.0.*(q0.*q2.-q1.*q3)))
# convert to degrees
theta_y_Nalu = theta_y_Nalu_rad.*(180.0./pi)

# Convert time to full scale
NaluWave_time = NaluWAVE_data[:,1].*sqrt(rescaling)

# Convert x, z, and theta y to full scale after skipping transient portion of array
# Also subtract center of gravity to get displacements, not position
NaluWave_x = (NaluWAVE_data[:,2] .- x_CG).*(rescaling)
NaluWave_z = (NaluWAVE_data[:,4] .- z_CG).*(rescaling)
NaluWave_thy = theta_y_Nalu

# Get wave elevation from standalone simulation
amrwind_data = readdlm(file_amrwind)
amrwind_data = Float64.(amrwind_data)
amrwind_time = amrwind_data[:,1]
amrwind_wave = amrwind_data[:,2]
AMRWind_WAVE_time_interp = Interpolations.deduplicate_knots!(amrwind_time)
itp_waves_amrwind =  interpolate((AMRWind_WAVE_time_interp,), amrwind_wave, Gridded(Linear()))
amr_time_evenly_spaced = collect(amrwind_data[1,1]:dt:amrwind_data[end,1])
amrwind_wave_evenly_spaced = itp_waves_amrwind.(amr_time_evenly_spaced)
AMRWind_time = amr_time_evenly_spaced.*sqrt(rescaling)
AMRWind_Wave = amrwind_wave_evenly_spaced.*rescaling

# Create unified data object
nt = size(NaluWave_time)[1]
data = Array{Float64}(undef, nt, 5)
for n in 1:nt
    data[n,1] = NaluWave_time[n]
    data[n,2] = NaluWave_x[n]
    data[n,3] = NaluWave_z[n]
    data[n,4] = NaluWave_thy[n]
    data[n,5] = AMRWind_Wave[n]
end

# Write to file together
header = ["time(s)" "surge(m)" "heave(m)" "pitch(deg)" "elevation(m)"]
writedlm("raw_displacements_elevation.txt",[header; data] , '\t')

# Define shift, offset
shift = 50
offset = 89200 + 5000
t_offset = offset * dt * sqrt(rescaling)

# Redefine data
nt_off = nt - offset
data = Array{Float64}(undef, nt_off, 5)
for n in 1:nt_off
    data[n,1] = NaluWave_time[n + offset] - t_offset
    data[n,2] = NaluWave_x[n + offset]
    data[n,3] = NaluWave_z[n + offset]
    data[n,4] = NaluWave_thy[n + offset]
    data[n,5] = AMRWind_Wave[n + offset - shift]
end

# Write to file together
writedlm("adjusted_displacements_elevation.txt",[header; data] , '\t')