# postprocessing_OC6
Post-processing tools for the computation of PSD

To run the codes, first ensure you have your data saved in convenient locations.
Then,
From terminal, type: ``julia WaveLoadsOnFixedStructure.jl"

From the Julia REPL (i.e. after having typed ``julia" in your terminal), type: include("WaveLoadsOnFixedStructure.jl")

VSCODE-specific: ctrl/cmd + A the contents of WaveLoadsOnFixedStructure.jl and then press shift+enter (good for debugging purposes)

### Abridged motion comparison

Run `julia MotionOfFloatingStructureAbridged.jl` after configuring its input paths.
This script requires only two sets of files, not
the StarCCM+ MAT file. It creates `outputs_abridged/` if needed and saves:

- `displacement_time_series.pdf`: surge, heave, and pitch comparisons.
- `x_PSD.pdf`, `z_PSD.pdf`, `thy_PSD.pdf`: smoothed displacement spectra.

Plots and spectra use the overlapping interval after the configured 
transient cutoff (`t_skip` and `t_final`, in model-scale seconds). All comparisons
are at full scale. The terminal reports bias, mean
absolute error (MAE), root mean square error (RMSE), and maximum absolute error,
in metres for surge/heave and degrees for pitch. 
