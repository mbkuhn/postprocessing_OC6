# postprocessing_OC6
Post-processing tools for the computation of PSD

To run the codes, first ensure you have your data saved in convenient locations.
Then,
From terminal, type: ``julia WaveLoadsOnFixedStructure.jl"

From the Julia REPL (i.e. after having typed ``julia" in your terminal), type: include("WaveLoadsOnFixedStructure.jl")

VSCODE-specific: ctrl/cmd + A the contents of WaveLoadsOnFixedStructure.jl and then press shift+enter (good for debugging purposes)

### Abridged motion comparison

Run `julia MotionOfFloatingStructureAbridged.jl` after configuring its input paths.
The script compares a current ExaWind run against a new run, currently represented
by the experimental text file as an explicitly labelled placeholder. It requires
that file and the current ExaWind restart files, not the StarCCM+ MAT file.
It creates `outputs_abridged/` if needed and saves:

- `displacement_time_series.pdf`: surge, heave, and pitch comparisons.
- `x_PSD.pdf`, `z_PSD.pdf`, `thy_PSD.pdf`: smoothed displacement spectra.

All comparisons are in model scale. Before use, the experimental placeholder's
time is divided by `sqrt(rescaling)` and surge/heave by `rescaling`; pitch remains
in degrees. The current run is already at model scale, with its center of gravity
subtracted to obtain displacements. To use actual new-run data, replace the
placeholder-loading block with model-scale `new_run_time`, `new_run_x`,
`new_run_z`, and `new_run_thy` arrays (pitch in degrees).

Plots and spectra retain the intersection of both runs' available times and the
configured interval: `t_skip` defaults to 37 model-scale seconds, while
`t_final = Inf` uses all shared remaining data, regardless of which run ends
first. Set a finite `t_final` to impose an additional end limit. Frequency bands
and plot limits are also converted to model scale.

The terminal reports bias (new minus current), mean
absolute error (MAE), root mean square error (RMSE), and maximum absolute error,
in model-scale metres for surge/heave and degrees for pitch. The current run is
linearly interpolated onto the retained new-run timestamps without extrapolation,
phase alignment, or mean removal. Both spectra use their native sampling grids
within the shared interval. No bar graphs are generated.
