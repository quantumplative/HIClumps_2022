import sns_clump_phase_plots
model_type = ['all']
models = [['P0', 'P0noBHs', 'GM1', 'GM1noBHs', 'GM2', 'GM2noBHs', 'GM2SI1', 'GM3', 'GM3noBHs', 'GM3SI1']]
outputs = ['001536', '003456']
legend=True
show=True
phasex='size'
phasey='mass'
sns_clump_phase_plots.grid_phase(model_type, models, outputs, phasex, phasey, legend, show)
