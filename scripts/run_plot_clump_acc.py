import timescales 

#models=['P0']
models =  ['P0', 'GM1', 'GM1noBHs', 'GM2']
timescales.plot_mean_clump_acc(models)
timescales.plot_median_clump_acc(models)
