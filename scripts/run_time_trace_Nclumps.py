import sf_rates

outfn = 'GMs_track_halos_v2.dat'
models=['GM1']
#'GM2noBHs', "GM2SI1", "GM3", "GM3noBHs", "GM3SI1", "P0", "P0noBHs","GM1", "GM1noBHs"]
#models=['GM2noBHs', "GM2SI1", "GM3", "GM3noBHs", "GM3SI1", "P0", "P0noBHs","GM1", "GM1noBHs"]
sf_rates.time_trace_Nclumps(models,  outfn)
