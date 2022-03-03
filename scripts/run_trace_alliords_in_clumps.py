import pyfof_HI_clumps
#models=['GM2', 'GM2noBHs', "GM3", "GM3noBHs", "GM3SI1", "P0", "P0noBHs","GM1", "GM1noBHs"]
models=['P0']
for m in models:
	print('working on model ' + m)
	pyfof_HI_clumps.trace_alliords_ever_inclumps(m) 
	print("done with model " + m)
