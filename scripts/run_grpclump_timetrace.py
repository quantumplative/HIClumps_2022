import pyfof_HI_clumps
models=['GM2', 'GM2noBHs', "GM3", "GM3noBHs", "GM3SI1", "P0", "P0noBHs","GM1", "GM1noBHs"]
for m in models:
	pyfof_HI_clumps.timetrace_grpclumps(m, 1, 4e-08, 2)
	print("done with model " + m)
