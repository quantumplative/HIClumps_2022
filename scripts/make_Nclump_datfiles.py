import test_timetrace_clumpgrps

models   = ['P0', 'GM1noBHs', 'GM2', 'GM2noBHs', 'GM1', 'GM3SI1', 'P0noBHs']
si_array = [6, 0, 5, 3, 5, 0, 4]
ei_array = [19, 19,  37, 35, 56, 27, 38]
for i in range(len(models)):
	try:
		test_timetrace_clumpgrps.make_clump_count_datfile(models[i], si_array[i], ei_array[i])
	except: 
		print('check model ', models[i])
