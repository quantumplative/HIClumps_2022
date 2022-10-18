import test_timetrace_clumpgrps
#models   = ['P0', 'GM1noBHs', 'GM2', 'GM2noBHs', 'GM1', 
models = ['GM3SI1', 'P0noBHs']
#si_array = [6, 0, 5, 3, 5, 
si_array = [0, 4]
#ei_array = [19, 19,  37, 35, 56, 
ei_array = [27, 38]
for i in range(len(models)):
	test_timetrace_clumpgrps.make_migration_free_fall_time(models[i], si_array[i], ei_array[i]) 
