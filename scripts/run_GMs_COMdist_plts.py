import test_timetrace_clumpgrps
models   = ['P0', 'P0noBHs', 'GM2', 'GM1noBHs', 'GM2noBHs']
si_array = [6, 4, 5, 0, 3]
ei_array = [19, 38, 37, 19, 35] 
cmaps    = ['mako', 'crest', 'rocket_r', 'Blues', 'flare'] 
show      = False
for i in range(len(models)): 
	test_timetrace_clumpgrps.COMdist_kde(models[i],  si_array[i], ei_array[i], cmaps[i], show)
	#test_timetrace_clumpgrps.COMdist_hist(models[i],  si_array[i], ei_array[i], cmaps[i], show)
