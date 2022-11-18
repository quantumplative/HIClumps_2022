import test_timetrace_clumpgrps

models = ['GM1', 'GM1noBHs', 'GM2', 'GM2noBHs', 'GM3SI1']
output_index = [55, 17, 35, 34, 22] 

test_timetrace_clumpgrps.GMs_metallicity_plot(models, output_index, True, False, False)	
