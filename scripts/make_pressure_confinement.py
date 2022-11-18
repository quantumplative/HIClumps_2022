import test_timetrace_clumpgrps

#models = ['P0', 'P0noBHs', 'GM1', 'GM1noBHs', 'GM2', 'GM2noBHs'
#models = ['GM3noBHs', 'GM3SI1']
models = ['GM2']
#output 3456 z = 0.17 
#output_index = [18, 35, 55, 17, 35, 34] 
#output_index = [35, 22]
output_index = [6] 

#test_timetrace_clumpgrps.allGMs_pressure_confinement_plot(models, output_index, True, False)	
#messing up bc of GM1 step 3456
for i in range(len(models)):
	#try: 
	test_timetrace_clumpgrps.pressure_confinement_plot(models[i], output_index[i], True, False)	
	#	print('done with model ' + models[i])
	#except:
	#	print('something up with model ' + models[i])
