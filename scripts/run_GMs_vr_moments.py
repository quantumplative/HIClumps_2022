import test_timetrace_clumpgrps

#models   = ['P0', 'GM1noBHs', 'GM2', 'GM2noBHs']
#si_array = [6, 0, 5, 3]
#ei_array = [19, 19,  37, 35]
#models = ['GM1', 'GM3SI1']
#si_array = [5, 0]
#ei_array = [56, 27]
#models = ['GM1']
#si_array = [5]
#ei_array = [56]
models = ['GM3SI1']
si_array = [0]
ei_array = [27]
test_timetrace_clumpgrps.make_radcomvel_moments(models, si_array, ei_array)
