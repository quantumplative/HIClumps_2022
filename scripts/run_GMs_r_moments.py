import test_timetrace_clumpgrps

models   = ['P0', 'GM1noBHs', 'GM2', 'GM2noBHs', 'GM1', 'GM3SI1', 'P0noBHs']
si_array = [6, 0, 5, 3, 5, 0, 4]
ei_array = [19, 19,  37, 35, 56, 27, 38]
#models   = ['GM2noBHs', 'GM1', 'GM3SI1', 'P0noBHs']
#si_array = [3, 5, 0, 4]
#ei_array = [35, 56, 27, 38]
test_timetrace_clumpgrps.make_comdist_moments(models, si_array, ei_array)
