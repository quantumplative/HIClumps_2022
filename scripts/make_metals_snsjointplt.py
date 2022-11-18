import test_timetrace_clumpgrps

models = ['P0noBHs', 'GM1', 'GM1noBHs', 'GM2', 'GM2noBHs', 'GM3SI1']
output_index = [35, 55, 17, 35, 34, 22] 
test_timetrace_clumpgrps.make_metal_profile_jointplts(models, output_index)
