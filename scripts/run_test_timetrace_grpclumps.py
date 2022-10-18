import test_timetrace_clumpgrps
import pyfof_HI_clumps
m = 'GM3noBHs'
halo, output = pyfof_HI_clumps.get_model_outputs_halos(m)
si = 28
for l in range(len(output)-si): 
	try:
		#test_timetrace_clumpgrps.timetrace_grpclumps_metals_vr(m, 1, 4e-08, 2, l+si)
		test_timetrace_clumpgrps.timetrace_grpclumps(m, 1, 4e-08, 2, l+si)
	except: 
		print('check on index l = ', l)
