import tangos as db 
import numpy as np 
import sys 

def new_halo_num_tracker(model, sim_name, ini_stepnum, ini_halo_num):
	halo = db.get_halo(sim_name + '/%' + str(ini_stepnum) + '/halo_' + str(ini_halo_num))
	sim = db.get_simulation(sim_name) 

	#halo_cat = [] 
	#step_num = [] 
	outfn   = '../datfiles/' +  model + '_halo_track_info.dat'	
	outfile = open(outfn, 'w')
	outfile.write('model output z time halo\n')
	
	for i in range(len(halo.calculate_for_progenitors('z()')[0])-1):
		halo_num = halo.calculate('earlier('+str(i+1)+')').halo_number
		step_num = halo.calculate('earlier('+str(i+1)+')').timestep
		z = halo.calculate('earlier('+str(i+1)+')').calculate('z()')
		t = halo.calculate('earlier('+str(i+1)+')').calculate('t()')
		outfile.write('%s %i %.2f %.2f %i\n' % (model, int(step_num.extension[-6:]), z, t, halo_num))
		#print('sim_name length = ' + str(len(sim_name)))
		#print('tangos time step out length = ' + str(len(str(halo.calculate('earlier('+str(i+1)+')').timestep))))
		#np.savetxt(output_base_f + '_halo_cat.txt', halo_cat) 
		#np.savetxt(output_base_f + 'timesteps.txt', step_num)
	outfile.close()
	#return halo_cat, step_num 
