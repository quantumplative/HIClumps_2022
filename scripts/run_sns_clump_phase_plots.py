import numpy as np
import sns_clump_phase_plots

dict_outputs={"GM2": '003840', "GM2noBHs": '003840', "GM2SI1": '003968', "GM3": '003840', "GM3noBHs": '003968', "GM3SI1": '003840', "P0": '003840', "P0noBHs": '003840', "GM1": '003840', "GM1noBHs": '003840'}
#all_models = ['P0', 'P0noBHs', 'GM1', 'GM1noBHs']
all_models=['GM2', 'GM2noBHs', 'GM2SI1', 'GM3', 'GM3noBHs', 'GM3SI1']
#all_models = ['P0']
outputs = [dict_outputs['GM2'],dict_outputs['GM2noBHs'], dict_outputs['GM2SI1'], dict_outputs['GM3'], dict_outputs['GM3noBHs'], dict_outputs['GM3SI1']]

all_models=np.flip(all_models)
outputs = np.flip(outputs)

#sns_clump_phase_plots.scatter_clump_phase(all_models, '003456')
outputs = ['001536', '003456']
for m in outputs:
	#sns_clump_phase_plots.joint_clump_phase(all_models, m, 'vdisp_temp_phase_DMQ_output_' + m + '.pdf')
	sns_clump_phase_plots.joint_clump_phase(all_models, m, 'vdisp_pressure_phase_DMQ_output_' + m + '.pdf', 'vdisp', 'pressure', True, False)
