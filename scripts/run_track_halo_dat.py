import numpy as np
import track_halo_number
#outputs=['003840', '003840', '003968', '003840', '003968', '003840', '003840', '003840', '003840', '003840' ]
#models=['GM2', 'GM2noBHs', "GM2SI1", "GM3", "GM3noBHs", "GM3SI1", "P0", "P0noBHs","GM1", "GM1noBHs"]
outputs=["003840"]
models=["GM3"]
halos = np.ones(len(models)).astype(int)
track_halo_number.halo_num_tracker(models, outputs, halos, 'GM3_track_halos.dat')
