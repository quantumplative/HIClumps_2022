import pyfof_HI_clumps
model    = 'h148_6144' 
output   = '003456'
link_len = 1 
HI_cut   = 4e-08
min_m    = 2
pyfof_HI_clumps.pyfof_clump_data(model, output, link_len, HI_cut, min_m)
