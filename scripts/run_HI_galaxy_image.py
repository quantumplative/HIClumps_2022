import pyfof_HI_clumps
model = 'P0'
HIcut = 4e-08
output = '003456'
width =600
res = 1000
vmin = 1e-11
vmax = 5e-2
show = True 
pyfof_HI_clumps.HI_galaxy_image(model, output, HIcut, width, res, vmin, vmax, show)
