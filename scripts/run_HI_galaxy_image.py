import pyfof_HI_clumps
model = 'h148_6144'
HIcut = 4e-08
output = '003456'
width =800
res = 1000
vmin = 1e-11
vmax = 5e-2
show = False 
pyfof_HI_clumps.HI_galaxy_image(model, output, HIcut, width, res, vmin, vmax, show)
