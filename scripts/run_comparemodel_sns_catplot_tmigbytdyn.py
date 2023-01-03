import timescales
import matplotlib.pyplot as plt 
models =  ['P0', 'GM1', 'GM1noBHs', 'GM2']
outputs = [384, 512, 972, 1152, 1536, 2554, 3456]
timescales.comparemodel_sns_catplot_tmigbytdyn(models, outputs)
