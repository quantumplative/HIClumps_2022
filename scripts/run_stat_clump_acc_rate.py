import timescales

models = ['P0', 'P0noBHs', 'GM1', 'GM1noBHs', 'GM2']
for m in models:
	timescales.stat_clump_acc_rate(m)
