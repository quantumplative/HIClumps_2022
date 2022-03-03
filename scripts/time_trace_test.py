pyfof_HI_clumps.get_fn('GM2')
data = pd.read_table(fn +'timetrace_grps.dat', sep='\s+')
data
data[0]
data[0:]
data[0:1]
data[0:1]['output1']
data_output1 = data[data['output1'] == int(data[0:1]['output1'])]
data_output1
data_output1_tracegrps = data_output1[data_output1['grpmatch']>0]
data_output1_tracegrps
data_output1_tracegrps = data_output1[data_output1['grpmatch']>0 and data_output1['fraction'] >= .50]
data_output1_tracegrps
data_output1_tracegrps[data_output1_tracegrps['fraction'] >=0.5]
data_output1_tracegrps = data_output1_tracegrps[data_output1_tracegrps['fraction'] >=0.5]
np.unique(data_outpu1_tracegrps['grpmatch'])
np.unique(data_output1_tracegrps['grpmatch'])
uni_tracegrps = np.unique(data_output1_tracegrps['grpmatch'])
len(uni_tracegrps)
history
data[0:1]['output1']
data[1:2]['output1']
data['output1']
outputs = np.unique(data['output1'])
outputs
outputs = np.flip(np.unique(data['output1']))
outputs
data[0:1]['output1']
data_output2 = data[data['output1'] == outputs[1]]
data_output2
Data_output1mo
