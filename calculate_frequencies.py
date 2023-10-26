'''
3193 (expected position based on crt lookup in final table) vs 3132 (actual position in output list)
'''
reference='/nfs/jbailey5/baileyweb/bailey_share/processed_data/afola/jason_paris_mip_project/mergeddata/variant/reference_AA_table.csv'
alternate='/nfs/jbailey5/baileyweb/bailey_share/processed_data/afola/jason_paris_mip_project/mergeddata/variant/alternate_AA_table.csv'

output_file=open('alternate_freqs.csv', 'w')

def get_counts(input_file):
	count_dict={}
	header=[]
	for line_number, line in enumerate(open(input_file)):
		line=line.strip().split(',')
		if line_number==2:
			mutations=line[1:]
		if line_number>5:
			counts=line[1:]
			sample=line[0]
			for column_number, count in enumerate(counts):
				if sample not in count_dict:
					count_dict[sample]={}
				if column_number not in count_dict[sample]:
					count_dict[sample][column_number]=float(count)
		else:
			header.append(','.join(line)+'\n')
	return header, count_dict, mutations

def plot_heatmap(graphing_list, x_values, y_values, x_title, y_title, count_title, output_path, width=2000, height=4000):
	import plotly.express as px
	fig = px.imshow(graphing_list, aspect='auto', labels=dict(x=x_title, y=y_title,
	color=count_title), x=x_values, y=y_values)
	fig.update_xaxes(side="top")
	fig.update_layout(width=width, height=height, autosize=False)
	fig.write_html(output_path)

header, ref_counts, mutations=get_counts(reference)
header, alt_counts, mutations=get_counts(alternate)

for line in header:
	output_file.write(line)

samples=[]
graphing_list=[]
printing=False
for sample_number, sample in enumerate(ref_counts):
	samples.append(sample)
	output_line=[sample]
	graphing_line=[]
	for column_number in ref_counts[sample]:
		alt_count=alt_counts[sample][column_number]
		ref_count=ref_counts[sample][column_number]
		if alt_count+ref_count>0:
			fraction=alt_count/(alt_count+ref_count)
		else:
			fraction=-1
		graphing_line.append(fraction)
		output_line.append(str(fraction))
	graphing_list.append(graphing_line)
	output_file.write(','.join(output_line)+'\n')

plot_heatmap(graphing_list, mutations, samples, 'mutations', 'samples', 'mutation frequency', 'frequency_heatmap.html')
