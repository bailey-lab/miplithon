reference='/home/dpereus/big_data/tz_miseq_run/230915_variant/analysis/reference_AA_table.csv'
alternate='/home/dpereus/big_data/tz_miseq_run/230915_variant/analysis/alternate_AA_table.csv'

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
			for count_number, count in enumerate(counts):
				mutation=mutations[count_number]
				if sample not in count_dict:
					count_dict[sample]={}
				if mutation not in count_dict[sample]:
					count_dict[sample][mutation]=float(count)
		else:
			header.append(','.join(line)+'\n')
	return header, count_dict

header, ref_counts=get_counts(reference)
header, alt_counts=get_counts(alternate)

for line in header:
	output_file.write(line)

samples, mutations=[],[]
graphing_list=[]
for sample_number, sample in enumerate(ref_counts):
	samples.append(sample)
	output_line=[sample]
	graphing_line=[]
	for mutation in ref_counts[sample]:
		if sample_number==0:
			mutations.append(mutation)
		ref_count=ref_counts[sample][mutation]
		alt_count=alt_counts[sample][mutation]
		if alt_count+ref_count>0:
			fraction=alt_count/(alt_count+ref_count)
		else:
			fraction=-1
		graphing_line.append(fraction)
		output_line.append(str(fraction))
	graphing_list.append(graphing_line)
	output_file.write(','.join(output_line)+'\n')



def plot_heatmap(graphing_list, x_values, y_values, x_title, y_title, count_title, output_path, width=2000, height=4000):
	import plotly.express as px
#	print(graphing_list)
	fig = px.imshow(graphing_list, aspect='auto', labels=dict(x=x_title, y=y_title,
	color=count_title), x=x_values, y=y_values)
	fig.update_xaxes(side="top")
	fig.update_layout(width=width, height=height, autosize=False)
	fig.write_html(output_path)

plot_heatmap(graphing_list, mutations, samples, 'mutations', 'samples', 'mutation frequency', 'frequency_heatmap.html')
