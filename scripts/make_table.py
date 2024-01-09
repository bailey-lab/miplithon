prevalence_files=snakemake.input.desired_files
metadata_sheet=snakemake.input.metadata_sheet
summary=snakemake.output.summary
heirarchy=snakemake.params.heirarchy
heirarchy=heirarchy.split(':')
sample_column=snakemake.params.sample_column

def make_metadata_dict(metadata_sheet):
	metadata_dict={}
	h={}
	for line_number, line in enumerate(open(metadata_sheet)):
		line=line.strip().split(',')
		if line_number==0:
			for column_number, column in enumerate(line):
				h[column]=column_number
		else:
			sample=line[h[sample_column]]
			metadata_dict[sample]=line
	return metadata_dict, h

metadata_dict, header_dict=make_metadata_dict(metadata_sheet)

#prevalences/Region:all_3_1_summary.tsv

def get_counts(sample_file, category, filter_type):
	parsed_counts={}
	samples=[line.strip() for line in open(sample_file)]
	category_index=heirarchy.index(category)
	if len(heirarchy)>1:
		subsidiary=heirarchy[category_index+1]
	else:
		filter_type='all'
	for sample in samples:
		metadata=metadata_dict[sample]
		category_value=metadata[header_dict[category]]
		if filter_type=='all':
			parsed_counts[category_value]=parsed_counts.setdefault(category_value, 0)+1
		else:
			subsidiary_value=metadata[header_dict[subsidiary]]
			parsed_counts[subsidiary_value]=parsed_counts.setdefault(subsidiary_value, 0)+1
		parsed_counts['overall']=parsed_counts.setdefault('overall', 0)+1
	return parsed_counts

def format_line(prevalence_dict, category_value, filtered_files, output_table):
		output_line=[category_value]
		for prevalence_file in filtered_files:
			cov_dict, alt_dict=prevalence_dict[prevalence_file]
			cov_count=cov_dict.setdefault(category_value, 0)
			alt_count=alt_dict.setdefault(category_value, 0)
			if cov_count>0:
				output_line.append(f'{round(alt_count/cov_count, 4)} ({alt_count}/{cov_count})')
			else:
				output_line.append(f'0.0 ({alt_count}/{cov_count})')
		output_table.write('\t'.join(output_line)+'\n')

def format_table(prevalence_dict, category_values, output_header, filtered_files):
	output_table=open(summary, 'w')
	output_table.write('\t'.join(output_header)+'\n')
	for category_value in category_values:
		format_line(prevalence_dict, category_value, filtered_files, output_table)
	format_line(prevalence_dict, 'overall', filtered_files, output_table)
#k13-Arg561His_Region:all_3_1_cov_samples.txt

prevalence_dict={}
category_type, filter_type=summary.split('/')[-1].split('_')[0].split(':')
categories=set([])
output_header=[category_type]
filtered_files=[]
for prevalence_file in prevalence_files:
	prefix='_'.join(prevalence_file.split('_')[:4])
	mutation, region, cov, alt=prevalence_file.split('_')[:4]
	obs_category, obs_filter=region.split(':')
	if obs_filter==filter_type and obs_category==category_type:
		filtered_files.append(prevalence_file)
		output_header.append(mutation.split('/')[-1])
		cov_counts=get_counts(prevalence_file, category_type, filter_type)
		alt_counts=get_counts(prefix+'_alt_samples.txt', category_type, filter_type)
		prevalence_dict[prevalence_file]=[cov_counts, alt_counts]
		categories=categories | set(cov_counts.keys()) | set(alt_counts.keys())
categories=categories-set(['overall'])
format_table(prevalence_dict, categories, output_header, filtered_files)
