import pickle
UMI_counts=snakemake.input.UMI_counts
UMI_counts=pickle.load(open(UMI_counts, 'rb'))
metadata=snakemake.input.metadata
coverage_files=snakemake.output.coverage_files
alternate_files=snakemake.output.alternate_files
sample_column=snakemake.params.sample_column

def check_thresholds(sample):
	num, denom=False, False
	found_mutation=False
	if sample in UMI_counts:
		for UMI_mutation in UMI_counts[sample]:
			if mutation in UMI_mutation:
				found_mutation=True
				cov, ref, alt=UMI_counts[sample][UMI_mutation]
#				if cov!=alt+ref:
#					print(sample, mutation, cov, ref, alt)
				if cov>=cov_threshold:
					denom=True
				if alt>=alt_threshold and cov>=cov_threshold:
					num=True
	return num, denom, found_mutation

h={}
for sample_number, coverage_sample in enumerate(coverage_files):
	print('working on', coverage_sample)
	cov_file=open(coverage_sample, 'w')
	alt_file=open(alternate_files[sample_number], 'w')
	mutation,region,cov,alt=coverage_sample.split('/')[-1].split('_')[:-2]
	cov_threshold, alt_threshold=int(cov), int(alt)
	region_type, region_value=region.split(':')
	overall_found=False
	for line_number, line in enumerate(open(metadata)):
		line=line.strip().split(',')
		if line_number==0:
			for column_number, column in enumerate(line):
				h[column]=column_number
		else:
			sample=line[h[sample_column]]
			if region_value=='all' or region_value==line[h[region_type]]:
				num, denom, found=check_thresholds(sample)
				if found:
					overall_found=True
				if denom:
					cov_file.write(sample+'\n')
				if num:
					alt_file.write(sample+'\n')
	if not overall_found:
		print(mutation, 'not a found mutation')
