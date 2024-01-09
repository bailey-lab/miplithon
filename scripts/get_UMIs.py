import pickle
cov_counts=snakemake.input.cov_counts
alt_counts=snakemake.input.alt_counts
ref_counts=snakemake.input.ref_counts
UMI_counts=open(snakemake.output.UMI_counts, 'wb')
UMI_suffix=snakemake.params.UMI_suffix

def get_AA_counts(input_file):
	mutation_dict={}
	for line_number, line in enumerate(open(input_file)):
		if line_number==2:
			mutations=line.strip().split(',')[1:]
		if line_number>5:
			line=line.strip().split(',')
			sample=line[0].replace(UMI_suffix, '')
			counts=line[1:]
			mutation_dict[sample]=list(map(int, (map(float, counts))))
	return mutations, mutation_dict

def combine_counts(cov_counts, ref_counts, alt_counts, AA_mutations):
	combined_dict={}
	for sample in cov_counts:
		combined_dict.setdefault(sample, {})
		for mutation_number, mutation in enumerate(cov_counts[sample]):
			cov_count=cov_counts[sample][mutation_number]
			ref_count=ref_counts[sample][mutation_number]
			alt_count=alt_counts[sample][mutation_number]
			mutation_name=str(mutation_number)+': '+AA_mutations[mutation_number]
			combined_dict[sample][mutation_name]=[cov_count, ref_count, alt_count]
	return combined_dict

AA_mutations, cov_counts=get_AA_counts(cov_counts)
junk, ref_counts=get_AA_counts(ref_counts)
junk, alt_counts=get_AA_counts(alt_counts)

#final=[cov_counts, ref_counts, alt_counts, AA_mutations]
final=combine_counts(cov_counts, ref_counts, alt_counts, AA_mutations)
pickle.dump(final, UMI_counts)
