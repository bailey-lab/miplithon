import pickle
UMI_counts=snakemake.input.UMI_counts
UMI_counts=pickle.load(open(UMI_counts, 'rb'))
metadata=snakemake.input.metadata
sample_column=snakemake.params.sample_column
column_names=open(metadata).readline().strip().split(',')
header_dict={}
for column_number, column_name in enumerate(column_names):
	header_dict[column_name]=column_number
sample_column=header_dict[sample_column]
valid_samples=set([line.strip().split(',')[sample_column] for line in open(metadata)][1:])



N_ref=snakemake.output.N_ref
D_ref=snakemake.output.D_ref
N_cov=snakemake.output.N_cov
D_cov=snakemake.output.D_cov

N_mutation=snakemake.params.N_mutation
D_mutation=snakemake.params.D_mutation


refs=[N_ref, D_ref]
covs=[N_cov, D_cov]
mutations=[N_mutation, D_mutation]

NFD_cov_file=open(snakemake.output.NFD_cov, 'w')
NFD_alt_file=open(snakemake.output.NFD_alt, 'w')
F_alt=set([line.strip() for line in open(snakemake.input.F_alt)])
F_cov=set([line.strip() for line in open(snakemake.input.F_cov)])



def check_thresholds(valid_samples, mutation):
	cov_samples, ref_samples=set([]), set([])
	for sample in UMI_counts:
		if sample in valid_samples:
			for UMI_mutation in UMI_counts[sample]:
				if mutation in UMI_mutation:
					cov, ref, alt=UMI_counts[sample][UMI_mutation]
					if cov>=cov_threshold:
						cov_samples.add(sample)
					if ref>=ref_threshold and cov>=cov_threshold: #unlike most of these, in this case we're interested in a non-zero ref count
						ref_samples.add(sample)
	return cov_samples, ref_samples

NFD_cov, NFD_alt=F_cov, F_alt
for mutation_number, mutation in enumerate(mutations):
	print('working on', mutation)
	mutation,cov,ref=mutation.split('/')[-1].split('_')[:3]
	cov_threshold, ref_threshold=int(cov), int(ref)
	cov_samples, ref_samples=check_thresholds(valid_samples, mutation)
	ref_file=open(refs[mutation_number], 'w')
	cov_file=open(covs[mutation_number], 'w')
	ref_file.write('\n'.join(ref_samples))
	cov_file.write('\n'.join(cov_samples))
	NFD_cov=NFD_cov & cov_samples
	NFD_alt=NFD_alt & ref_samples
NFD_cov_file.write('\n'.join(NFD_cov))
NFD_alt_file.write('\n'.join(NFD_alt))
