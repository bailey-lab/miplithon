'''
unlike get_prevalences, which retrieves the prevalences of specific mutations,
this program retrieves all mutations (either within a specific gene or all
genes) that meet a threshold coverage and a threshold alternate count. This
version imposes an additional filter by requiring at least x samples to contain
the mutation of interest in order to report it.
'''
import pickle
filtered_mutations=snakemake.output.filtered_mutations
sample_column=snakemake.params.sample_column
search_term, cov_threshold, alt_threshold, count_threshold=filtered_mutations.split('/')[-1].split('_')[:4]
cov_threshold, alt_threshold, count_threshold=map(int, [cov_threshold, alt_threshold, count_threshold])
filtered_mutations=open(filtered_mutations, 'w')
metadata=snakemake.input.metadata
header_dict={}
column_names=open(metadata).readline().strip().split(',')
for column_number, column_name in enumerate(column_names):
	header_dict[column_name]=column_number
sample_column=header_dict[sample_column]
valid_samples=set([line.strip().split(',')[sample_column] for line in open(metadata)][1:])
#print('valid samples are', valid_samples)

amino_acid_set=set(['Ala', 'Arg', 'Asn', 'Asp', 'Cys', 'Glu', 'Gln', 'Gly',
'His', 'Ile', 'Leu', 'Lys', 'Met', 'Phe', 'Pro', 'Ser', 'Thr', 'Trp', 'Tyr',
'Val'])

UMI_counts=snakemake.input.UMI_counts
UMI_counts=pickle.load(open(UMI_counts, 'rb'))

def check_thresholds(valid_samples, search_term, cov_threshold, alt_threshold, count_threshold):
	cov_dict={}
	if search_term=='all':
		search_term=''
	for sample in UMI_counts:
		if sample in valid_samples:
			for mutation in UMI_counts[sample]:
				if search_term in mutation:
					parsed_mutation=mutation.split('-')[-1]
					start, end=parsed_mutation[:3], parsed_mutation[-3:]
					middle=parsed_mutation[3:-3]
					if middle.isdigit() and start in amino_acid_set and end in amino_acid_set and start!=end:
						cov_dict.setdefault(mutation, [set([]), set([])])
						cov, ref, alt=UMI_counts[sample][mutation]
						if cov>=cov_threshold:
							cov_dict[mutation][0].add(sample)
						if alt>=alt_threshold and cov>=cov_threshold:
							cov_dict[mutation][1].add(sample)
	for mutation in list(cov_dict.keys()):
		if len(cov_dict[mutation][1])<count_threshold:
			cov_dict.pop(mutation)
	return cov_dict



cov_dict=check_thresholds(valid_samples, search_term, cov_threshold, alt_threshold, count_threshold)
filtered_mutations.write('\n'.join(list(cov_dict.keys())))
cov_threshold, alt_threshold, count_threshold=map(str, [cov_threshold, alt_threshold, count_threshold])
for mutation in cov_dict:
	mutation_str=mutation.replace(' ', '_')
	cov_file=open('prevalences_by_threshold/'+'_'.join([mutation_str, cov_threshold, alt_threshold, count_threshold, 'cov.txt']), 'w')
	alt_file=open('prevalences_by_threshold/'+'_'.join([mutation_str, cov_threshold, alt_threshold, count_threshold, 'alt.txt']), 'w')
	cov_file.write('\n'.join(list(cov_dict[mutation][0])))
	alt_file.write('\n'.join(list(cov_dict[mutation][1])))
