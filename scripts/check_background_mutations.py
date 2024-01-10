'''
checks a series of background mutations to see if they're present in R561H mutants
'''
background_mutations=snakemake.input.background_mutations
foreground_mutation=snakemake.input.foreground_mutation
foreground_samples=set([line.strip() for line in open(foreground_mutation)])
overlapping_mutations=open(snakemake.output.background_mutations, 'w')

overlapping_mutations.write('mutation_name_+_thresholds\t#_of_bg_samples\t#_of_bg_samples_in_561H\tshared_sample_names\n')
for background_mutation in background_mutations:
	background_samples=set([line.strip() for line in open(background_mutation)])
	overlapping_mutations.write(f'{background_mutation}\t{len(background_samples)}\t{len(background_samples & foreground_samples)}\t{background_samples & foreground_samples}\n')
	
