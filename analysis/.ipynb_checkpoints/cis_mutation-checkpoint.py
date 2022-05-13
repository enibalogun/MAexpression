#!/usr/bin/env python
import pandas as pd
import numpy as np
import gffpandas.gffpandas as gffpd

# +
#### CREATING DATAFRAMES AND DICTIONARIES ####
CC2344 = ["CC2344-L1", "CC2344-L2", "CC2344-L3", "CC2344-L4", "CC2344-L5", "CC2344-L6", "CC2344-L7", "CC2344-L8", "CC2344-L9", "CC2344-L10", "CC2344-L11", "CC2344-L12", "CC2344-L13", "CC2344-L14", "CC2344-L15"]
CC2931 = ["CC2931-L1", "CC2931-L2", "CC2931-L3", "CC2931-L4", "CC2931-L5", "CC2931-L6", "CC2931-L7", "CC2931-L9", "CC2931-L10", "CC2931-L11", "CC2931-L13", "CC2931-L14", "CC2931-L15"]
total = CC2344 + CC2931

DEGs_dic = {}
enriched_DEGs = pd.DataFrame()

### OPENING FILES ####
mutations = pd.read_csv('/research/projects/chlamydomonas/MAexpression/data/genome_info/mutation_info/all_mutations.csv', delimiter = '\t')
DEGs = pd.read_csv('/research/projects/chlamydomonas/MAexpression/analysis/DEGs/total_genes1.csv', delimiter = ',')
DEGs = DEGs.replace('intergenic', np.nan, regex = True)

for i in DEGs.columns:
    DEGs_dic[i] = DEGs[i].dropna().values.tolist()
    
#### IDENTIFYING MUTATIONS ENRICHED IN OBSERVED DEGS ####
for i in CC2344 + CC2931:
    gene_list = DEGs_dic[i]
    spec_mut = mutations.loc[mutations['sample'] == i].groupby('gene').count().reset_index()
    unique_mut = spec_mut['gene'].values.tolist()
    intersect = list(set(gene_list).intersection(set(unique_mut)))
    a = spec_mut.loc[spec_mut['gene'].isin(intersect)]
    enriched_DEGs.at['observed', i] = a['chromosome'].sum()
    
#### IMPORTING THE ANNOTATION FOR GENES IN VERSION 6 ####
v6_genes = pd.read_csv('/research/projects/chlamydomonas/MAexpression/analysis/annotation/v6_genes.csv', delimiter = '\t')

#### IMPORTING TOTAL BASE COUNT PER CHROMOSOME ####
base_counts = pd.read_csv('/research/projects/chlamydomonas/MAexpression/data/genome_info/chrom_base_count', skiprows=2, delimiter = '\t', index_col = 'chromosome')

for trials in range(1000):
    simulated_DEGs = pd.DataFrame()
    for sample in CC2344 + CC2931:
        #### CREATING SIMULATED GENES ####
        simulated_genes = pd.DataFrame(index = [e for e in range(len(DEGs_dic[sample]))], columns = ['start', 'end', 'chromosome', 'gene_start', 'gene_end', 'gene_length', 'gene_midpoint', 'min_distance_to_mutations'])

        #### RANDOMLY CHOOSING GENES ACROSS THE GENOME ####
        for b in range(len(DEGs_dic[sample])):
            chrom = np.random.choice(list(base_counts.index.values), 1, replace = True)
            location = np.random.randint(1, base_counts.at[chrom[0], 'base_counts'])
            simulated_genes.at[b, 'chromosome'] = chrom[0]
            simulated_genes.at[b, 'start'] = location
            simulated_genes.at[b, 'end'] = location + 1

        #### ASSIGNING THE GENE NAME TO EACH SIMULATED GENE ####
        simulated_genes['gene'] = 'intergenic'
        for i in range(len(DEGs_dic[sample])):
            specific_chrom = v6_genes.loc[v6_genes['seq_id'] == simulated_genes.at[i, 'chromosome']]
            for a in list(specific_chrom.index.values):
                if specific_chrom.at[a, 'start'] <= simulated_genes.at[i, 'start'] and specific_chrom.at[a, 'end'] >= simulated_genes.at[i, 'start']:
                    simulated_genes.at[i, 'gene'] = specific_chrom.at[a, 'attributes']
                    simulated_genes.at[i, 'gene_start'] = specific_chrom.at[a, 'start']
                    simulated_genes.at[i, 'gene_end'] = specific_chrom.at[a, 'end']
        if simulated_genes.at[i, 'gene'] == 'intergenic':
            simulated_genes.at[i, 'gene_start'] = simulated_genes.at[i, 'start']
            simulated_genes.at[i, 'gene_end'] = simulated_genes.at[i, 'end']
        simulated_DEGs[sample] = simulated_genes['gene']
    simulated_DEGs = simulated_DEGs.replace('intergenic', np.nan, regex = True)

    # #### IDENTIFYING MUTATIONS ENRICHED IN SIMULATED DEGS ####
    for i in CC2344 + CC2931:
        gene_list = simulated_DEGs[i].dropna().values.tolist()
        spec_mut = mutations.loc[mutations['sample'] == i].groupby('gene').count().reset_index()
        unique_mut = spec_mut['gene'].values.tolist()
        intersect = list(set(gene_list).intersection(set(unique_mut)))
        a = spec_mut.loc[spec_mut['gene'].isin(intersect)]
        enriched_DEGs.at['simulated_' + str(trials), i] = a['chromosome'].sum()
enriched_DEGs.to_csv('/research/projects/chlamydomonas/MAexpression/analysis/enriched_DEGs.csv', sep = '\t', index = True, header = True)
