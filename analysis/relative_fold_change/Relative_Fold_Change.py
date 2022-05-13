import json
import pandas as pd
import numpy as np
import random
import math
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.ticker import MaxNLocator
from matplotlib.gridspec import GridSpec
import seaborn as sns
import scipy
from scipy import stats
from scipy.stats import mannwhitneyu
from scipy.stats import iqr
stats.junk = lambda chisq, df: stats.chi2.sf(chisq, df)
import csv
import gffpandas.gffpandas as gffpd
import scikit_posthocs as sp

#### LIST OF SAMPLE NAMES ####
CC2344 = ['CC2344-ANC', "CC2344-L1", "CC2344-L2", "CC2344-L3", "CC2344-L4", "CC2344-L5", "CC2344-L6", "CC2344-L7", "CC2344-L8", "CC2344-L9", "CC2344-L10", "CC2344-L11", "CC2344-L12", "CC2344-L13", "CC2344-L14", "CC2344-L15"]
CC2931 = ["CC2931-ANC", "CC2931-L1", "CC2931-L2", "CC2931-L3", "CC2931-L4", "CC2931-L5", "CC2931-L6", "CC2931-L7", "CC2931-L9", "CC2931-L10", "CC2931-L11", "CC2931-L13", "CC2931-L14", "CC2931-L15"]

CC2344_list = [ i + '-rep1' for i in CC2344] + [ i + '-rep2' for i in CC2344] + [ i + '-rep3' for i in CC2344]
CC2931_list = [ i + '-rep1' for i in CC2931] + [ i + '-rep2' for i in CC2931] + [ i + '-rep3' for i in CC2931]
    
CC2344_rep = [ i + '-rep1' for i in CC2344[1:]] + [ i + '-rep2' for i in CC2344[1:]] + [ i + '-rep3' for i in CC2344[1:]]
CC2931_rep = [ i + '-rep1' for i in CC2931[1:]] + [ i + '-rep2' for i in CC2931[1:]] + [ i + '-rep3' for i in CC2931[1:]]

#### DATAFRAME RECORDING THE GENERATION TIME PER SAMPLE ####
dic_gen = {'CC2344-L1': 912.356113, 'CC2344-L10': 917.129696, 'CC2344-L11': 889.5859554, 'CC2344-L12': 950.0552184, 'CC2344-L13': 961.4186064,
           'CC2344-L14': 931.447801, 'CC2344-L15': 946.6643063, 'CC2344-L2': 923.1078072, 'CC2344-L3': 1000.469526, 'CC2344-L4': 808.9505794,
           'CC2344-L5': 957.6380465, 'CC2344-L6': 970.6307256, 'CC2344-L7': 990.9451516, 'CC2344-L8': 1009.966123, 'CC2344-L9': 901.0619061, 
           'CC2931-L1': 1050.109001, 'CC2931-L10': 1097.978141, 'CC2931-L11': 1021.13559, 'CC2931-L13': 1041.362593, 
           'CC2931-L14': 1016.111493, 'CC2931-L15': 1052.540951, 'CC2931-L2': 1056.765369, 'CC2931-L3': 1000.399127, 'CC2931-L4': 1011.411706,
           'CC2931-L5': 993.8603657, 'CC2931-L6': 1083.095655, 'CC2931-L7': 1067.34507, 'CC2931-L9': 1079.236285}
generations = pd.Series(dic_gen)

#### DATAFRAME RECORDING THE NUMBER OF MUTATIONS PER SAMPLE ####
dic_mut = {'CC2344-L1': 396, 'CC2344-L10': 59, 'CC2344-L11': 46, 'CC2344-L12': 74, 'CC2344-L13': 49, 'CC2344-L14': 46, 'CC2344-L15': 53, 
           'CC2344-L2': 80, 'CC2344-L3': 63, 'CC2344-L4': 24, 'CC2344-L5': 68, 'CC2344-L6': 38, 'CC2344-L7': 45, 'CC2344-L8': 75, 'CC2344-L9': 27, 
           'CC2931-L1': 89, 'CC2931-L10': 87, 'CC2931-L11': 85, 'CC2931-L13': 97, 'CC2931-L14': 79, 'CC2931-L15': 141, 'CC2931-L2': 123,
           'CC2931-L3': 52, 'CC2931-L4': 100, 'CC2931-L5': 335, 'CC2931-L6': 84, 'CC2931-L7': 72, 'CC2931-L9': 113}
mutations = pd.Series(dic_mut)

#### TAKING THE AVERAGE GENERATION TIME AND NUMBER OF MUTATIONS ####
CC2344_gen_mean = sum([dic_gen[i] for i in CC2344[1:]])/len([dic_gen[i] for i in CC2344[1:]])
CC2931_gen_mean = sum([dic_gen[i] for i in CC2931[1:]])/len([dic_gen[i] for i in CC2931[1:]])

CC2344_mut_mean = sum([dic_mut[i] for i in CC2344[1:]])/len([dic_mut[i] for i in CC2344[1:]])
CC2931_mut_mean = sum([dic_mut[i] for i in CC2931[1:]])/len([dic_mut[i] for i in CC2931[1:]])

""
#### OPENING FILES ####
CC2931_rel = pd.read_csv('/research/projects/chlamydomonas/MAexpression/analysis/raw_counts/CC2931_raw', delimiter = '\t', index_col = 'Unnamed: 0')
CC2344_rel = pd.read_csv('/research/projects/chlamydomonas/MAexpression/analysis/raw_counts/CC2344_raw', delimiter = '\t', index_col = 'Unnamed: 0')

for i in ['CC2344_rel_1', 'CC2931_rel_1', 'CC2344_sum_per_mut', 'CC2931_sum_per_mut']:
    exec('{} = pd.DataFrame()'.format(i))
    
#### DROPPING ROWS WITH ZERO ANCESTRAL COUNTS ####
CC2344_rel = CC2344_rel[CC2344_rel['CC2344-ANC'] != 0]
CC2931_rel = CC2931_rel[CC2931_rel['CC2931-ANC'] != 0]

#### TAKING THE RELATIVE FOLD CHANGE ####
for i in list(CC2344_rel.columns) + list(CC2931_rel.columns):
    if 'CC2931' in i:
        CC2931_rel[i] = CC2931_rel[i]/CC2931_rel['CC2931-ANC']
    if 'CC2344' in i:
        CC2344_rel[i] = CC2344_rel[i]/CC2344_rel['CC2344-ANC']

for i in CC2344[1:] + CC2931[1:]:
    if 'CC2931' in i:
        CC2931_rel[i] = CC2931_rel[[i + '-rep1', i + '-rep2', i + '-rep3']].mean(axis = 1)
        CC2931_rel_1[i] = CC2931_rel[i]
    if 'CC2344' in i:
        CC2344_rel[i] = CC2344_rel[[i + '-rep1', i + '-rep2', i + '-rep3']].mean(axis = 1)
        CC2344_rel_1[i] = CC2344_rel[i]
        
CC2344_rel_1.set_index(CC2344_rel.index.values, inplace = True)
CC2931_rel_1.set_index(CC2931_rel.index.values, inplace = True)

#### EXPORTING RELATIVE FOLD VALUES ####
CC2344_rel.to_csv('/research/projects/chlamydomonas/MAexpression/analysis/raw_counts/CC2344_rel', sep = '\t', index = True, header = True)
CC2931_rel.to_csv('/research/projects/chlamydomonas/MAexpression/analysis/raw_counts/CC2931_rel', sep = '\t', index = True, header = True)
    
#### MEAN/MEDIAN/VARIANCE/SUM OF GENES PER SAMPLE ####
CC2344_rel_1['mean'] = CC2344_rel_1[CC2344[1:]].mean(axis = 1) 
CC2931_rel_1['mean'] = CC2931_rel_1[CC2931[1:]].mean(axis = 1)

CC2344_rel_1['median'] = CC2344_rel_1[CC2344[1:]].median(axis = 1) 
CC2931_rel_1['median'] = CC2931_rel_1[CC2931[1:]].median(axis = 1)

CC2344_rel_1['var'] = CC2344_rel_1[CC2344[1:]].var(axis = 1)
CC2931_rel_1['var'] = CC2931_rel_1[CC2931[1:]].var(axis = 1)

CC2344_rel_1['sum'] = CC2344_rel_1[CC2344[1:]].sum(axis = 1)
CC2931_rel_1['sum'] = CC2931_rel_1[CC2931[1:]].sum(axis = 1)

#### RELATIVE FOLD CHANGE PER OBSERVED AND EXPECTED MUTATION ####
for i in CC2344[1:]:
    CC2344_sum_per_mut.at[i, 'sum'] = CC2344_rel_1[i].sum()
    CC2344_sum_per_mut.at[i, 'sum_per_mut'] = CC2344_sum_per_mut.at[i, 'sum']/dic_mut[i]
    CC2344_sum_per_mut.at[i, 'sum_per_exp_mut'] = CC2344_sum_per_mut.at[i, 'sum']/CC2344_mut_chi.at[i, 'exp_mut']
for i in CC2931[1:]:
    CC2931_sum_per_mut.at[i, 'sum'] = CC2931_rel_1[i].sum()
    CC2931_sum_per_mut.at[i, 'sum_per_mut'] = CC2931_sum_per_mut.at[i, 'sum']/dic_mut[i]
    CC2931_sum_per_mut.at[i, 'sum_per_exp_mut'] = CC2931_sum_per_mut.at[i, 'sum']/CC2931_mut_chi.at[i, 'exp_mut']

#### CONVERTING DATAFRAME USING PD.MELT ####
CC2344_sum_per_mut = CC2344_sum_per_mut.reset_index()
CC2344_RELATIVE_FOLD = pd.melt(CC2344_sum_per_mut, id_vars = 'index', value_vars = ['sum_per_mut', 'sum_per_exp_mut'])
CC2344_RELATIVE_FOLD['sample'] = 'CC2344'
CC2344_RELATIVE_FOLD['fold'] = 'relative fold'

CC2931_sum_per_mut = CC2931_sum_per_mut.reset_index()
CC2931_RELATIVE_FOLD = pd.melt(CC2931_sum_per_mut, id_vars = 'index', value_vars = ['sum_per_mut', 'sum_per_exp_mut'])
CC2931_RELATIVE_FOLD['sample'] = 'CC2931'
CC2931_RELATIVE_FOLD['fold'] = 'relative fold'

# ## SEABORN
# a = sns.distplot(CC2931_sum_per_mut['sum_per_mut'])
# a = sns.distplot(CC2931_sum_per_mut['sum_per_exp_mut'])
# plt.xlim(0, 450)
# plt.title('CC2931 - RELATIVE FOLD SUM PER MUTATION')
# plt.savefig('/research/projects/chlamydomonas/MAexpression/analysis/DEGs/images/CC2931_RELATIVE_FOLD_SUM_PER_MUT.pdf', format = 'pdf', dpi = 150, bbox_inches = 'tight')

# b = sns.distplot(CC2344_sum_per_mut['sum_per_mut'])
# b = sns.distplot(CC2344_sum_per_mut['sum_per_exp_mut'])
# plt.xlim(0, 1000)
# plt.title('CC2344 - RELATIVE FOLD SUM PER MUTATION')
# plt.savefig('/research/projects/chlamydomonas/MAexpression/analysis/DEGs/images/CC2344_RELATIVE_FOLD_SUM_PER_MUT.pdf', format = 'pdf', dpi = 150, bbox_inches = 'tight')

#############################################################################################################################
# ### INVESTIGATING THE OUTLIER GENE ####
gene1 = CC2344_rel_1.loc[CC2344_rel_1['sum'] == CC2344_rel_1['sum'].max()]
gene2 = CC2344_rel_1.loc[CC2344_rel_1['var'] == CC2344_rel_1['var'].max()]
gene3 = CC2344_rel_1.loc[CC2344_rel_1['mean'] == CC2344_rel_1['mean'].max()]

geneID = pd.read_csv('/research/projects/chlamydomonas/MAexpression/data/genome_info/v5_to_v6_liftover/preliminary_v6_Cre_liftover.tsv', delimiter = '\t', header = None)
geneID.columns = ['v6', 'v5']
UNIPROT = pd.read_csv('/research/projects/chlamydomonas/MAexpression/data/genome_info/annotation_package/edited_anno/UNIPROT_BY_GENENAME', delimiter = '\t')

gene = geneID.loc[geneID['v6'] == gene1.index.item()]
gene['v5'].item()
# No function associated with the gene 

############################################################################################################################################
# ### MUTATIONAL VARIANCE 

for i in ['CC2931_Vm', 'CC2344_Vm']:
    exec('{} = pd.DataFrame()'.format(i))

#### FINDING THE RELATIVE FOLD CHANGE OVER MA GENERATIONS (VARIANCE PER GENE) ####
CC2344_Vm['Vm_per_gen'] = CC2344_rel[CC2344[1:]].var(axis = 1)/CC2344_gen_mean
CC2931_Vm['Vm_per_gen'] = CC2931_rel[CC2931[1:]].var(axis = 1)/CC2931_gen_mean
CC2344_Vm['log10_Vm_per_gen'] = np.log10(CC2344_Vm['Vm_per_gen'])
CC2931_Vm['log10_Vm_per_gen'] = np.log10(CC2931_Vm['Vm_per_gen'])

#### FINDING THE RELATIVE FOLD CHANGE OVER MUTATIONS (VARIANCE PER GENE) ####
CC2344_Vm['Vm_per_mut'] = CC2344_rel[CC2344[1:]].var(axis = 1)/CC2344_mut_mean
CC2931_Vm['Vm_per_mut'] = CC2931_rel[CC2931[1:]].var(axis = 1)/CC2931_mut_mean
CC2344_Vm['log10_Vm_per_mut'] = np.log10(CC2344_Vm['Vm_per_mut'])
CC2931_Vm['log10_Vm_per_mut'] = np.log10(CC2931_Vm['Vm_per_mut'])

###################################################################
# ### MUTATIONAL VARIANCE - ISOLATING HIGH/LOW EXPRESSION GENES
for i in ['Vm_per_mut', 'log10_Vm_per_mut', 'abval_Vm_per_mut', 'abval_log10_Vm_per_mut', 'Vm_per_gen', 'log10_Vm_per_gen', 'abval_Vm_per_gen', 'abval_log10_Vm_per_gen']:
    #### CREATING DATAFRAME ####
    for a in ['CC2344_high', 'CC2344_low', 'CC2931_high', 'CC2931_low']:
        exec('{} = pd.DataFrame()'.format(a))
    #### ISOLATING HIGH/LOW EXPRESSION GENES ####
    CC2344_high[i] = CC2344_Vm[i].loc[CC2344_hi_exp].values.tolist()
    CC2344_high['expression'] = 'high'
    CC2344_high['id'] = 'CC2344'
    CC2344_low[i] = CC2344_Vm[i].loc[CC2344_low_exp].values.tolist()
    CC2344_low['expression'] = 'low'
    CC2344_low['id'] = 'CC2344'

    CC2931_high[i] = CC2931_Vm[i].loc[CC2931_hi_exp].values.tolist()
    CC2931_high['expression'] = 'high'
    CC2931_high['id'] = 'CC2931'
    CC2931_low[i] = CC2931_Vm[i].loc[CC2931_low_exp].values.tolist()
    CC2931_low['expression'] = 'low'
    CC2931_low['id'] = 'CC2931'
    rel_per_gen = pd.concat([CC2344_high, CC2344_low, CC2931_high, CC2931_low], axis = 0)
    #### SEABORN ####
    plt.figure()
    sns.boxplot(data = rel_per_gen, x = 'expression', y = i, hue = 'id')
    plt.savefig('/research/projects/chlamydomonas/MAexpression/analysis/relative_fold_change/boxplots_mut_var/' + i + '.pdf', format = 'pdf', dpi = 150, bbox_inches = 'tight')
    ## MANN WHITNEY U TEST ##
    CC2344_ts, CC2344_p = mannwhitneyu(CC2344_high[i].values.tolist(), CC2344_low[i].values.tolist())
    CC2931_ts, CC2931_p = mannwhitneyu(CC2931_high[i].values.tolist(), CC2931_low[i].values.tolist())
    print(CC2344_ts, CC2344_p, CC2931_ts, CC2931_p)

############################################################################################################################################
# ### PLOTS

#### HISTOGRAM OF THE MEAN/MEDIAN ####
sns.histplot(CC2344_rel_1['mean'])
sns.histplot(CC2344_rel_1['median'])
plt.xlim(0, 12)
plt.title('CC2344 - The mean/median relative fold change of all samples')
plt.savefig('/research/projects/chlamydomonas/MAexpression/analysis/relative_fold_change/CC2344_allsamples_mean_and_median.pdf', format = 'pdf', dpi = 150, bbox_inches = 'tight')

plt.figure()
sns.histplot(CC2931_rel_1['mean'])
sns.histplot(CC2931_rel_1['median'])
plt.xlim(0, 12)
plt.title('CC2931 - The mean/median relative fold change of all samples')
plt.savefig('/research/projects/chlamydomonas/MAexpression/analysis/relative_fold_change/CC2931_allsamples_mean_and_median.pdf', format = 'pdf', dpi = 150, bbox_inches = 'tight')

#### HISTOGRAM OF THE RELATIVE FOLD OF EACH GENE PER MA LINE OF A GENOTYPE ####
CC2344_RF = CC2344_rel_1.drop(['mean', 'median', 'var', 'sum'], axis = 1)
CC2931_RF = CC2931_rel_1.drop(['mean', 'median', 'var', 'sum'], axis = 1)

for i in CC2344[1:]:
    plt.figure(i)
    sns.histplot(CC2344_RF[i])
    plt.xlim(0, 12)
    plt.axvline(1, color = 'red')
    plt.title(i)
    plt.xlabel('relative fold')
    plt.savefig('/research/projects/chlamydomonas/MAexpression/analysis/relative_fold_change/histogram/histogram_' + i + '.pdf', format = 'pdf', dpi = 150, bbox_inches = 'tight')

for i in CC2931[1:]:
    plt.figure(i)
    sns.histplot(CC2931_RF[i])
    plt.xlim(0, 12)
    plt.axvline(1, color = 'red')
    plt.title(i)
    plt.xlabel('relative fold')
    plt.savefig('/research/projects/chlamydomonas/MAexpression/analysis/relative_fold_change/histogram/histogram_' + i + '.pdf', format = 'pdf', dpi = 150, bbox_inches = 'tight')
    
#### GRAPHING THE RELATIVE FOLD CHANGE FOR ALL SAMPLES ####
sum_FC = pd.DataFrame(index = CC2344[1:] + CC2931[1:])
for i in list(sum_FC.index.values):
    if 'CC2344' in i:
        sum_FC.at[i, 'sum'] = CC2344_rel_1[i].sum()
    if 'CC2931' in i:
        sum_FC.at[i, 'sum'] = CC2931_rel_1[i].sum()
sum_FC = sum_FC.reset_index()
plt.figure()
sns.barplot(data = sum_FC, y = 'index', x = 'sum')
plt.ylabel('MA samples')
plt.title('RELATIVE CHANGE FOR ALL SAMPLES')
plt.savefig('/research/projects/chlamydomonas/MAexpression/analysis/relative_fold_change/all_samples_FC.pdf', format = 'pdf', dpi = 150, bbox_inches = 'tight')

############################################################################################################################################
# ### BOOTSTRAPPING 

#### CC2344 ####
for i in ['Vm_per_mut', 'log10_Vm_per_mut', 'Vm_per_gen', 'log10_Vm_per_gen']:
    CC2344_high[i] = CC2344_Vm[i].loc[CC2344_hi_exp].values.tolist()
    CC2344_low[i] = CC2344_Vm[i].loc[CC2344_low_exp].values.tolist()
    x = CC2344_high[i].values.tolist()
    y = CC2344_low[i].values.tolist()
    
    #### OBSERVED TEST STATISTIC ####
    obs_median_high = np.median(x)
    obs_median_low = np.median(y)
    obs_test_statistic = obs_median_low - obs_median_high

    test_statistics = []
    for e in range(50000):
        a = np.random.choice(x + y, len(x + y), replace=True, p=None)
        median_high = np.median(a[len(x)])
        median_low = np.median(a[len(x):len(x+y)])
        test_statistics.append(median_low - median_high)

    histplot = pd.DataFrame({'median_test_stat':test_statistics})
    plt.figure()
    mut_var_median = sns.histplot(data = histplot, x = 'median_test_stat')
    plt.axvline(obs_test_statistic, color = 'red')
    plt.savefig('/research/projects/chlamydomonas/MAexpression/analysis/relative_fold_change/bootstrapping/CC2344_' + i + '.pdf', format = 'pdf', dpi = 150, bbox_inches = 'tight')
    pvalue = (histplot.loc[histplot['median_test_stat'] >= obs_test_statistic].size)/histplot.size
    mut_var_median
    print(pvalue)
    
#### CC2931 ####
for i in ['Vm_per_mut', 'log10_Vm_per_mut', 'Vm_per_gen', 'log10_Vm_per_gen']:
    CC2931_high[i] = CC2931_Vm[i].loc[CC2931_hi_exp].values.tolist()
    CC2931_low[i] = CC2931_Vm[i].loc[CC2931_low_exp].values.tolist()
    x = CC2931_high[i].values.tolist()
    y = CC2931_low[i].values.tolist()
    
    #### OBSERVED TEST STATISTIC ####
    obs_median_high = np.median(x)
    obs_median_low = np.median(y)
    obs_test_statistic = obs_median_low - obs_median_high

    test_statistics = []
    for e in range(50000):
        a = np.random.choice(x + y, len(x + y), replace=True, p=None)
        median_high = np.median(a[len(x)])
        median_low = np.median(a[len(x):len(x+y)])
        test_statistics.append(median_low - median_high)

    histplot = pd.DataFrame({'median_test_stat':test_statistics})
    plt.figure()
    mut_var_median = sns.histplot(data = histplot, x = 'median_test_stat')
    plt.axvline(obs_test_statistic, color = 'red')
    plt.savefig('/research/projects/chlamydomonas/MAexpression/analysis/relative_fold_change/bootstrapping/CC2931_' + i + '.pdf', format = 'pdf', dpi = 150, bbox_inches = 'tight')
    pvalue = (histplot.loc[histplot['median_test_stat'] >= obs_test_statistic].size)/histplot.size
    mut_var_median
    print(pvalue)

###############################################################################
# ### RELATIVE FOLD CHANGE AMONG FPKM BINS

#### OPENING FILES ####
for i in ['CC2344_fpkm_bins', 'CC2931_fpkm_bins']:
    exec('{} = pd.DataFrame()'.format(i))
    
CC2344 = ['CC2344-ANC', "CC2344-L1", "CC2344-L2", "CC2344-L3", "CC2344-L4", "CC2344-L5", "CC2344-L6", "CC2344-L7", "CC2344-L8", "CC2344-L9", "CC2344-L10", "CC2344-L11", "CC2344-L12", "CC2344-L13", "CC2344-L14", "CC2344-L15"]
CC2931 = ["CC2931-ANC", "CC2931-L1", "CC2931-L2", "CC2931-L3", "CC2931-L4", "CC2931-L5", "CC2931-L6", "CC2931-L7", "CC2931-L9", "CC2931-L10", "CC2931-L11", "CC2931-L13", "CC2931-L14", "CC2931-L15"]
    
CC2931_raw = pd.read_csv('/research/projects/chlamydomonas/MAexpression/analysis/raw_counts/CC2931_rel', delimiter = '\t', index_col = 'Unnamed: 0')
CC2344_raw = pd.read_csv('/research/projects/chlamydomonas/MAexpression/analysis/raw_counts/CC2344_rel', delimiter = '\t', index_col = 'Unnamed: 0')

#### DETERMINING THE RELATIVE FOLD SUM OF GENES WITHIN EACH FPKM PERCENTILE BIN (EXCLUDING ANCESTOR) ####
for i in CC2344_bins.keys():
    genes = CC2344_bins[i]
    genes = set.intersection(set(genes), set(list(CC2344_raw.index.values)))
    CC2344_sum = CC2344_raw[CC2344[1:]].loc[genes].sum(axis = 1)
    CC2344_fpkm_bins.at[i, 'bins'] = CC2344_sum.sum()
    
for i in CC2931_bins.keys():
    genes = CC2931_bins[i]
    genes = set.intersection(set(genes), set(list(CC2931_raw.index.values)))
    CC2931_sum = CC2931_raw[CC2931[1:]].loc[genes].sum(axis = 1)
    CC2931_fpkm_bins.at[i, 'bins'] = CC2931_sum.sum()
    
#### PLOTTING THE SUM OF RELATIVE FOLD CHANGES PER FPKM PERCENTILE BINS #####
CC2344_fpkm_bins = CC2344_fpkm_bins.reset_index()
CC2931_fpkm_bins = CC2931_fpkm_bins.reset_index()
sns.barplot(data = CC2344_fpkm_bins, x = 'index', y = 'bins', color = 'blue')
sns.barplot(data = CC2931_fpkm_bins, x = 'index', y = 'bins', color = 'red', alpha = 0.5)
plt.xlabel('FPKM percentile bins')
plt.ylabel('Total sum of relative fold changes')
plt.axhline(0, color = 'black')
plt.title('Sum of relative fold changes per FPKM percentile bins [Blue = CC2344, Red = CC2931]')
plt.savefig('/research/projects/chlamydomonas/MAexpression/analysis/relative_fold_change/sum_of_relative_fold_per_FPKM_percentile_bins.pdf', format = 'pdf', dpi = 150, bbox_inches = 'tight')
