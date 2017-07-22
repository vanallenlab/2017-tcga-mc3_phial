def rename_gistic(gistic):
	gistic.loc[23312, 'Gene Symbol'] = 'MTG8' # RUNX1 -> MTG8
	gistic.loc[748, 'Gene Symbol'] = 'MYH' # MUTYH -> MYH
	gistic.loc[14388, 'Gene Symbol'] = 'CHK1' # CHEK1 -> CHK1
	gistic.loc[23688, 'Gene Symbol'] = 'CHK2' # CHEK2 -> CHK2
	gistic.loc[22041, 'Gene Symbol'] = 'POLD' # POLD1 -> POLD
	gistic.loc[5530, 'Gene Symbol'] = 'MMSET' # WHSC1 (NSD2) -> MMSET
	gistic.loc[23698, 'Gene Symbol'] = 'EWS' # EWSR1 -> EWS
	return(gistic)

def load_tissue(tissue):
	input_file = 'data/somatic_phial_output/' + tissue + '/' + tissue + '_complete_muts_indels_scna_detailed.txt'
	df = pd.read_csv(input_file, sep = '\t', low_memory = False, index_col = None, comment = '#')
	
	df['CODE'] = tissue
	return(df)

def subset_tissue_samples(df, df_summary, tissue):
	tmp = df_summary[df_summary['CODE'] == tissue]
	df_return = df[df['Tumor_Sample_Barcode'].isin(tmp['SNV sample full'].tolist())]

	n_var = len(df_return)
	n_samples = len(df_return['Tumor_Sample_Barcode'].unique())

	print '....Number of Samples for tissue:', str(n_samples)
	print '....Number of Variants for tissue:', str(n_var)
	print ''
	return(df_return)

def subset_target(df):
	target_categories = ['Potential TARGET Actionability', 'Investigate TARGET Actionability', 'High Priority']
	df_target = df[df['Score_bin'].isin(target_categories)]

	n_var = len(df_target)
	n_samples = len(df_target['Tumor_Sample_Barcode'].unique())
	n_samples_ = len(df['Tumor_Sample_Barcode'].unique())
	frac_target_samples = round((float(n_samples) / float(n_samples_)), 4)

	print '....Number of TARGET related variants:', str(n_var)
	print '....Number of Samples with TARGET related variants:', str(n_samples)
	print '....Fraction of Samples with TARGET related event:', str(frac_target_samples)
	
	return(df_target)

def genehits_across_samples(df, samples):
	genes = df['Gene'].unique().tolist()
	hits_across_samples = pd.DataFrame([], index = genes, columns = samples)

	for sample_ in samples:
		hits_across_samples[sample_] = df[df['Tumor_Sample_Barcode'] == sample_]['Gene'].value_counts()

	genehits_across_samples = hits_across_samples.fillna(0)
	genehits_across_samples['num_of_samples'] = (genehits_across_samples != 0).astype(int).sum(axis=1)
	genehits_across_samples['frac_of_samples'] = genehits_across_samples['num_of_samples'].astype(int) / len(samples)
	genehits_across_samples['percent_of_samples'] = (genehits_across_samples['frac_of_samples']*100).round(2)
	
	return(genehits_across_samples)

def target_hits_per_sample(df, tissue):
	n_samples = len(df.iloc[:,:-3].columns.tolist())
	mut_target_per_sample = (df.iloc[:,:-3] > 0).astype(int).sum()
	df_mut_target_per_sample = pd.DataFrame(mut_target_per_sample, columns = ['# TARGET muts'])
	mean_hits_per_sample = round(df_mut_target_per_sample['# TARGET muts'].mean(),2)

	print '....' + tissue, 'samples have a mean number of ' + str(mean_hits_per_sample) + ' TARGET-related SNV/InDels per sample.'
	return(df_mut_target_per_sample)

def gistic_sample_cut(gistic):
	gistic_short = gistic.iloc[:,3:].rename(columns = lambda x: str(x)[:-16])
	df_gistic = pd.concat([gistic.iloc[:,0:3], gistic_short], axis = 1)
	return(df_gistic)

def subset_samples_gistic(gistic, samples):
	samples_short = [sample_[:-16] for sample_ in samples]
	df_gistic = gistic.loc[:, sorted(samples_short)]
	df_gistic = pd.concat([gistic.iloc[:,0:3], df_gistic], axis = 1)
	
	missing_ids = list(set(samples_short) - set(gistic.iloc[:,3:].columns.tolist()))
	print '....Missing Individual IDs from GISTIC:', missing_ids
	frac_missing = float(len(missing_ids))/ float(len(samples))
	pcnt_missing = round(frac_missing*100, 2)
	print '....Percent of samples missing:', str(pcnt_missing), '%'
	print '....Number of GISTIC samples for tissue:', str(len(samples) - len(missing_ids))

	df_gistic = df_gistic.drop(missing_ids, axis = 1)
	return(df_gistic)

def subset_threshold_gistic(gistic):
	samples = gistic.iloc[:,3:].columns.tolist()
	n_samples = len(samples)
	n_amp = (gistic.iloc[:,3:] == +2).sum(axis = 1).sum()
	n_del = (gistic.iloc[:,3:] == -2).sum(axis = 1).sum()
	n_total = n_amp + n_del

	print '\n'
	print '...', str(n_total), 'TARGET related amplifications or deletions in cohort.'
	print '....', str(n_del), 'TARGET related deletions (-2) in cohort'
	print '....', str(n_amp), 'TARGET related amplifications (+2) in cohort'

	gistic.loc[:,'num_of_samples'] = (gistic.iloc[:,3:(n_samples+3)].abs() == 2).sum(axis=1)
	gistic.loc[:,'frac_of_samples'] = gistic['num_of_samples'].astype(int) / len(samples)
	gistic.loc[:,'percent_of_samples'] = (gistic['frac_of_samples']*100).round(2)

	gistic.loc[:,'num_of_samples_amp'] = (gistic.iloc[:,3:(n_samples+3)] == 2).sum(axis=1)
	gistic.loc[:,'frac_of_samples_amp'] = gistic.loc[:,'num_of_samples_amp'].astype(int) / len(samples)
	gistic.loc[:,'percent_of_samples_amp'] = (gistic.loc[:,'frac_of_samples_amp']*100).round(2)

	gistic.loc[:,'num_of_samples_del'] = (gistic.iloc[:,3:(n_samples+3)] == -2).sum(axis=1)
	gistic.loc[:,'frac_of_samples_del'] = gistic.loc[:,'num_of_samples_del'].astype(int) / len(samples)
	gistic.loc[:,'percent_of_samples_del'] = (gistic.loc[:,'frac_of_samples_del']*100).round(2)

	n_samples_hit = n_samples - (df_gistic.iloc[:,3:-9].abs() == 2).astype(int).sum().value_counts()[0]
	frac_samples_hit = round(float(n_samples_hit)/float(n_samples),4)

	print 'Number of samples with TARGET related CNVs:', str(n_samples_hit)
	print 'Fraction of CNV samples with TARGET related CNVs:', str(frac_samples_hit)

	return(gistic)

def gistic_hits_per_sample(gistic, tissue):
	gistic_hits_per_sample = (gistic.iloc[:,3:-9].abs() == 2).astype(int).sum()
	df_gistic_hits_per_sample = pd.DataFrame(gistic_hits_per_sample, columns = ['# TARGET cnv'])

	mean_hits_per_sample = round(df_gistic_hits_per_sample['# TARGET cnv'].mean(),2)
	print '....' + tissue, 'samples have a mean number of ' + str(mean_hits_per_sample) + ' TARGET-related CNVs per sample.'
	return(df_gistic_hits_per_sample)