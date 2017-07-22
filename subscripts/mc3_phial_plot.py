def plot_top_muts(df, tissue, num_display, samples):
	n_samples = len(samples)
	df_top_genes = df['Gene'].value_counts().head(num_display)
	genes = df_top_genes.index.tolist()
	hits = df_top_genes.value_counts().head(num_display).index.tolist()

	fig = plt.figure(figsize = (9,6))
	ax = plt.subplot()
	ax.spines["top"].set_visible(False)    
	ax.spines["bottom"].set_visible(False)    
	ax.spines["right"].set_visible(False)    
	ax.spines["left"].set_visible(False)
	plt.tick_params(axis="both", which="both", bottom="off", top="off",
    	labelbottom="on", left="off", right="off", labelleft="on")  

	df_top_genes.plot.barh(stacked=True, width = 0.75, color = tableau10[6])

	ax.set_yticklabels(genes, style = 'italic', size = 13)
	plt.xlabel('Overall number of mutants in gene', fontsize = 13)
	title_ = 'Most Common TARGET-related SNVs / InDels: TCGA mc3 - ' + tissue + ' (n=' + str(n_samples) +')' + '\n Putatively Actionable or Biologically Relevant'
	plt.title(title_ , fontsize = 14)

	for i, v in enumerate(df_top_genes.get_values().tolist()):
		ax.text(v + 0.2, i - 0.2, str(v), color='black')

	return(fig)


def plot_top_muts_cohort(df, tissue, num_display):
	n_samples = len(df.iloc[:,:-3].columns.tolist())
	df_top = df['percent_of_samples'].sort_values(ascending = False).head(num_display)
	genes = df_top.index.tolist()
	percent = df_top.get_values().tolist()

	fig = plt.figure(figsize=(9, 6))

	ax = plt.subplot()
	ax.spines["top"].set_visible(False)    
	ax.spines["bottom"].set_visible(False)    
	ax.spines["right"].set_visible(False)    
	ax.spines["left"].set_visible(False)
	plt.tick_params(axis="both", which="both", bottom="off", top="off",
    	labelbottom="on", left="off", right="off", labelleft="on")  

	df_top.plot.barh(stacked=True, width = 0.75, color = tableau10[6])

	ax.set_yticklabels(genes, style = 'italic', size = 13)
	plt.xlabel('Percentage of samples with mutant', size = 13)
	title_ = 'Most Common TARGET-related SNVs / InDels: TCGA mc3 - ' + tissue + ' (n=' + str(n_samples) +')' + '\n Putatively Actionable or Biologically Relevant'
	plt.title(title_, fontsize = 14)

	for i, v in enumerate(percent):
		ax.text(v + 0.2, i - 0.2, str(v) + ' %', color='black')

	return(fig)


def plot_muts_per_sample(df, tissue):
	n_samples = len(df)
	fig = plt.figure(figsize=(9, 2))

	ax = plt.subplot()
	ax.spines["top"].set_visible(False)    
	ax.spines["bottom"].set_visible(False)    
	ax.spines["right"].set_visible(False)    
	ax.spines["left"].set_visible(False)
	plt.tick_params(axis="both", which="both", bottom="off", top="off",
    	labelbottom="on", left="off", right="off", labelleft="off")  

	bp = pylab.boxplot(df['# TARGET muts'], 
                   vert = False, notch = False, patch_artist = True)
	pylab.setp(bp['boxes'], color='black')
	pylab.setp(bp['whiskers'], color='black')
	pylab.setp(bp['fliers'], marker='*')

	colors = [tableau10[5]]
	for patch, color in zip(bp['boxes'], colors):
		patch.set_facecolor(color)

	for median in bp['medians']:
		median.set(color=tableau10[1], linewidth=2)
    
	plt.ylim([0.80,1.20])

	title_ = 'Putatively Actionable or Biologically Relevant \n SNVs/InDels per Sample - TCGA mc3 ' + tissue + ' (n=' + str(n_samples) +')'
	plt.title(title_, fontsize = 14)
	plt.xlabel('Number of relevant alterations per sample', fontsize = 12)

	plt.grid(True)
	return(fig)

def plot_top_cnv(df, tissue, num_display):
	n_samples = len(df_gistic.iloc[:,3:-9].columns.tolist())
	index_ = df['num_of_samples'].sort_values(ascending = False).head(num_display).index.tolist()
	genes = df.loc[index_, 'Gene Symbol']
	num = df['num_of_samples'].sort_values(ascending = False).head(num_display).get_values().tolist()
	num_amp = df.loc[index_, 'num_of_samples_amp']
	num_del = df.loc[index_, 'num_of_samples_del']

	fig = plt.figure(figsize=(9, 6))

	ax = plt.subplot()
	ax.spines["top"].set_visible(False)    
	ax.spines["bottom"].set_visible(False)    
	ax.spines["right"].set_visible(False)    
	ax.spines["left"].set_visible(False)
	plt.tick_params(axis="both", which="both", bottom="off", top="off",
		labelbottom="on", left="off", right="off", labelleft="on")  

	num_amp.plot(kind='barh', stacked = True, width = 0.75, color = tableau10[2],
		label = 'Amplification')
	num_del.plot(kind='barh', stacked = True, width = 0.75, color = tableau10[0], 
		label = 'Deletion', left = num_amp)

	ax.set_yticklabels(genes, style = 'italic', size = 12)
	plt.xlabel('Number of samples with mutant', fontsize = 12)
	title_ = 'Most Common TARGET-related CNVs: TCGA mc3 - ' + tissue + ' (n=' + str(n_samples) +')' + '\n Putatively Actionable or Biologically Relevant'
	plt.title(title_, fontsize = 14)

	for i, v in enumerate(num):
		ax.text(v + 0.2, i - 0.2, str(v), color='black')

	plt.legend(frameon = False)
    
	plt.grid(False)
	return(fig)

def plot_top_cnv_cohort(df, tissue, num_display):
	n_samples = len(df_gistic.iloc[:,3:-9].columns.tolist())
	index_ = df['percent_of_samples'].sort_values(ascending = False).head(num_display).index.tolist()
	genes = df.loc[index_, 'Gene Symbol']
	percent = df['percent_of_samples'].sort_values(ascending = False).head(num_display).get_values().tolist()
	percent_amp = df.loc[index_, 'percent_of_samples_amp']
	percent_del = df.loc[index_, 'percent_of_samples_del']

	fig = plt.figure(figsize=(9, 6))

	ax = plt.subplot()
	ax.spines["top"].set_visible(False)    
	ax.spines["bottom"].set_visible(False)    
	ax.spines["right"].set_visible(False)    
	ax.spines["left"].set_visible(False)
	plt.tick_params(axis="both", which="both", bottom="off", top="off",
		labelbottom="on", left="off", right="off", labelleft="on")  

	percent_amp.plot(kind='barh', stacked = True, width = 0.75, color = tableau10[2],
		label = 'Amplification')
	percent_del.plot(kind='barh', stacked = True, width = 0.75, color = tableau10[0], 
		label = 'Deletion', left = percent_amp)

	ax.set_yticklabels(genes, style = 'italic', size = 12)
	plt.xlabel('Percentage of samples with mutant', fontsize = 12)
	title_ = 'Most Common TARGET-related CNVs: TCGA mc3 - ' + tissue + ' (n=' + str(n_samples) +')' + '\n Putatively Actionable or Biologically Relevant'
	plt.title(title_, fontsize = 14)

	for i, v in enumerate(percent):
		ax.text(v + 0.2, i - 0.2, str(v) + '%', color='black')

	plt.legend(frameon = False)
    
	plt.grid(False)
	return(fig)

def plot_cnvs_per_sample(df, tissue):
	n_samples = len(df)
	fig = plt.figure(figsize=(9, 2))

	ax = plt.subplot()
	ax.spines["top"].set_visible(False)    
	ax.spines["bottom"].set_visible(False)    
	ax.spines["right"].set_visible(False)    
	ax.spines["left"].set_visible(False)
	plt.tick_params(axis="both", which="both", bottom="off", top="off",
    	labelbottom="on", left="off", right="off", labelleft="off")  

	bp = pylab.boxplot(df['# TARGET cnv'], 
		vert = False, notch = False, patch_artist = True)
	pylab.setp(bp['boxes'], color='black')
	pylab.setp(bp['whiskers'], color='black')
	pylab.setp(bp['fliers'], marker='*')

	colors = [tableau10[4]]
	for patch, color in zip(bp['boxes'], colors):
		patch.set_facecolor(color)

	for median in bp['medians']:
		median.set(color=tableau10[1], linewidth=2)
    
	plt.ylim([0.80,1.20])

	title_ = 'Putatively Actionable or Biologically Relevant \n CNVs per Sample - TCGA mc3 ' + tissue + ' (n=' + str(n_samples) +')'
	plt.title(title_, fontsize = 14)
	plt.xlabel('Number of relevant alterations per sample', fontsize = 12)

	plt.grid(True)
	return(fig)