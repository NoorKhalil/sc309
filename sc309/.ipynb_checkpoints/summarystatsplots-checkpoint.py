def summarystatsplots(SingleCellData):
    
    sc.pl.violin(SingleCellData, ['n_genes', 'n_counts', 'percent_mito'],
             jitter=0.4, multi_panel=True) # Violin plots of number of genes, total counts and mitochondrial percent
    
    sc.pl.scatter(SingleCellData, x='n_counts', y='percent_mito') # Scatter plot between total counts and mitochondrial percent
    sc.pl.scatter(SingleCellData, x='n_counts', y='n_genes') # Scatter plot between total counts and number of genes
    
    sc.pl.highly_variable_genes(SingleCellData) # Plotting highly variable genes (mean expression and dispersion)