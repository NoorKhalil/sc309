def filterRegress(SingleCellData):
    """
    Filters genes and cells, normalizes the data, determines highly variable genes, calculates percentage of mitochondrial genes, and regresses out the total counts and mitochondrial percentage.
    """
    sc.pp.filter_cells(SingleCellData, min_genes=200) # Filtering cells with gene count < 200
    sc.pp.filter_genes(SingleCellData, min_cells=3) ## Filtering genes expressed in < 3 cells
     
# for each cell compute fraction of counts in mito genes vs. all genes
# the `.A1` is only necessary as X is sparse (to transform to a dense array after summing)
    mito_genes = SingleCellData.var_names.str.startswith('MT-') # Selecting mitochondrial genes; starting with MT
    SingleCellData.obs['percent_mito'] = np.sum(
    SingleCellData[:, mito_genes].X, axis=1) / np.sum(SingleCellData.X, axis=1) # Calculating the percentage of mitochondrial genes in each cell

    SingleCellData.obs['n_counts'] = SingleCellData.X.sum(axis=1) # add the total counts per cell as observations-annotation to adata
    
    SingleCellData = SingleCellData[SingleCellData.obs['n_genes'] < 2500, :] # Filtering cells with gene counts > 2500
    SingleCellData = SingleCellData[SingleCellData.obs['percent_mito'] < 0.05, :] # Filtering cells with mitochondrial percentage > 5%
    
    sc.pp.normalize_per_cell(SingleCellData, counts_per_cell_after=1e4) # Normalization; may smoothen out biological variation
    
    sc.pp.log1p(SingleCellData) # Log-transformation
    
    sc.pp.highly_variable_genes(SingleCellData, min_mean=0.0125, max_mean=3, min_disp=0.5) # Finding the highly variable genes based on level of expression and dispersion
    
    SingleCellData = SingleCellData[:, SingleCellData.var['highly_variable']] # Selecting and keeping only highly variable genes
    sc.pp.regress_out(SingleCellData, ['n_counts', 'percent_mito']) # Regressing out based on total RNA counts and mitochondrial percentage