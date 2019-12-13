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
    
def summarystatsplots(SingleCellData):
    """
    Draws violin plots for the number of genes, total RNA counts and mitochondrial percent for each cell, scatter plots for relationships between these parameters and calculate and plot highly variable genes
    """
    
    sc.pl.violin(SingleCellData, ['n_genes', 'n_counts', 'percent_mito'],
             jitter=0.4, multi_panel=True) # Violin plots of number of genes, total counts and mitochondrial percent
    
    sc.pl.scatter(SingleCellData, x='n_counts', y='percent_mito') # Scatter plot between total counts and mitochondrial percent
    sc.pl.scatter(SingleCellData, x='n_counts', y='n_genes') # Scatter plot between total counts and number of genes
    
    sc.pl.highly_variable_genes(SingleCellData) # Plotting highly variable genes (mean expression and dispersion)
    
def topclustergenes(SingleCellData)
    """
   Performs PCA, compute UMAP and do louvain clustering. Find the genes for each cluster by comparing with the remaining clusters and plots the top 25 genes for each cluster.
    """

    sc.tl.pca(SingleCellData, svd_solver='arpack') # PCA
    sc.pp.neighbors(SingleCellData, n_neighbors=10, n_pcs=10) # Finding neighbors
    
    sc.tl.umap(SingleCellData) # Computing UMAP
    sc.tl.louvain(SingleCellData) ## Louvain clustering
    
    sc.tl.rank_genes_groups(SingleCellData, 'louvain', method='wilcoxon') # Finding top genes for each cluster
    sc.pl.rank_genes_groups(SingleCellData, n_genes=25, sharey=False) # Plotting top 25 genes
    
def visualizegenes(SingleCellData)
    """
   Plots variance explained by each PC, draws expression gradient plots, violin and dot plots for the canonical genes of the kidney and immune cells
    """

    sc.pl.pca_variance_ratio(SingleCellData, log=True) # Variance explained by PCs
    sc.pl.umap(SingleCellData, color=['louvain','NPHS1', 'NKG7', 'WT1', 'FCER1G', 'TYROBP', 'APOE', 'IL1B', 'CD3D', 'GZMA'], size = 50)
    # Checking the expression of canonical genes for kidney and immune cells
    
    sc.pl.violin(SingleCellData, ['NPHS1', 'NKG7', 'WT1', 'FCER1G', 'TYROBP', 'APOE', 'IL1B', 'CD3D', 'GZMA'], groupby='louvain')
    # Violin plots
    
    ax = sc.pl.dotplot(SingleCellData, ['NPHS1', 'NKG7', 'WT1', 'FCER1G', 'TYROBP', 'APOE', 'IL1B', 'CD3D', 'GZMA'], groupby='louvain')
    # Dot plot