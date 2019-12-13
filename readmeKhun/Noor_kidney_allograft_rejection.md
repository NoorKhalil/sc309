
## Installing package


```python
!pip install -e ../
```


```python
!pwd
```

    /Users/khunzawlatt/Desktop/Python_course/sc309/Notebooks



```python
from sc309 import scanpy_workflow
```


```python
import scanpy as sc
```


```python
sc.settings.verbosity = 3
sc.logging.print_versions()
```

    scanpy==1.4.4.post1 anndata==0.6.22.post1 umap==0.3.10 numpy==1.17.4 scipy==1.3.0 pandas==0.25.3 scikit-learn==0.21.2 statsmodels==0.10.0 python-igraph==0.7.1 louvain==0.6.1


## First, import your single cell data from a text file


```python
#gzip -d is the same as gunzip, but works on Windows

!gzip -d ../Kidney_allo_rejection_python_input.txt.gz
adata = sc.read_text(
    '../Kidney_allo_rejection_python_input.txt')  # the directory with the `.mtx` file
    # use gene symbols for the variable names (variables-axis index)
#Github cannot handle large files, so zipping compresses the large file.
!gzip ../Kidney_allo_rejection_python_input.txt
```


```python

```

#### Do not run unless for GSE dataset


```python
# import gzip
# import shutil
# with gzip.open('../GSE109564_Kidney.biopsy.dge.txt.gz', 'rb') as f_in:
#     with open('../GSE109564_Kidney.biopsy.dge.txt', 'wb') as f_out:
#         shutil.copyfileobj(f_in, f_out)
# # Transporsing the text files in python
# dataframe = pd.read_csv("../GSE109564_Kidney.biopsy.dge.txt",delimiter="\t")
# dataframe = pd.DataFrame.transpose(dataframe)
# dataframe.to_csv("Kidney.biopsy.dge.csv", encoding='utf-8', index=False)
# # Starting analysis
# adata = sc.read_csv(
#     './Kidney.biopsy.dge.csv')  # the directory with the .mtx file
#     # use gene symbols for the variable names (variables-axis index)
```


```python
adata
```




    AnnData object with n_obs × n_vars = 4487 × 20477 



## You can use scanpy to look at the top genes to make sure that your data has been imported correctly


```python
sc.pl.highest_expr_genes(adata, n_top=20)
```

    normalizing by total count per cell
        finished (0:00:00): normalized adata.X and added    'n_counts', counts per cell before normalization (adata.obs)



![png](output_13_1.png)


## Filtering cells and Normalization


```python
scanpy_workflow.filterRegress(adata)
```

    normalizing by total count per cell


    Trying to set attribute `.obs` of view, making a copy.


        finished (0:00:00): normalized adata.X and added    'n_counts', counts per cell before normalization (adata.obs)
    extracting highly variable genes
        finished (0:00:01)
    --> added
        'highly_variable', boolean vector (adata.var)
        'means', float vector (adata.var)
        'dispersions', float vector (adata.var)
        'dispersions_norm', float vector (adata.var)
    regressing out ['n_counts', 'percent_mito']


    //anaconda3/lib/python3.7/site-packages/statsmodels/compat/pandas.py:23: FutureWarning: The Panel class is removed from pandas. Accessing it from the top-level namespace will also be removed in the next version
      data_klasses = (pandas.Series, pandas.DataFrame, pandas.Panel)


        finished (0:00:18)



```python
scanpy_workflow.summarystatsplots(adata)
```


![png](output_16_0.png)



![png](output_16_1.png)



![png](output_16_2.png)



    ---------------------------------------------------------------------------

    AttributeError                            Traceback (most recent call last)

    <ipython-input-9-3bde175908d3> in <module>
    ----> 1 scanpy_workflow.summarystatsplots(adata)
    

    ~/Desktop/Python_course/sc309/sc309/scanpy_workflow.py in summarystatsplots(SingleCellData)
         42     sc.pl.scatter(SingleCellData, x='n_counts', y='n_genes') # Scatter plot between total counts and number of genes
         43 
    ---> 44     sc.pl.highly_variable_genes(SingleCellData) # Plotting highly variable genes (mean expression and dispersion)
         45 
         46 def topclustergenes(SingleCellData):


    ~/.local/lib/python3.7/site-packages/scanpy/plotting/_preprocessing.py in highly_variable_genes(adata_or_result, log, show, save, highly_variable_genes)
         32         result = adata_or_result
         33     if highly_variable_genes:
    ---> 34         gene_subset = result.highly_variable
         35     else:
         36         gene_subset = result.gene_subset


    //anaconda3/lib/python3.7/site-packages/pandas/core/generic.py in __getattr__(self, name)
       5177             if self._info_axis._can_hold_identifiers_and_holds_name(name):
       5178                 return self[name]
    -> 5179             return object.__getattribute__(self, name)
       5180 
       5181     def __setattr__(self, name, value):


    AttributeError: 'DataFrame' object has no attribute 'highly_variable'


## PCA and UMAP


```python
scanpy_workflow.topclustergenes(adata)
```

    computing PCA with n_comps = 50
        finished (0:00:06)
    computing neighbors
        using 'X_pca' with n_pcs = 10
        finished: added to `.uns['neighbors']`
        'distances', distances for each pair of neighbors
        'connectivities', weighted adjacency matrix (0:00:03)
    computing UMAP
        finished: added
        'X_umap', UMAP coordinates (adata.obsm) (0:00:08)
    running Louvain clustering
        using the "louvain" package of Traag (2017)
        finished: found 19 clusters and added
        'louvain', the cluster labels (adata.obs, categorical) (0:00:01)
    ranking genes
        finished: added to `.uns['rank_genes_groups']`
        'names', sorted np.recarray to be indexed by group ids
        'scores', sorted np.recarray to be indexed by group ids
        'logfoldchanges', sorted np.recarray to be indexed by group ids
        'pvals', sorted np.recarray to be indexed by group ids
        'pvals_adj', sorted np.recarray to be indexed by group ids (0:00:42)



![png](output_18_1.png)


## Visualizing genes


```python
scanpy_workflow.visualizegenes(adata)
```


![png](output_20_0.png)



![png](output_20_1.png)



![png](output_20_2.png)



![png](output_20_3.png)



```python

```
