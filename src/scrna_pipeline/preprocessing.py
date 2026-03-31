import scanpy as sc

def normalise_data (adata, target_sum = 1e4): # Cell counts are scaled to add up to 10,000
    sc.pp.normalize_total(adata, target_sum=target_sum)
    print ("Normalisation complete!")
    return adata

def log_transform_data (adata): # Log transform all values.
    sc.pp.log1p(adata)
    print ("Log transformation complete!")
    return adata

def hvg_selection (adata, n_top_genes = 2000): # Selects the top 2000 highly variable genes.
    sc.pp.highly_variable_genes(adata, n_top_genes=n_top_genes, subset=False) # Does not subset them, just flags the top 2000 in a boolean "highly_variable" row.
    print (f"{n_top_genes} highly variable genes selected!")
    return adata

def run_preprocessing(adata): # Runs all the functions sequentially. 
    print ("Starting normalisation...")
    adata = normalise_data(adata)
    print ("Starting log transformation...")
    adata = log_transform_data(adata)
    print ("Starting highly variable gene selection...")
    adata = hvg_selection(adata)
    return adata

if __name__ == "__main__": # Runs the whole pipeline up until this point.
    from .data_loading import load_pbmc
    from .qc import run_qc
    adata = load_pbmc() 
    adata = run_qc(adata)
    run_preprocessing(adata)

