import matplotlib.pyplot as plt
import scanpy as sc
from pathlib import Path

def run_pca (adata, n_comps = 50):
    sc.tl.pca(adata, n_comps = n_comps, mask_var="highly_variable")
    print ("PCA complete!")
    return adata

def plot_elbow (adata, save_dir="results/figures"):
    Path(save_dir).mkdir(parents=True, exist_ok=True) # Makes the directory if it doesnt already exist
    sc.pl.pca_variance_ratio(adata, n_pcs=50, log=True, show=False)
    plt.savefig(f"{save_dir}/pca_plot.png") # Saves it in the save_dir.
    print (f"PCA plot created at {save_dir}!")
    return

def compute_neighbours(adata, n_pcs = 10, n_neighbors = 10):
    sc.pp.neighbors(adata=adata, n_neighbors=n_neighbors, n_pcs=n_pcs)
    print (f"{n_neighbors} neighbours computed!")
    return adata

def run_umap(adata):
    sc.tl.umap(adata=adata, random_state=42)
    print ("UMAP ran!")
    return adata

def run_dimensionality_reduction(adata): # Runs all the functions sequentially. 
    print ("Starting PCA...")
    adata = run_pca(adata)
    print ("Starting elbow plot...")
    plot_elbow(adata)
    print ("Computing neighbours...")
    adata = compute_neighbours(adata)
    print ("Running UMAP...")
    adata = run_umap(adata)
    return adata

if __name__ == "__main__": # Runs the whole pipeline up until this point.
    from .data_loading import load_pbmc
    from .qc import run_qc
    from .preprocessing import run_preprocessing
    adata = load_pbmc() 
    adata = run_qc(adata)
    adata = run_preprocessing(adata)
    adata = run_dimensionality_reduction(adata)
