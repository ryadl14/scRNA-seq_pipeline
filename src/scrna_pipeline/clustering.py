import scanpy as sc
import matplotlib.pyplot as plt
from pathlib import Path
import igraph


def run_leiden (adata, resolution = 0.5):
    sc.tl.leiden(adata=adata, resolution=resolution, flavor="igraph", n_iterations=2, directed=False)
    print ("Leiden algorithm successfully ran!")
    return adata

def plot_clusters(adata, save_dir="results/figures"):
    Path(save_dir).mkdir(parents=True, exist_ok=True) # Makes the directory if it doesnt already exist
    sc.pl.umap(adata=adata, color=["leiden"], show=False)
    plt.savefig(f"{save_dir}/leiden_plot.png") # Saves it in the save_dir.
    print (f"Leiden plot created at {save_dir}!")
    return
    

def run_clustering(adata):
    print("Starting Leiden algorithm...")
    adata = run_leiden(adata)
    print("Plotting clusters...")
    plot_clusters(adata)
    return adata

if __name__ == "__main__": # Runs the whole pipeline up until this point.
    from .data_loading import load_pbmc
    from .qc import run_qc
    from .preprocessing import run_preprocessing
    from .dimensionality import run_dimensionality_reduction
    adata = load_pbmc() 
    adata = run_qc(adata)
    adata = run_preprocessing(adata)
    adata = run_dimensionality_reduction(adata)
    adata = run_clustering (adata)