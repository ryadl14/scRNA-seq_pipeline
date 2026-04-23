import scanpy as sc
import matplotlib.pyplot as plt
from pathlib import Path


def find_marker_genes (adata):
    # Groups by the leiden column, using raw counts to preserve expression differences for use use in Wilcoxons test.
    sc.tl.rank_genes_groups(adata=adata, groupby="leiden", method="wilcoxon")
    print ("Marker genes found!")
    return adata


def annotate_clusters (adata, cluster_to_celltype):
    adata.obs["cell_type"] = adata.obs["leiden"].map(cluster_to_celltype)
    print ("Annotation complete!")
    return adata

def plot_annotation (adata, save_dir="results/figures"):
    Path(save_dir).mkdir(parents=True, exist_ok=True) # Makes the directory if it doesnt already exist
    sc.pl.umap(adata=adata, color=["cell_type"], show=False)
    plt.savefig(f"{save_dir}/annotation_plot.png", bbox_inches="tight") # Saves it in the save_dir.
    print (f"Annotation plot created at {save_dir}!")
    return

def run_annotation (adata):
    default_annotation = { # A default PBMC dictionary.
    "0": "CD4+ T cells",
    "1": "CD14+ Monocytes",
    "2": "CD4+ T cells",
    "3": "B cells",
    "4": "CD8+ T cells",
    "5": "NK cells",
    "6": "FCGR3A+ Monocytes",
    "7": "Dendritic cells",
    "8": "Megakaryocytes",
    }
    print ("Starting to find marker genes...")
    adata = find_marker_genes(adata)
    print ("Starting to annotate clusters...")
    adata = annotate_clusters(adata, cluster_to_celltype=default_annotation)
    print ("Starting to plot annotations...")
    plot_annotation(adata)
    return adata

if __name__ == "__main__": # Runs the whole pipeline up until this point.
    from .data_loading import load_pbmc
    from .qc import run_qc
    from .preprocessing import run_preprocessing
    from .dimensionality import run_dimensionality_reduction
    from .clustering import run_clustering
    adata = load_pbmc() 
    adata = run_qc(adata)
    adata = run_preprocessing(adata)
    adata = run_dimensionality_reduction(adata)
    adata = run_clustering (adata)
    adata = run_annotation(adata)