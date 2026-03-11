import scanpy as sc # imports pandas internally
import matplotlib.pyplot as plt
from pathlib import Path

def compute_qc_metrics(adata):
    mt_gene = adata.var_names.str.startswith("MT-") # Creates a true or false if var_names begins with MT-
    adata.var["mt"] = mt_gene # Creates a column called mt with the boolean values.
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True) # Calculates metrics, selects variables with True in the mt column, and updates the adata.obs object.
    return adata

def plot_qc_metrics(adata, save_dir="results/figures"): # Save directory defaults to results/figures.
    Path(save_dir).mkdir(parents=True, exist_ok=True) # Makes the directory if it doesnt already exist
    sc.pl.violin(adata, keys=["n_genes_by_counts", "total_counts", "pct_counts_mt"], jitter=0.4, multi_panel=True, show=False) # Makes a violin plot with n_gene_by_counts, total_counts, percentage_counts of mt, with some added jitter to make the points clear instead of overlayed.
    plt.savefig(f"{save_dir}/violin_plot.png") # Saves it in the save_dir.
    print (f"Violin plots created at {save_dir}!")
    return adata
    
def filter_cells(adata, min_n_gene_counts = 200, max_n_gene_counts = 2500, max_pct_mt = 20):
    previous_cells = adata.n_obs
    sc.pp.filter_cells(adata, min_genes=min_n_gene_counts)
    sc.pp.filter_cells(adata, max_genes= max_n_gene_counts) # In future, will add dedicated doublet detection tools: e.g. Scrublet 
    adata = adata[adata.obs["pct_counts_mt"] < max_pct_mt].copy() # Removes cells where % mitochondrial DNA > 20.
    removed_cells = previous_cells - adata.n_obs
    print(f"Total cells before filtering: {previous_cells} \n Total cells removed: {removed_cells}. \n Total cells after filtering: {adata.n_obs}")
    return adata

def filter_genes(adata, min_cells = 3):
    previous_genes = adata.n_vars
    sc.pp.filter_genes(adata, min_cells=min_cells) # Filters out genes who show up in fewer than 3 cells.
    print(f"Number of genes before filtering: {previous_genes}. \n Number of genes after filtering: {adata.n_vars}.")
    return adata

def run_qc(adata): # Runs the whole QC pipeline at once.
    print("Running QC pipeline...")
    adata = compute_qc_metrics(adata)
    print("QC metrics computed...")
    adata = plot_qc_metrics(adata)
    print("Plots produced...")
    adata = filter_cells(adata)
    print("Cells filtered...")
    adata = filter_genes(adata)
    print("Genes filtered...")
    print("QC pipeline complete!")
    return adata


if __name__ == "__main__":
    from .data_loading import load_pbmc
    adata = load_pbmc() 
    run_qc(adata)

# To run, python -m src.scrna_pipeline.qc in your bash terminal.