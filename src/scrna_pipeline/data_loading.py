import scanpy as sc

def load_pbmc(data_dir: str = "data/raw"): # Default location is under data/raw, but user can state where they want the data loaded.
    adata = sc.datasets.pbmc3k()
    adata.raw = adata.copy()
    print (f"Data has been loaded at {data_dir}!")
    print (f"Loaded: {adata.n_obs} cells x {adata.n_vars} genes")
    return adata


if __name__ == "__main__":
    adata = load_pbmc()