import subprocess, re
import numpy as np
import pandas as pd
import scanpy as sc


class Preprocessor():
    def __init__(self, raw: str, mode: str) -> None:
        self.raw = raw
        self.mode = mode
        self.matrix = None

    def count(self, ref_dir: str) -> None:
        # Working directory is /mnt/data/cellar
        # User uploaded raw fastq data sits in /data/raw
        # Prepared reference data sits in /data/ref
        # Temporary unpack location is /data/tmp
        base_name = self.raw.split('/')[-1]
        tmp_location = '/mnt/data/cellar/data/tmp'
        subprocess.run(['mkdir', tmp_location])
        try:
            subprocess.run(['tar', '-xvf', self.raw, '-C', tmp_location], check=True, capture_output=True)
        except subprocess.CalledProcessError as e:
            print(f'Cannot extract tar file, error:\n {e.stderr}')
            raise e
        try:
            # Extract the base name of the tar file
            base_dir = re.match(r'^([^.]*)\..+$', base_name).group(1)
            subprocess.run(['cellranger', 'count', f'--id=run_count_{base_dir}',
            f'--fastqs={tmp_location}/{base_dir}',
            # TODO: Implement sample selection with "--sample"
            f'--transcriptome={ref_dir}'], check=True, capture_output=True)
            self.matrix = f'/mnt/data/cellar/run_count_{base_dir}/outs/filtered_feature_bc_matrix/'
        except subprocess.CalledProcessError as e:
            print(f'Cannot make fastq, error:\n {e.stderr}')
            raise e
        subprocess.run(['rm', '-rf', tmp_location])

    def read(self):
        if self.matrix is None:
            raise ValueError('Matrix is not initialized')
        adata = sc.read_10x_mtx(
            self.matrix,                             # the directory with the `.mtx` file
            var_names='gene_symbols',                # use gene symbols for the variable names (variables-axis index)
            cache=True)                              # write a cache file for faster subsequent reading
        if self.mode == 'rna_seq':
            adata.var_names_make_unique()       # this is unnecessary if using `var_names='gene_ids'` in `sc.read_10x_mtx`
            sc.pp.filter_cells(adata, min_genes=200)
            sc.pp.filter_genes(adata, min_cells=3)
            adata.var['mt'] = adata.var_names.str.startswith('MT-')  # annotate the group of mitochondrial genes as 'mt'
            sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
            adata = adata[adata.obs.n_genes_by_counts < 2500, :]
            adata = adata[adata.obs.pct_counts_mt < 5, :]
            sc.pp.normalize_total(adata, target_sum=1e4)
            sc.pp.log1p(adata)
            sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
            adata = adata[:, adata.var.highly_variable]
            sc.pp.scale(adata, max_value=10)
            adata.write('/mnt/data/cellar/preprocessed/pbmc_1k_v3_filtered.h5ad')


preprocessor = Preprocessor(raw='/mnt/data/cellar/data/raw/pbmc_1k_v3_fastqs.tar', mode='rna_seq')
# TODO: Choose reference data to use according to the kind of raw data
preprocessor.count('/mnt/data/cellar/data/ref/refdata-cellranger-GRCh38-3.0.0')
preprocessor.read()