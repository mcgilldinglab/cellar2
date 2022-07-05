import subprocess, re
import numpy as np
import pandas as pd
import scanpy as sc


class Preprocessor():
    def __init__(self, raw: str, mode: str) -> None:
        self.base_dir = re.match(r'^([^.]*)\..+$', raw.split('/')[-1]).group(1)
        self.raw = raw
        self.mode = mode
        self.matrix = None

    def count(self, ref_dir: str, probe_dir: str=None, image_dir: str=None, lib_dir: str=None) -> None:
        # Working directory is /mnt/data/cellar
        # User uploaded raw fastq data sits in /data/raw
        # Prepared reference data sits in /data/ref
        # Temporary unpack location is /data/tmp
        tmp_location = '/mnt/data/cellar/data/tmp'
        subprocess.run(['mkdir', tmp_location])
        try:
            subprocess.run(['tar', '-xvf', self.raw, '-C', tmp_location], check=True)
        except subprocess.CalledProcessError as e:
            print(f'Cannot extract tar file, error:\n {e.stderr}')
            raise e
        try:
            # Extract the base name of the tar file
            if self.mode == 'rna_seq':
                subprocess.run(['cellranger', 'count', f'--id=output_{self.base_dir}',
                f'--fastqs={tmp_location}/{self.base_dir}',
                # TODO: Implement sample selection with "--sample"
                f'--transcriptome={ref_dir}'], check=True, capture_output=True)
                self.matrix = f'/mnt/data/cellar/output_{self.base_dir}/outs/filtered_feature_bc_matrix/'
            elif self.mode == 'spatial_ffpe':
                subprocess.run(['spaceranger', 'count', f'--id=output_{self.base_dir}',
                f'--transcriptome={ref_dir}',
                f'--probe-set={probe_dir}',
                f'--fastqs={tmp_location}/{self.base_dir}',
                f'--image={image_dir}',
                # TODO: Find a way to get `slide` and `area` either from the raw data or from the user input
                f'--slide=V10L13-021',
                f'--area=B1', '--reorient-images'], check=True)
                self.matrix = f'/mnt/data/cellar/output_{self.base_dir}/outs/filtered_feature_bc_matrix/'
            elif self.mode == 'spatial_ff':
                subprocess.run(['spaceranger', 'count', f'--id=output_{self.base_dir}',
                f'--transcriptome={ref_dir}',
                f'--fastqs={tmp_location}/{self.base_dir}',
                f'--image={image_dir}',
                # TODO: Find a way to get `slide` and `area` either from the raw data or from the user input
                f'--slide=V10N30-322',
                f'--area=A1'], check=True, capture_output=True)
                self.matrix = f'/mnt/data/cellar/output_{self.base_dir}/outs/filtered_feature_bc_matrix/'

            elif self.mode == 'multi':
                subprocess.run(['cellranger-arc', 'count', f'--id=output_{self.base_dir}',
                f'--reference={ref_dir}',
                f'--libraries={lib_dir}'], check=True, capture_output=True)
                self.matrix = f'/mnt/data/cellar/output_{self.base_dir}/outs/filtered_feature_bc_matrix/'
                
        except subprocess.CalledProcessError as e:
            print(f'Cannot run count pipeline, error:\n {e.stderr}')
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
            adata = adata[adata.obs.pct_counts_mt < 20, :]
            sc.pp.normalize_total(adata, target_sum=1e4)
            sc.pp.log1p(adata)
            sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
            adata = adata[:, adata.var.highly_variable]
            sc.pp.scale(adata, max_value=10)
            adata.write(f'/mnt/data/cellar/preprocessed/{self.base_dir}.h5ad')
        elif self.mode == 'spatial_ffpe' or self.mode == 'spatial_ff':
            adata.var_names_make_unique()
            adata.var["mt"] = adata.var_names.str.startswith("MT-")
            sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True)
            sc.pp.filter_cells(adata, min_counts=5000)
            sc.pp.filter_cells(adata, max_counts=35000)
            adata = adata[adata.obs["pct_counts_mt"] < 20]
            sc.pp.filter_genes(adata, min_cells=10)
            sc.pp.normalize_total(adata, inplace=True)
            sc.pp.log1p(adata)
            sc.pp.highly_variable_genes(adata, flavor="seurat", n_top_genes=2000)
            # TODO: Change output file name
            adata.write(f'/mnt/data/cellar/preprocessed/{self.base_dir}.h5ad')

        elif self.mode == 'multi':
            adata.var_names_make_unique()
            adata.var["mt"] = adata.var_names.str.startswith("MT-")
            sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True)
            sc.pp.filter_cells(adata, min_counts=int(input("Please enter the minimum count value for filtering: ")))
            sc.pp.filter_cells(adata, max_counts=int(input("Please enter the maximum count value for filtering: ")))
            
            adata = adata[adata.obs["pct_counts_mt"] < int(input("Please enter the minimum percentage of mitochondrial genes expressed for filtering: "))]
            sc.pp.filter_genes(adata, min_cells=int(input("Please enter the minimum number of cells the gene is detected in for filtering: ")))
            sc.pp.normalize_total(adata, inplace=True)
            sc.pp.log1p(adata)
            sc.pp.highly_variable_genes(adata, flavor="seurat", n_top_genes=int(input("Please enter the number of top genes to be considered in filtering: ")))
            # TODO: Change output file name
            adata.write(f'/mnt/data/cellar/preprocessed/{self.base_dir}.h5ad')
            
        elif self.mode == 'codex':
            pass



preprocessor = Preprocessor(raw='/mnt/data/cellar/data/raw/pbmc_granulocyte_sorted_3k_fastqs.tar', mode='multi')
# TODO: Choose reference data to use according to the kind of raw data
preprocessor.count(ref_dir='/mnt/data/cellar/data/ref/refdata-cellranger-arc-atac-human',
                   lib_dir='/mnt/data/cellar/data/raw/pbmc_granulocyte_sorted_3k_library.csv')
preprocessor.read()
