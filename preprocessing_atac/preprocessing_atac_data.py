import subprocess, re
import numpy as np
import pandas as pd
import scanpy as sc
import episcanpy as epi
import anndata as ad

class Preprocessor():

    def __init__(self, raw:str, mode:str) -> None:
        
        self.base_dir = re.match(r'^([^.]*)\..+$', raw.split('/')[-1]).group(1)
        self.raw = raw
        self.mode = mode
        self.matrix = None


    def count(self, ref_dir: str, probe_dir: str=None, image_dir: str= None) -> None:

        tmp_location = '/mnt/data/cellar/data/tmp'
        subprocess.Popen(['mkdir', tmp_location])


        try:
            subprocess.run(['tar', '-xvf', self.raw, '-C', tmp_location], check = True)

        except subprocess.CalledProcessError as e:
            print(f'Cannot extract tar file, error:\n {e.stderr}')
            raise e

        try:
            if self.mode == 'rna_seq':
                subprocess.run(['cellranger', 'count', f'--id=output_{self.base_dir}',
                f'--fastqs={tmp_location}/{self.base_dir}',
                # TODO: Implement sample selection with "--sample"
                f'--transcriptome={ref_dir}'], check=True)
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
                f'--area=A1'], check=True)
                self.matrix = f'/mnt/data/cellar/output_{self.base_dir}/outs/filtered_feature_bc_matrix      

            elif self.mode == 'ATAC':

                subprocess.Popen(['cellranger-atac', 'count', f'--id=output_{self.base_dir}', f'--fastqs={tmp_location}/{self.base_dir}', f'--reference={ref_dir}'])
                self.matrix = f'/mnt/data/cellar/output_{self.base_dir}/outs/filtered_peak_bc_matrix/matrix.mtx'

        except subprocess.CalledProcessError as e:
            print(f'Cannot run count pipeline, error:\n {e.stderr}')
            raise e

        #subprocess.run(['rm', '-rf', tmp_location])

    def read(self):
        file_path = self.matrix
        name_cells = f'/mnt/data/cellar/output_{self.base_dir}/outs/filtered_peak_bc_matrix/barcodes.tsv'
        names_var = f'/mnt/data/cellar/output_{self.base_dir}/outs/filtered_peak_bc_matrix/peaks.bed'
        adata = epi.pp.read_ATAC_10x(matrix = file_path, cell_names = name_cells, var_names = names_var)
        epi.pp.filter_cells(adata, min_features=1)
        epi.pp.filter_features(adata, min_cells=1)
        adata.obs['log_nb_features'] = [np.log10(x) for x in adata.obs['nb_features']]
        min_features = 1000
        epi.pp.filter_cells(adata, min_features=min_features)
        min_cells = 5
        epi.pp.filter_features(adata, min_cells=min_cells)
        epi.pp.filter_cells(adata, min_features=2000)
        epi.pp.filter_cells(adata, max_features=25000)
        adata.layers['binary'] = adata.X.copy()
        adata.write(f'/mnt/data/cellar/preprocessed/atac_v1_hgmm_1k.h5ad')


        

        




preprocessor = Preprocessor(raw = '/mnt/data/cellar/data/raw/atac_v1_E18_brain_flash_5k_fastqs.tar', mode = 'ATAC')
preprocessor.count(ref_dir='/mnt/data/cellar/data/ref/refdata-cellranger-arc-atac-human')
preprocessor.read()

        
