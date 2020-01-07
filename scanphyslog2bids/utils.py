import click
import os.path as op
import pandas as pd
from glob import glob
from tqdm import tqdm

# from https://github.com/translationalneuromodeling/tapas/blob/7371f120610a25b233f9da6e51013cb077d7879b/PhysIO/code/model/tapas_physio_create_retroicor_regressors.m
NAMES = [
    # modality, cos/sin, order
    'cardiac_cos_00', 'cardiac_sin_00',
    'cardiac_cos_01', 'cardiac_sin_01',
    'cardiac_cos_02', 'cardiac_sin_02',
    'resp_cos_00', 'resp_sin_00',
    'resp_cos_01', 'resp_sin_01',
    'resp_cos_02', 'resp_sin_02',
    'resp_cos_03', 'resp_sin_03',
    # additive (cardiac + resp) or difference (cardiac - resp), cos/sin, order
    'interaction_add_cos_00', 'interaction_add_sin_00',
    'interaction_diff_cos_00', 'interaction_diff_sin_00',
    'hrv', 'rvt'
]

@click.command()
@click.option('--directory', required=True, help='Physiology derivative dir')
def add_header_to_physio_tsv(directory):
    files = sorted(glob(op.join(directory, '**', '*.tsv'), recursive=True))
    for tsv in tqdm(files):
        df = pd.read_csv(tsv, sep='\t', header=None)
        if df.iloc[0, 0] == 'cardiac_cos_00':
            continue

        if df.iloc[:, -1].isna().all():
            # Remove empty col
            df = df.iloc[:, :-1]

        if len(NAMES) != df.shape[1]:
            raise ValueError("Number of columns is not the same as the assumed regressors!")

        df.columns = NAMES
        df.to_csv(tsv, sep='\t', index=False)
