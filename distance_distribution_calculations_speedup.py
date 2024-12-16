import numpy as np
from scipy.spatial import distance
import linecache
from scipy.stats import mannwhitneyu
import pandas as pd
from pybedtools import BedTool
import os
from functools import partial
import sys



def get_coordinates(file_path, file_format='cif'):
    '''
    It returns the corrdinate matrix V (N,3) of a .pdb or .mmcif file.
    The main problem of this function is that coordiantes are not always in
    the same column position of a .pdb or .mmcif file. Do changes appropriatelly,
    in case that the data aren't stored correctly.

    Input:
    file (str): the path of the .pdb or .mmcif file.

    Ouput:
    V (np.array): the matrix of coordinates
    '''
    V = list()

    with open(file_path, 'r') as f:
        lines = f.readlines()
        for line in lines:
            if line.startswith("CONNECT") or line.startswith("END") or line.startswith("TER"):
                break
            if line.startswith("ATOM"):
                col = line.split(" ")
                col = list(filter(None, col))

                x, y, z = float(col[10]), float(col[11]), float(col[12])
                V.append([x, y, z])

    return np.array(V)

def average_binned_distance_matrix(dir_path, Nmodels, file_format='cif'):
    '''
    Modified to return coordinates_array instead of distances_matrix
    '''
    coordinates_list = []
    cif_models = [os.path.join(dir_path, fname) for fname in os.listdir(dir_path) if fname.endswith('cif')]
    
    if len(cif_models) != Nmodels:
        raise ValueError(f"Number of .cif files ({len(cif_models)}) of: \n {dir_path} \ndoes not match the expected number of models ({Nmodels}).")
    
    for model_path in cif_models:
        V = get_coordinates(model_path, file_format)
        coordinates_list.append(V)
    
    coordinates_array = np.stack(coordinates_list, axis=0)  # Shape: (Nmodels, Nbeads, 3)
    first_bin = int(linecache.getline(cif_models[0].replace('.cif', '.txt'), 2).strip().split()[-1])
    last_bin = first_bin + (coordinates_array.shape[1] - 1) * 1000
    return coordinates_array, first_bin, last_bin

def set_bin_number_for_models(region_start_ref, region_start_mod, genes, enhancers, df_enhancers_region_ref,df_genes_region_ref, df_enhancers_region_mod, df_genes_region_mod):
    
    full_genes_for_region = genes.loc[(genes.index.isin(df_genes_region_ref['name'])) & (genes.index.isin(df_genes_region_mod['name']))].astype({'gene_start_mod': 'int32', 'gene_end_mod': 'int32'})
    full_enhancers_for_region = enhancers.loc[(enhancers.index.isin(df_enhancers_region_ref['name'])) & (enhancers.index.isin(df_enhancers_region_mod['name']))].astype({'enh_start_mod': 'int32', 'enh_end_mod': 'int32'})
    full_enhancers_for_region['enh_center_position_ref'] = (full_enhancers_for_region['enh_start_ref'] + full_enhancers_for_region['enh_end_ref']) // 2
    full_enhancers_for_region['enh_center_position_mod'] = (full_enhancers_for_region['enh_start_mod'] + full_enhancers_for_region['enh_end_mod']) // 2

    full_enhancers_for_region['enh_model_position_ref'] = (full_enhancers_for_region['enh_center_position_ref'] - region_start_ref) // 1000 + 1
    full_enhancers_for_region['enh_model_position_mod'] = (full_enhancers_for_region['enh_center_position_mod'] - region_start_mod) // 1000 + 1

    full_enhancers_for_region['enh_model_coloring_start_ref'] = (full_enhancers_for_region['enh_start_ref'] - region_start_ref) // 1000 + 1
    full_enhancers_for_region['enh_model_coloring_end_ref'] = (full_enhancers_for_region['enh_end_ref'] - region_start_ref) // 1000 + 1
    full_enhancers_for_region['enh_model_coloring_start_mod'] = (full_enhancers_for_region['enh_start_mod'] - region_start_mod) // 1000 + 1
    full_enhancers_for_region['enh_model_coloring_end_mod'] = (full_enhancers_for_region['enh_end_mod'] - region_start_mod) // 1000 + 1

    full_genes_for_region.loc[full_genes_for_region['gene_strand'] == '+', 'gene_model_position_ref'] = (full_genes_for_region['gene_start_ref'] - region_start_ref) // 1000 + 1
    full_genes_for_region.loc[full_genes_for_region['gene_strand'] == '-', 'gene_model_position_ref'] = (full_genes_for_region['gene_end_ref'] - region_start_ref) // 1000 + 1

    full_genes_for_region['gene_model_coloring_start_ref'] = (full_genes_for_region['gene_start_ref'] - region_start_ref) // 1000 + 1
    full_genes_for_region['gene_model_coloring_end_ref'] = (full_genes_for_region['gene_end_ref'] - region_start_ref) // 1000 + 1

    full_genes_for_region.loc[full_genes_for_region['gene_strand'] == '+', 'gene_model_position_mod'] = (full_genes_for_region['gene_start_mod'] - region_start_mod) // 1000 + 1
    full_genes_for_region.loc[full_genes_for_region['gene_strand'] == '-', 'gene_model_position_mod'] = (full_genes_for_region['gene_end_mod'] - region_start_mod) // 1000 + 1

    full_genes_for_region['gene_model_coloring_start_mod'] = (full_genes_for_region['gene_start_mod'] - region_start_mod) // 1000 + 1
    full_genes_for_region['gene_model_coloring_end_mod'] = (full_genes_for_region['gene_end_mod'] - region_start_mod) // 1000 + 1

    full_genes_for_region = full_genes_for_region.astype({'gene_model_position_ref': 'int32', 'gene_model_position_mod': 'int32'})

    return full_enhancers_for_region, full_genes_for_region

def find_loci_distance_opt(gene_pos, enh_pos, coordinates_array):
    gene_coords = coordinates_array[:, gene_pos - 1, :]  # Shape: (Nmodels, 3)
    enh_coords = coordinates_array[:, enh_pos - 1, :]    # Shape: (Nmodels, 3)
    distances = np.linalg.norm(gene_coords - enh_coords, axis=1)  # Shape: (Nmodels,)
    return distances


def MWtest_distances_opt(df):
    grouped = df.groupby('sample')['distance']
    data_ref = grouped.get_group('reference')
    data_mut = grouped.get_group('mutated')
    result = mannwhitneyu(data_ref, data_mut, alternative='two-sided')
    return result.pvalue, result.statistic

def process_model_data(path, models_in_ensemble):
    return average_binned_distance_matrix(path, models_in_ensemble)

def run_distance_for_region(path_to_models_ref, path_to_models_mod, models_in_ensemble):
    coordinates_array_ref, first_bin_ref, last_bin_ref = process_model_data(path_to_models_ref, models_in_ensemble)
    coordinates_array_mod, first_bin_mod, last_bin_mod = process_model_data(path_to_models_mod, models_in_ensemble)
    return coordinates_array_ref, coordinates_array_mod, first_bin_ref, last_bin_ref, first_bin_mod, last_bin_mod


def prepare_genes_enhancer_per_region(region_chr, first_bin_coord_ref, last_bin_coord_ref, first_bin_coord_mod,
                                      last_bin_coord_mod, enhancers_bed_ref, enhancers_bed_mod, genes_bed_ref,
                                      genes_bed_mod):
    # wziete z plikow smoot.txt - pierwszy i ostatni+1000 (smooth=1000)
    region_chr_ref = region_chr
    region_start_ref = first_bin_coord_ref
    region_end_ref = last_bin_coord_ref
    region_chr_mod = region_chr
    region_start_mod = first_bin_coord_mod
    region_end_mod = last_bin_coord_mod

    region_bed_ref = BedTool([(region_chr_ref, region_start_ref, region_end_ref)])
    region_bed_mod = BedTool([(region_chr_mod, region_start_mod, region_end_mod)])
    enhancers_region_ref = enhancers_bed_ref.intersect(region_bed_ref, f=1)
    genes_region_ref = genes_bed_ref.intersect(region_bed_ref, f=1)
    enhancers_region_mod = enhancers_bed_mod.intersect(region_bed_mod, f=1)
    genes_region_mod = genes_bed_mod.intersect(region_bed_mod, f=1)

    df_enhancers_region_ref = enhancers_region_ref.to_dataframe()
    df_genes_region_ref = genes_region_ref.to_dataframe()

    df_enhancers_region_mod = enhancers_region_mod.to_dataframe()
    df_genes_region_mod = genes_region_mod.to_dataframe()

    return df_enhancers_region_ref, df_genes_region_ref, df_enhancers_region_mod, df_genes_region_mod

def path_to_data_region(data_dir_path, c, p, region_ref, region_mod):
    return f'{data_dir_path}/{c}/{p}/models3D_{c}_{p}/models3D_{c}_{p}/results_{c}_{p}_{region_ref}', f'{data_dir_path}/{c}/{p}/models3D_{c}_{p}/models3D_{c}_{p}_mod/results_{c}_{p}_{region_ref}'


def process_row(row_gene, coordinates_array_ref, coordinates_array_mod, region_chr_ref, region_start_ref, region_end_ref, region_chr_mod, region_start_mod, region_end_mod):
    distances_gene_enh = find_loci_distance_opt(row_gene['gene_model_position_ref'], row_gene['enh_model_position_ref'], coordinates_array_ref)
    distances_gene_enh_mut = find_loci_distance_opt(row_gene['gene_model_position_mod'], row_gene['enh_model_position_mod'], coordinates_array_mod)
    df_ref = pd.DataFrame(distances_gene_enh, columns=['distance'])
    df_mut = pd.DataFrame(distances_gene_enh_mut, columns=['distance'])
    df_ref['sample'] = 'reference'
    df_mut['sample'] = 'mutated'
    df = pd.concat([df_ref, df_mut])
    mwh_pvalue, mwh_statistic = MWtest_distances_opt(df)
    avg_dist_ref = distances_gene_enh.mean()
    avg_dist_mut = distances_gene_enh_mut.mean()
    avg_dist_sub = abs(avg_dist_ref - avg_dist_mut)
    number_bins_ref = coordinates_array_ref.shape[1]
    number_bins_mod = coordinates_array_mod.shape[1]
    return [region_chr_ref, region_start_ref, region_end_ref, region_chr_mod, region_start_mod, region_end_mod] + list(row_gene) + [avg_dist_ref, avg_dist_mut, avg_dist_sub, mwh_pvalue, mwh_statistic, number_bins_ref, number_bins_mod]



def run_calculate_distance(c, p, region_ref, region_mod, models_in_ensemble, data_dir_path):
    

    dir_path_ref = f'{data_dir_path}/{c}/{p}/models3D_{c}_{p}/models3D_{c}_{p}/results_{c}_{p}_{region_ref}'
    dir_path_mod = f'{data_dir_path}/{c}/{p}/models3D_{c}_{p}/models3D_{c}_{p}_mod/results_{c}_{p}_{region_mod}'



    enhancers = pd.read_csv(f'{data_dir_path}/{c}/{p}/enhanceratlas2_liftovered_hg38_filtered_by_chrom_in_regions_{c}_{p}_with_converted_regions.tsv',sep='\t')
    genes = pd.read_csv(f'{data_dir_path}/{c}/{p}/gencode.v40.annotation_genes_converted_in_regions_{c}_{p}_with_mod_regions_labelled.tsv',sep='\t')

    
    genes_bed_ref = BedTool.from_dataframe(genes.loc[genes['gene_affected_by_svs'].isna()].reset_index()[['gene_chr_ref', 'gene_start_ref', 'gene_end_ref', 'index']])
    genes_bed_mod = BedTool.from_dataframe(genes.loc[genes['gene_affected_by_svs'].isna()].reset_index()[['gene_chr_mod', 'gene_start_mod', 'gene_end_mod', 'index']].astype({'gene_start_mod': 'int32', 'gene_end_mod': 'int32'}))

    
    
    enhancers.rename(columns={'affected_by': 'enh_affected_by_svs', 'chr': 'enh_chr_ref', 'start': 'enh_start_ref', 'end': 'enh_end_ref', 'chr_mod': 'enh_chr_mod', 'start_mod': 'enh_start_mod', 'end_mod': 'enh_end_mod'}, inplace=True)

    
    enhancers_bed_ref = BedTool.from_dataframe(enhancers.loc[enhancers['enh_affected_by_svs'].isna()].reset_index()[['enh_chr_ref', 'enh_start_ref', 'enh_end_ref', 'index']])
    enhancers_bed_mod = BedTool.from_dataframe(enhancers.loc[enhancers['enh_affected_by_svs'].isna()].reset_index()[['enh_chr_mod', 'enh_start_mod', 'enh_end_mod', 'index']].astype({'enh_start_mod': 'int32', 'enh_end_mod': 'int32'}))

    path_to_models_ref, path_to_models_mod = path_to_data_region(data_dir_path, c, p, region_ref, region_mod)

    print('distance calculations')
    # avg_mat_ref, distances_matrix_ref, avg_mat_mod, distances_matrix_mod, first_bin_coord_ref, last_bin_coord_ref, first_bin_coord_mod, last_bin_coord_mod = run_distance_for_region(
    #     path_to_models_ref, path_to_models_mod, models_in_ensemble)
    
    coordinates_array_ref, coordinates_array_mod, first_bin_coord_ref, last_bin_coord_ref, first_bin_coord_mod, last_bin_coord_mod = run_distance_for_region(
        path_to_models_ref, path_to_models_mod, models_in_ensemble)
    
    region_chr, region_start, region_end = region_ref.split('_')
    
    
    df_enhancers_region_ref, df_genes_region_ref, df_enhancers_region_mod, df_genes_region_mod = prepare_genes_enhancer_per_region(
    region_chr, first_bin_coord_ref, last_bin_coord_ref, first_bin_coord_mod, last_bin_coord_mod, enhancers_bed_ref,
    enhancers_bed_mod, genes_bed_ref, genes_bed_mod)

    
    full_enhancers_for_region, full_genes_for_region = set_bin_number_for_models(first_bin_coord_ref,
                                                                             first_bin_coord_mod, genes, enhancers,
                                                                             df_enhancers_region_ref,
                                                                             df_genes_region_ref,
                                                                             df_enhancers_region_mod,
                                                                             df_genes_region_mod)
    
    
    region_chr_ref, region_start_ref, region_end_ref = region_ref.split('_')
    region_chr_mod, region_start_mod, region_end_mod = region_mod.split('_')

    # Dodaj kolumny dla pozycji TSS genu i odległości do enhancera przed pętlą
    full_genes_for_region['gene_TSS_pos'] = np.where(full_genes_for_region['gene_strand'] == '+', full_genes_for_region['gene_start_ref'], full_genes_for_region['gene_end_ref'])
    full_enhancers_for_region['enh_center_pos'] = full_enhancers_for_region['enh_center_position_ref']


    full_enhancers_for_region['enh_loci'] = (
    full_enhancers_for_region['enh_chr_ref'].astype(str) + "_" +
    full_enhancers_for_region['enh_start_ref'].astype(str) + "_" +
    full_enhancers_for_region['enh_end_ref'].astype(str))
    
    full_genes_for_region['key'] = 1
    full_enhancers_for_region['key'] = 1
    df_merged = pd.merge(full_genes_for_region, full_enhancers_for_region, on='key').drop('key', axis=1)
    
    # Obliczenie absolutnej wartości różnicy pozycji
    df_merged['enh_tSS_distance'] = (df_merged['gene_TSS_pos'] - df_merged['enh_center_pos']).abs()

    # Filtracja DataFrame z wieloma warunkami
    df_filtered = df_merged[
    (df_merged['enh_tSS_distance'] <= 1000000) &
    (df_merged['gene_model_position_ref'] != df_merged['enh_model_position_ref']) &
    (df_merged['gene_model_position_mod'] != df_merged['enh_model_position_mod'])
    ]

    df_results = pd.DataFrame(columns=['region_chr_ref', 'region_start_ref', 'region_end_ref', 'region_chr_mod', 'region_start_mod','region_end_mod'] + list(df_filtered.columns) + ['average_dist_ref', 'average_dist_mut', 'average_dist_sub', 'mwh_pvalue', 'mwh_statistics','number_bins_ref', 'number_bins_mod'])

    custom_process_row = partial(
        process_row,
        coordinates_array_ref=coordinates_array_ref,
        coordinates_array_mod=coordinates_array_mod,
        region_chr_ref=region_chr_ref,
        region_start_ref=region_start_ref,
        region_end_ref=region_end_ref,
        region_chr_mod=region_chr_mod,
        region_start_mod=region_start_mod,
        region_end_mod=region_end_mod)
    
    print('now is time for distance distributions tests')
    results = df_filtered.apply(custom_process_row, axis=1)
    
    
    df_results = pd.DataFrame([r for r in results if r is not None], 
                          columns=['region_chr_ref', 'region_start_ref', 'region_end_ref', 'region_chr_mod', 'region_start_mod','region_end_mod'] + list(df_filtered.columns) + ['average_dist_ref', 'average_dist_mut', 'average_dist_sub', 'mwh_pvalue', 'mwh_statistics','number_bins_ref', 'number_bins_mod'])
    
    
    df_results.rename(columns={'gene_id': 'gene_ensemble_ID', 'gene_name': 'gene_hgnc_ID', 'score': 'enh_score'}, inplace=True)


    data_dir_path_to_results = f'{data_dir_path}/{c}/{p}/distance_distributions_{c}_{p}/distance_distributions_{c}_{p}_{region_ref}'
    
    os.makedirs(data_dir_path_to_results, exist_ok = True)
    
    
    df_results[['gene_chr_mod', 'gene_start_mod', 'gene_end_mod', 'gene_model_coloring_start_mod',
            'gene_model_coloring_end_mod', 'gene_hgnc_ID']].drop_duplicates().to_csv(f'{data_dir_path_to_results}/genes_beads_annotation_mod.tsv', sep='\t', index=False)


    df_results_filtered = df_results[
        ['gene_chr_ref', 'gene_start_ref', 'gene_end_ref', 'gene_strand', 'gene_ensemble_ID', 'gene_hgnc_ID',
         'gene_type', 'enh_chr_ref', 'enh_start_ref',
         'enh_end_ref', 'enh_score', 'enh_loci',
         'average_dist_ref', 'average_dist_mut', 'average_dist_sub', 'mwh_pvalue', 'mwh_statistics']]

    df_results.to_csv(f'{data_dir_path_to_results}/results_{region_ref}_loci_distance_distribution.tsv', sep='\t', index=False)

    df_results_filtered.to_csv(f'{data_dir_path_to_results}/results_{region_ref}_loci_distance_distribution_info.tsv', sep='\t',
                               index=False)
    




if __name__ == "__main__":
    
    c = sys.argv[1]
    p = sys.argv[2]
    region_ref = sys.argv[3]
    region_mod = sys.argv[4]
    models_in_ensemble = int(sys.argv[5])
    data_dir_path = sys.argv[6] 
    
#     c = 'H1ESC'
#     p = 'Nean'
#     region_ref = 'chr8_49410110_53686185'
#     region_mod = 'chr8_49410110_53686185'
#     models_in_ensemble = 100
    
#     data_dir_path = '/mnt/raid/altai_faramir/postdoc_Cambridge_24/data/preprocessed_interaction'
    print(f"Time for: {c}, {p} in region: {region_ref}")
    run_calculate_distance(c, p, region_ref, region_mod, models_in_ensemble, data_dir_path)
    






