import scanpy as sc
import logging
import pandas as pd
import anndata as ad
import squidpy as sq
import scipy
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns


logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'  # Custom date format
)

def calc_fov_qc_metrics(adata, obs_column, obs_value, fov_size=(0.51, 0.51)):
    # Filter the AnnData object based on the given obs value
    sample_data = adata[adata.obs[obs_column] == obs_value].copy()
    fov_area_mm2 = fov_size[0] * fov_size[1]

    total_cells = sample_data.n_obs
    num_fovs = len(sample_data.obs['fov'].unique())  # Assuming fov_count column exists
    fovs = str(sample_data.obs['fov'].unique())
    cell_density = total_cells / (num_fovs * fov_area_mm2)

    metrics = pd.DataFrame({
        'Sample': [obs_value],
        'FOVs': [fovs],
        'Total Cells': [total_cells],
        'Cell Density (cells/mm²)': [cell_density]
    })
    return metrics

def preprocess(adata, res_dir, nm, save_intermediate):
    logging.info("### Preprocessing Function - 6k variant ###")
    logging.info("## res_dir = %s", res_dir)

    adata.var["SystemControl"] = adata.var_names.str.startswith("SystemControl")
    adata.var["Negative"] = adata.var_names.str.startswith("Negative")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["SystemControl", "Negative"], inplace=True, log1p=False)
    adata.obs["pct_counts_SystemControl"] = adata.obs["pct_counts_SystemControl"].fillna(0)
    adata.obs["pct_counts_Negative"] = adata.obs["pct_counts_Negative"].fillna(0)
    if (save_intermediate):
        adata.write(res_dir+ nm + "_adata_raw_qc.h5ad")
    
    logging.info("## QC Metrics calulated for the raw data")
    logging.info("## Chkpt 1")

    sc.pp.filter_cells(adata, inplace=True, min_counts=100)
    sc.pp.filter_genes(adata, inplace=True, min_cells=400)
    adata = adata[adata.obs['pct_counts_SystemControl'] < 5].copy()
    adata = adata[adata.obs['pct_counts_Negative'] < 5].copy()

    metrics_list = []
    for obs_value in adata.obs["tma_sid"].unique():
        metrics = calc_fov_qc_metrics(adata, "tma_sid", obs_value, (0.51,0.51))
        metrics_list.append(metrics)
    basic_metrics = pd.concat(metrics_list, ignore_index=True)
    basic_metrics.to_csv(res_dir+"basic_metrics_tma.csv", index=False)

    tma_sids_to_keep = basic_metrics[basic_metrics['Total Cells'] >= 225]['Sample']
    tma_sids_to_rem = basic_metrics[basic_metrics['Total Cells'] < 225]['Sample']
    print(tma_sids_to_rem)
    adata = adata[adata.obs['tma_sid'].isin(tma_sids_to_keep)].copy()

    neg_genes = [gene for gene in adata.var_names if gene.startswith('Negative')]
    ctlr_genes = [gene for gene in adata.var_names if gene.startswith('SystemControl')]
    adata = adata[:, ~adata.var_names.isin(ctlr_genes)]
    adata = adata[:, ~adata.var_names.isin(neg_genes)]
    if (save_intermediate):
        adata.write(res_dir + nm + "_adata_filtered.h5ad")

    logging.info("## Cells who have less than 100 counts are removed")
    logging.info("## Cells where the pct of system control counts and negative counts are > 5 are removed")
    logging.info("## Genes who are expressed in less than 400 cells are removed")
    logging.info("## Chkpt 2")

    adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(adata, target_sum=1e4, inplace=True)
    sc.pp.log1p(adata)
    if (save_intermediate):
        adata.write(res_dir + nm + "_adata_lognormalized.h5ad")
    logging.info("## Counts layer saved, counts total normalized to 10000 and log1p transformed")
    logging.info("## Chkpt 3")

    sc.pp.highly_variable_genes(adata, n_top_genes=6000, flavor="cell_ranger", subset=False)
    if (save_intermediate):
        adata.write(res_dir + nm + "_adata_hvg.h5ad")

    adata_hvg = adata[:, adata.var['highly_variable']]
    sc.pp.calculate_qc_metrics(adata_hvg, inplace=True, log1p=False)
    logging.info("## Top 5k highly variable genes are identified using cell_ranger")
    logging.info("## Chkpt 4")
    if (save_intermediate):
        adata_hvg.write(res_dir + nm + "_adata_hvg_subset.h5ad")
    logging.info("## Negative genes and system control genes removed")
    logging.info("## Chkpt 5")
    adata_hvg.raw = adata_hvg
    return adata_hvg

def batch_correction(adata, method, batch_factor, npcs = None, dis_met = "canberra"):
    #Note: UMAP uses the neighbours calculated by the sc.pp.neighbours..but tsne calculates it from scarch
    logging.info("### Batch Correction ###")
    logging.info("## Method = %s, Batch factor = %s", method, batch_factor)
    n_samples, n_features = adata.X.shape
    if (npcs is None):
        npcs = min(50, n_samples-1, n_features-1) 
    print(npcs, n_samples, n_features)

    if method == 'nobc':
        logging.info("### No Batch Correction NPC: npcs ###")
        adata_no_bc = adata.copy()
        #sc.pp.scale(adata_no_bc)
        sc.pp.pca(adata_no_bc, svd_solver = "arpack", n_comps=npcs)
        sc.pp.neighbors(adata_no_bc, n_pcs = npcs, metric=dis_met, random_state=123)
        sc.tl.umap(adata_no_bc, random_state=123)
        sc.tl.tsne(adata_no_bc, perplexity = npcs, metric=dis_met)
        return adata_no_bc
    
    elif method == 'bbknn_pynn':
        logging.info("### BBKNN - PYNN Variant NPC: npcs ###")
        adata_pynn = adata.copy()
        #sc.pp.scale(adata_pynn)
        sc.pp.pca(adata_pynn, svd_solver = "arpack", n_comps=npcs)
        sc.external.pp.bbknn(adata_pynn, batch_key=batch_factor, approx=True, use_annoy=False, n_pcs=npcs, metric=dis_met, random_state=123)
        sc.tl.umap(adata_pynn, random_state=123)
        sc.tl.tsne(adata_pynn, metric=dis_met)
        return adata_pynn
    
    elif method == 'bbknn_annoy':
        logging.info("### BBKNN - ANNOY Variant NPC: npcs ###")
        adata_annoy = adata.copy()
        #sc.pp.scale(adata_annoy)
        sc.pp.pca(adata_annoy, svd_solver = "arpack", n_comps=npcs)
        sc.external.pp.bbknn(adata_annoy, batch_key=batch_factor, approx=True, use_annoy=True, metric=dis_met, random_state=123)
        sc.tl.umap(adata_annoy, random_state=123)
        sc.tl.tsne(adata_annoy, metric=dis_met)
        return adata_annoy
    
    elif method == 'harmony':
        logging.info("### Harmony - NPC: npcs ###")
        adata_harm = adata.copy()
        #sc.pp.scale(adata_harm)
        sc.pp.pca(adata_harm, svd_solver = "arpack", n_comps=npcs)
        logging.info("### PCA Done ###")
        sc.external.pp.harmony_integrate(adata_harm, key=batch_factor, max_iter_harmony=20, plot_convergence=True)
        logging.info("### Harmony Done ###")
        sc.pp.neighbors(adata_harm, use_rep="X_pca_harmony", metric=dis_met, random_state=123)
        logging.info("### Neighbors Done ###")
        sc.tl.umap(adata_harm, random_state=123)
        logging.info("### UMAP Done ###")
        sc.tl.tsne(adata_harm, use_rep="X_pca_harmony", metric=dis_met)
        logging.info("### tSNE Done ###")
        return adata_harm

    elif method == 'combat':
        logging.info("### Combat - NPC: npcs ###")
        adata_combat = adata.copy()
        sc.pp.combat(adata_combat, key = batch_factor)
        #sc.pp.scale(adata_combat)
        sc.pp.pca(adata_combat, svd_solver = "arpack", n_comps=npcs)
        sc.pp.neighbors(adata_combat, metric=dis_met)
        sc.tl.umap(adata_combat)
        sc.tl.tsne(adata_annoy, metric=dis_met)
        return adata_combat

    elif method == 'scvi':
        logging.info("### SCVI - NPC: npcs ###")
        #No Scaling needed for SCVI
        import scvi
        adata_scvi = adata.copy()
        
        scvi.model.SCVI.setup_anndata(adata_scvi,layer="counts",batch_key=batch_factor)
        #model = scvi.model.SCVI(adata_scvi, n_hidden=128, n_latent=12, n_layers=1, dispersion='gene-batch')
        model = scvi.model.SCVI(adata_scvi)
        model.train()
        #model.save(res_dir+"scvi_model/")
        #model = scvi.model.SCVI.load(res_dir+"/scvi_model/", adata=adata_scvi)
        latent = model.get_latent_representation()
        adata_scvi.obsm["X_scVI"] = latent
        adata_scvi.layers["scvi_normalized"] = model.get_normalized_expression(library_size=10e4)
        sc.tl.pca(adata_scvi, n_comps = 50)
        sc.pp.neighbors(adata_scvi,use_rep="X_scVI", metric=dis_met, random_state=123)
        sc.tl.umap(adata_scvi, random_state=123)
        sc.tl.tsne(adata_scvi,use_rep="X_scVI", metric=dis_met)
        return adata_scvi
    else:
        logging.error(f"Unknown method: {method}")
        raise ValueError(f"Unknown method: {method}")

def cluster(adata, res):
    logging.info("### Clustering ###")
    logging.info("## leiden clustering, resolution = %s", res)

    sc.tl.leiden(adata, resolution=res, random_state=24) 
    return adata

def deg(adata, group_obs):
    logging.info("### DEG ###")
    logging.info("## res_dir = %s")
    logging.info("## wilcoxon test")
    # Need to use raw data for DEGs
    sc.tl.rank_genes_groups(adata, groupby=group_obs, method="wilcoxon", key_added="deg_"+group_obs)
    tcum = len(adata.obs[group_obs].unique())
    #cluster_names = adata.obs[group_obs].unique().tolist()
    cluster_names = list(adata.uns["deg_"+group_obs]['names'].dtype.names)
    #print(cluster_names)
    fdf = []
    for ii in range(tcum):
            print(ii)
            print(cluster_names[ii])
            df = pd.DataFrame({
            'gene': pd.DataFrame(adata.uns["deg_"+group_obs]['names']).iloc[:,ii],
            'logFC': pd.DataFrame(adata.uns["deg_"+group_obs]['logfoldchanges']).iloc[:,ii],
            'pvals_adj': pd.DataFrame(adata.uns["deg_"+group_obs]['pvals_adj']).iloc[:,ii],
            'Cluster': ii,
            'Name': cluster_names[ii]})
            fdf.append(df)
    fdf = pd.concat(fdf, ignore_index=True)
    return adata, fdf

def annotate_broad_clusters(deg_df, marker_list, pval_cutoff=10):
    df = deg_df.copy()
    df['is_marker'] = df['gene'].isin(marker_list['Markers'])
    df = df.sort_values(by=['logFC', 'Cluster'], ascending=[False, True])
    df['PredictedCellType'] = ""
    df['FinalPredictedCellType'] = ""  # NEW: final cell type per cluster repeated per row
    #print("Loading")
    annot_clust = {}

    for cnum in range(len(df['Cluster'].unique())):
        #print(cnum)
        score_here = pd.Series(0.0, index=marker_list['CellType'].unique())
        df_here = df[(df['is_marker']) & (df['Cluster'] == cnum) & (df['logFC'] >= 0) & (df['pvals_adj'] <= pval_cutoff)]
        #print(df_here)

        if df_here.empty:
            #df_here = df[(df['is_marker']) & (df['Cluster'] == cnum) & (df['logFC'] >= 0)]
            annot_clust[cnum] = "Unknown"
        else:
            for i in range(len(df_here)):
                row = df_here.iloc[i]  # Get the row using iloc
                gn = row['gene']
                ml = marker_list[marker_list['Markers'].isin([gn])]

                if not ml.empty:
                    cell_types = ml['CellType'].unique()
                    for cell_type in cell_types:
                        score_here[cell_type] += row['logFC']
                        df.loc[(df['Cluster'] == cnum) & (df['gene'] == gn), 'PredictedCellType'] = ', '.join(cell_types)

                        #cell_type = ml['CellType'].values[0]  # Get the first matching cell type
                        #score_here[cell_type] += row['logFC']
                        #df.loc[(df['Cluster'] == cnum) & (df['gene'] == gn), 'PredictedCellType'] = cell_type

            max_score = score_here.max()
            final_ct = score_here[score_here == max_score].index.values[0]
            df.loc[df['Cluster'] == cnum, 'FinalPredictedCellType'] = final_ct
            #print(cnum)
            #print(score_here)
            annot_clust[cnum] = score_here[score_here == max_score].index.values[0]

    return {'data': df, 'annot': annot_clust}

def hier_cluster_marker_annotate(adata, deg_df, markerlist, newobs, pvalcutoff=10):
    logging.info("### Hierachical Annotation ###")
    logging.info("### Broad Annotations based on the leiden clustering but the T-Cell annotation based on markers ###")
    #markerlist = markerlist[~markerlist['CellType'].isin(['T Cells'])]
    annres = annotate_broad_clusters(deg_df, markerlist, pvalcutoff)
    adata_ann = adata.copy()
    adata_ann.obs[newobs] = adata_ann.obs['leiden']
    adata_ann.obs[newobs] = adata_ann.obs[newobs].astype(int).map(annres['annot'])

    if "TumorType" in adata_ann.obs:
        adata_ann.obs.loc[((adata_ann.obs["TumorType"] == "Post-NIVO Normal") &
                        (adata_ann.obs[newobs] == "Tumor")), newobs] = "Epithelial Cells"
    
    Tmarkers = ['CD3D', 'CD3E', 'CD3G', 'CD247']
    #Tmarkers = ['CD3']
    #Tmarkers = ['CD3D', 'CD3E', 'CD3G', 'CD247', 'TRAC', "TRBC1", "TRBC2"]
    Tmask = adata_ann.X[:,adata_ann.var_names.isin(Tmarkers)].toarray() > 0.1
    Tmask = adata_ann.X[:,adata_ann.var_names.isin(Tmarkers)] > 0.1
    Tmask_dense = Tmask.toarray()  # convert sparse → dense
    #adata_ann.obs.loc[Tmask_dense.any(axis=1), newobs] = "T cells" no need for fibro
    #adata_ann.obs.loc[Tmask.any(axis=1), newobs] = "T Cells"
    return adata_ann, annres

def annotate_tcells_subtypes(deg_df, marker_list, pval_cutoff=0.05):
    df = deg_df.copy()
    df['is_marker'] = df['gene'].isin(marker_list['Markers'])
    df = df.sort_values(by=['logFC', 'Cluster'], ascending=[False, True])
    df = df[df['logFC'] > 0]
    df['PredictedCellType'] = ""
    df['MarkerGenesInCluster'] = ""  # list of marker genes per cluster repeated per row
    df['FinalPredictedCellType'] = ""  # NEW: final cell type per cluster repeated per row

    print("Loading")
    annot_clust = {}
    scores_list = []

    for cnum in range(len(df['Cluster'].unique())):
        #print(f"Processing cluster {cnum}")
        score_here = pd.Series(0.0, index=marker_list['CellType'].unique())
        df_here = df[(df['is_marker']) & (df['Cluster'] == cnum)]

        if df_here.empty:
            annot_clust[cnum] = "Unknown"
            for ct in marker_list['CellType'].unique():
                scores_list.append({'cluster': cnum, 'celltype': ct, 'score': 0.0})
            df.loc[df['Cluster'] == cnum, 'MarkerGenesInCluster'] = ""
            df.loc[df['Cluster'] == cnum, 'FinalPredictedCellType'] = "Unknown"
        else:
            df_tmp = df_here[(df_here['logFC']) > 0].copy()
            df_here = df_tmp[(df_tmp['pvals_adj']) <= pval_cutoff].copy()
            marker_genes = df_here['gene'].unique()
            marker_genes_str = ', '.join(marker_genes)
            df.loc[df['Cluster'] == cnum, 'MarkerGenesInCluster'] = marker_genes_str
            if df_here.empty:
                annot_clust[cnum] = "Unknown"
                for ct in marker_list['CellType'].unique():
                    scores_list.append({'cluster': cnum, 'celltype': ct, 'score': 0.0})
                df.loc[df['Cluster'] == cnum, 'MarkerGenesInCluster'] = ""
                df.loc[df['Cluster'] == cnum, 'FinalPredictedCellType'] = "Unknown"
                #df_here = df_tmp
                #marker_genes = df_here['gene'].unique()
                #df.loc[df['Cluster'] == cnum, 'MarkerGenesInCluster'] = marker_genes_str
            for i in range(len(df_here)):
                row = df_here.iloc[i]
                gn = row['gene']
                ml = marker_list[marker_list['Markers'].isin([gn])]
                if not ml.empty:
                    cell_types = ml['CellType'].unique()
                    for _, subrow in ml.iterrows():
                        score_here[subrow['CellType']] += row['logFC']

                    df.loc[(df['Cluster'] == cnum) & (df['gene'] == gn), 'PredictedCellType'] = ', '.join(cell_types)

            #print(score_here)
            max_score = score_here.max()
            final_ct = score_here[score_here == max_score].index.values[0]
            annot_clust[cnum] = final_ct

            # Assign final predicted cell type to all rows of cluster
            df.loc[df['Cluster'] == cnum, 'FinalPredictedCellType'] = final_ct

            for ct in score_here.index:
                scores_list.append({'cluster': cnum, 'celltype': ct, 'score': score_here[ct]})

    scores_df = pd.DataFrame(scores_list)

    return df, annot_clust, scores_df

def score_tcell_clusters(adata, marker_list, cluster_key='leiden', scale='zscore'):
    from scipy.stats import zscore

    expr_matrix = adata.X.copy()
    expr_df = pd.DataFrame(expr_matrix.toarray(), index=adata.obs_names, columns=adata.var_names)

    markers_in_data = marker_list['Markers'][marker_list['Markers'].isin(expr_df.columns)].unique()
    expr_df = expr_df[markers_in_data]
    #expr_df = expr_df.apply(lambda x: zscore(x, ddof=0), axis=0).fillna(0.0)
    clusters = adata.obs[cluster_key].unique()
    sigs = marker_list['Signatures'].unique()
    scores = []
    for cl in clusters:
        cells_in_cluster = adata.obs_names[adata.obs[cluster_key] == cl]
        cluster_expr = expr_df.loc[cells_in_cluster]
        for ss in sigs:
            markers_ss = marker_list[marker_list['Signatures'] == ss]['Markers'].values
            markers_ss = [g for g in markers_ss if g in cluster_expr.columns]
            if len(markers_ss) == 0:
                score = 0.0
            else:
                score = cluster_expr[markers_ss].mean().mean()
            scores.append({'cluster': cl, 'signatures': ss, 'score': score})
    scores_df = pd.DataFrame(scores)

    scores_df['scaled_score'] = scores_df.groupby('cluster')['score'].transform(lambda x: zscore(x, ddof=0))
    #scores_df['scaled_score'] = scores_df.groupby('cluster')['score'].transform(lambda x: x)
    scores_matrix = scores_df.pivot(index='signatures', columns='cluster', values='scaled_score')

    return scores_matrix

def cluster_subtype_annotate(adata_ann, deg_df, markerlist, newobs):
    logging.info("### Annotation ###")
    logging.info("## marker list provided#")
    annres = annotate_tcells_subtypes(deg_df, markerlist)
    adata_ann.obs[newobs] = adata_ann.obs['leiden']
    adata_ann.obs[newobs] = adata_ann.obs[newobs].astype(int).map(annres['annot'])
    return adata_ann, annres

def deg_dotplot(adata, obs):
    sc.tl.rank_genes_groups(adata, groupby=obs, method="wilcoxon", key_added= "deg_" + obs)
    sc.tl.dendrogram(adata, groupby=obs)
    sc.pl.rank_genes_groups_dotplot(adata, n_genes = 5, groupby=obs, key="deg_" + obs,
        values_to_plot="logfoldchanges", cmap='bwr', colorbar_title='Log Fold Change',
        title='Log FC of top 5 genes', save = "deg_dotplot.png", show = False)
    
def deg_dotplot_summary(adata, obs, n_top, top_genes = None):
    print(top_genes)
    adata, fdf = deg(adata, obs)
    fdf_top = (
    fdf.sort_values(["Cluster", "logFC",  "pvals_adj"], ascending=[True, False, True])
    .groupby("Cluster", group_keys=False)
    .head(n_top))
    if (top_genes == None):
        print("Here")
        top_genes = fdf_top["gene"].unique().tolist()
    clusters = sorted(adata.obs[obs].unique().tolist())

    print(f"Selected {len(top_genes)} unique top genes across {len(clusters)} clusters.")
    source_map = dict(zip(fdf_top["gene"], fdf_top["Name"]))
    X = adata[:, top_genes].X
    if scipy.sparse.issparse(X):
        X = X.toarray()

    # --- Step 3: Compute per-cluster average and % expression ---
    summary_list = []

    for cl in clusters:
        mask = adata.obs[obs] == cl
        X_cl = X[mask, :]
        X_rest = X[~mask, :]  # all other cells

        avg_expr = np.asarray(X_cl.mean(axis=0)).flatten()
        pct_expr = np.asarray((X_cl > 0).sum(axis=0) / X_cl.shape[0] * 100).flatten()
        logFC = np.log2((avg_expr + 1e-9) / (np.asarray(X_rest.mean(axis=0)).flatten() + 1e-9))


        df = pd.DataFrame({
            "cluster": cl,
            "gene": top_genes,
            "avg_expr": avg_expr,
            "pct_expr": pct_expr,
            "logFC": logFC,
            "source_cluster": [source_map.get(g, None) for g in top_genes]
        })
        summary_list.append(df)

    expr_summary = pd.concat(summary_list, ignore_index=True)
    return fdf_top, expr_summary

def select_top_genes(fdf, top_n=15, pval_cutoff=0.05):
    # Keep only significant genes
    sig = fdf.query(f"pvals_adj <= {pval_cutoff}").copy()
    
    # Compute score
    sig["score"] = sig["logFC"].abs() * -np.log10(sig["pvals_adj"])
    
    # Collapse multiple entries per gene: take largest absolute logFC
    sig_gene = sig.groupby("gene").agg(
        logFC=("logFC", lambda x: x.loc[x.abs().idxmax()]),
        score=("score", "max")
    ).reset_index()
    # Top upregulated genes
    top_up = (
        sig_gene.query("logFC > 0")
                .sort_values("score", ascending=False)
                .head(top_n)
                .gene
                .tolist()
    )
    # Top downregulated genes
    top_down = (
        sig_gene.query("logFC < 0")
                .sort_values("score", ascending=False)
                .head(top_n)
                .gene
                .tolist()
    )
    top_genes = top_up + top_down
    return top_genes


def make_heatmap_matrix(adata, gene_list, group_cols, meta_cols, agg_func='mean'):
    gene_mask = adata.var_names.isin(gene_list)
    gene_order = adata.var_names[gene_mask]  # preserves order in adata.var_names
    X = adata.X[:, gene_mask]
    df = pd.DataFrame(X.toarray(), columns=gene_order, index=adata.obs_names)
    for col in group_cols + meta_cols:
        df[col] = adata.obs[col].values
    df_agg = df.groupby(group_cols)[gene_list].agg(agg_func).reset_index()
    meta_map = adata.obs[group_cols + meta_cols].drop_duplicates()
    df_agg = df_agg.merge(meta_map, on=group_cols)
    df_agg['sample'] = df_agg[group_cols + meta_cols].astype(str).agg('_'.join, axis=1)
    df_agg = df_agg[['sample'] + group_cols + meta_cols + gene_list]
    return df_agg

from scipy.spatial import distance_matrix
def calc_distance_matrix_per(adata_sam):
    coords = adata_sam.obsm['spatial_fov']
    obs_df = adata_sam.obs.copy()
    obs_df["centroid_x"] = coords[:, 0]
    obs_df["centroid_y"] = coords[:, 1]
    cells = obs_df['cell'].to_numpy()
    centroids = obs_df[["centroid_x", "centroid_y"]].to_numpy()
    dist_mat = distance_matrix(centroids, centroids) #Minkowski with p = 2 same as euclidean
    dist_df = pd.DataFrame(dist_mat, index=cells, columns=cells)
    return dist_df

def calc_distance_matrix_per(adata_sam):
    # Check if 'spatial_fov' exists in obsm
    if 'spatial_fov' in adata_sam.obsm:
        coords = adata_sam.obsm['spatial_fov']
        # If it's a DataFrame, convert to array
        if isinstance(coords, pd.DataFrame):
            coords = coords.to_numpy()
        centroid_x = coords[:, 0]
        centroid_y = coords[:, 1]
    else:
        # fallback to obs
        centroid_x = adata_sam.obs['centroid_x'].to_numpy()
        centroid_y = adata_sam.obs['centroid_y'].to_numpy()
    
    # Build obs dataframe
    obs_df = adata_sam.obs.copy()
    obs_df['centroid_x'] = centroid_x
    obs_df['centroid_y'] = centroid_y
    
    # Cells and centroids
    cells = obs_df['cell'].to_numpy()
    centroids = obs_df[['centroid_x', 'centroid_y']].to_numpy()
    
    # Distance matrix
    dist_mat = distance_matrix(centroids, centroids)
    dist_df = pd.DataFrame(dist_mat, index=cells, columns=cells)
    
    return dist_df

def calc_query_dist_target_per(adata_sam, obs, query_ct, target_ct):
    #logging.info("### Calculating Distance of query cell to the target cell ###")
    #logging.info("### Calculation is happening per obs..Be careful ###")
    #logging.info("## Obs = %s, Query Cell Type = %s, Target Cell Type = %s", obs, query_ct, target_ct)

    adjMtx = calc_distance_matrix_per(adata_sam)
    target_ct_adata = adata_sam[adata_sam.obs[obs] == target_ct].copy()
    target_ct_cells = target_ct_adata.obs['cell'].to_numpy()

    query_ct_adata = adata_sam[adata_sam.obs[obs] == query_ct].copy()
    query_ct_cells = query_ct_adata.obs['cell'].to_numpy()

    dis = []
    for query_cell_i in query_ct_cells:
        distances = adjMtx.loc[query_cell_i, target_ct_cells]
        if distances.sum() == 0:
            max_distance = 25000
            dis.append(max_distance)
        else:
            min_distance = distances[distances > 0].min() 
            dis.append(min_distance)  
    query_ct_adata.obs['Distance_to_' + target_ct] = dis
    return query_ct_adata

def get_query_dist_target_all(adata, obs, query_ct, target_ct):
    sam_obs_list = adata.obs['tma_sid'].unique()
    results = []
    for i in sam_obs_list:
        logging.info("## SAMID_FOV: %s", i)
        sam1 = adata[adata.obs['tma_sid'] == i].copy()

        #num_fovs = sam1.obs['fov'].nunique()
        num_query = (sam1.obs[obs] == query_ct).sum()
        num_target  = (sam1.obs[obs] == target_ct).sum()
        
        print(f"\n🧬 Sample: {i}")
        #print("-" * 40)
        #print(f"Unique FOVs      : {num_fovs}")
        print(f"Query count     : {num_query}")
        print(f"Target count : {num_target}")
        print("=" * 40)

        classified_cells = calc_query_dist_target_per(sam1, obs, query_ct, target_ct)
        results.append(classified_cells)

    query_with_dist = ad.concat(results)
    return (query_with_dist)

def get_query_tn_target_all(adata, obs, query_ct, target_ct):
    df_all = pd.read_csv("/Users/lakshmi_nccs/Desktop/NCCS/Projects/Cosmx6k/processed_data/tma_all_touching_neighbors.csv")

    query_ct_adata = adata[adata.obs[obs] == query_ct].copy()
    query_ct_cells = query_ct_adata.obs['cell'].to_numpy()
    edges_query = df_all[df_all["cell"].isin(query_ct_cells)].copy()

    target_ct_adata = adata[adata.obs[obs] == target_ct].copy()
    target_ct_cells = target_ct_adata.obs['cell'].to_numpy()
    edges_query_target = edges_query[edges_query['neighbor'].isin(target_ct_cells)]

    query_cells_touching_target = np.unique(edges_query_target['cell'].to_numpy())
    query_ct_adata.obs['TN_status:' + target_ct] = "Not TN of " + target_ct
    query_ct_adata.obs.loc[
        query_ct_adata.obs['cell'].isin(query_cells_touching_target),
        'TN_status:' + target_ct
    ] = "TN of " + target_ct
    return (query_ct_adata)

def get_spat_fingerprint_query(adata_sam, query_obs, query_voi, res_obs, rad):
    logging.info("### Get Neighbours per FOV ###")
    logging.info("## Query Obs = %s, Radius = %s", query_obs, rad)
    
    if (rad == "TN"):
        df_all = pd.read_csv("/Users/lakshmi_nccs/Desktop/NCCS/Projects/Cosmx6k/processed_data/tma_all_touching_neighbors.csv")
        valid_cells = adata_sam.obs['cell'].values
        df_all = df_all[df_all['cell'].isin(valid_cells) & df_all['neighbor'].isin(valid_cells)].copy()

        query_ct_adata = adata_sam[adata_sam.obs[query_obs] == query_voi].copy()
        query_ct_cells = query_ct_adata.obs['cell'].to_numpy()
        neighbor_ci_list = []
        for ci in query_ct_cells:
            edges_query = df_all[df_all["cell"] == ci].copy()
            neighb_cell = np.unique(edges_query[['neighbor']].to_numpy())
            neighb_adata = adata_sam[adata_sam.obs['cell'].isin(neighb_cell)].copy()
            neighbor_ct_counts = neighb_adata.obs[res_obs].value_counts()
            row = adata_sam.obs.loc[adata_sam.obs['cell'] == ci].iloc[0]
            counts_dict = neighbor_ct_counts.to_dict()
            counts_dict.update({
                "cell": ci,
                "tma_sid": row['tma_sid'],
                "tt": row['TumorType'],
                "patient": row['Patient'],
                "tcelltype": row['tcell_type_resolved'],
                "query_obs": row[query_obs],
                "res_obs": row[res_obs],
            })
            neighbor_ci_list.append(counts_dict)
        neighbor_proportions_df = pd.DataFrame(neighbor_ci_list).fillna(0)
        neighbor_proportions_df.index = query_ct_cells
        return neighbor_proportions_df
    else:
        adjMtx = calc_distance_matrix_per(adata_sam)
        query_ct_cells = adata_sam[adata_sam.obs[query_obs] == query_voi].copy()
        ct_ci = query_ct_cells.obs['cell'].to_numpy()
        #ct_ci = adata_sam.obs.index
        neighbor_ci_list = []
        for ci in ct_ci:
            #neighbors_ci = adjMtx.index[adjMtx[ci] > 0]
            neighbors_ci = adjMtx.index[(adjMtx[ci] > 0) & (adjMtx[ci] <= rad)]
            neighbors_ci = [c for c in neighbors_ci if c in adata_sam.obs['cell'].values]
            neighb_adata = adata_sam[adata_sam.obs['cell'].isin(neighbors_ci)].copy()
            neighbor_ct_counts = neighb_adata.obs[res_obs].value_counts()
            row = adata_sam.obs.loc[adata_sam.obs['cell'] == ci].iloc[0]
            counts_dict = neighbor_ct_counts.to_dict()
            counts_dict.update({
                "cell": ci,
                "tma_sid": row['tma_sid'],
                "tt": row['TumorType'],
                "patient": row['Patient'],
                "tcelltype": row['tcell_type_resolved'],
                "query_obs": row[query_obs],
                "res_obs": row[res_obs],
            })
            neighbor_ci_list.append(counts_dict)
        neighbor_proportions_df = pd.DataFrame(neighbor_ci_list).fillna(0)
        neighbor_proportions_df.index = ct_ci
        return neighbor_proportions_df

def make_fingerprint_adata (full_df, obs_cols):
    obs_df = full_df[obs_cols]
    x_df = full_df.drop(columns=obs_cols)
    x_df.index = obs_df['cell']
    fing_adata = ad.AnnData(X=x_df.values, obs=obs_df)
    fing_adata.var['feature_names'] = x_df.columns
    fing_adata.var_names = x_df.columns
    fing_adata.obs_names = x_df.index   
    return fing_adata

def cluster_sf (sf_adata, npca, nn, res, filter = True):
    if (filter == True):
        print(sf_adata.var_names)
        keep_genes = (sf_adata.X >= 20).sum(axis=0) >= 200
        print(sf_adata.var_names[keep_genes])
        sf_adata = sf_adata[:, keep_genes].copy()

    sf_adata.layers["counts"] = sf_adata.X.copy()
    sc.pp.normalize_total(sf_adata, target_sum=100, inplace=True)
    #sc.pp.log1p(sf_adata)
    #sc.pp.scale(sf_adata)
    sc.tl.pca(sf_adata, n_comps=npca)
    #sc.external.pp.harmony_integrate(sf_adata, key="patient", max_iter_harmony=20, plot_convergence=True)
    #sc.pp.neighbors(sf_adata, n_pcs = 10, use_rep="X_pca_harmony", metric="euclidean")
    sc.pp.neighbors(sf_adata, n_pcs = nn, use_rep='X_pca', metric = "canberra")
    sc.tl.umap(sf_adata)
    sc.tl.tsne(sf_adata, metric="canberra" )
    sc.tl.leiden(sf_adata, resolution=res) 
    from sklearn.cluster import KMeans
    kmeans = KMeans(n_clusters=5, random_state=0).fit(sf_adata.X)
    sf_adata.obs['kmeans'] = kmeans.labels_.astype(str)
    return(sf_adata)


def assign_dominant_contact_niche(adata):
    df = pd.DataFrame(adata.X, index=adata.obs_names, columns=adata.var_names)
    df = df.fillna(0)
    df_max = df.max(axis=1)
    dominant = df.idxmax(axis=1)
    dominant[df_max == 0] = 'No-contact'
    adata.obs['contact_niche'] = pd.Categorical(dominant.where(dominant == "No-contact", dominant + "_rich_niche"))
    return adata


def assign_niches_from_leiden(adata, pure_threshold=90, dominant_thresh=60, max_top_neighbors=3):
    """
    Assigns niche labels to Leiden clusters based on mean neighbour-type fractions.
    Assumes adata.X are neighbour-type counts or frequencies.
    """
    # Normalize neighbour-type fractions per cell
    X = adata.X.copy()
    X = np.nan_to_num(X, nan=0.0)
    df = pd.DataFrame(X, index=adata.obs_names, columns=adata.var_names)
    leiden_means = df.groupby(adata.obs['leiden']).mean()
    
    # Visualize
    plt.figure(figsize=(12, 6))
    sns.heatmap(leiden_means, cmap="viridis", annot=True, fmt=".2f")
    plt.title("Mean Neighbour Fractions per Leiden Cluster")
    plt.ylabel("Leiden Cluster")
    plt.xlabel("Neighbour Type")
    plt.show()
    
    # Assign niche names based on dominant fractions
    niche_labels = {}
    for cluster, row in leiden_means.iterrows():
        sorted_types = row.sort_values(ascending=False)
        top_neighbors = sorted_types.index[:max_top_neighbors]
        top_neighbors_clean = [n.split()[0] for n in top_neighbors]
        top_fracs = sorted_types.values[:max_top_neighbors]
        top1, top2 = top_fracs[0], top_fracs[1] if len(top_fracs) > 1 else 0.0
        print(top_neighbors_clean)
        # Rule 1: pure niche (strongly dominated by one type)
        print(top1)
        if top1 >= pure_threshold:
            niche = f"{top_neighbors_clean[0]}_rich_niche"
        # Rule 2: clear dominant even if not pure
        elif (top1 - top2) >= 0.75 * top1 or top1 >= dominant_thresh:
            niche = f"{top_neighbors_clean[0]}_rich_niche"
        # Rule 3: mixed niche with top contributors
        else:
            mixed = "_".join(top_neighbors_clean)
            niche = f"{mixed}_mixed"
        
        print(f"Leiden {cluster} → {niche} (top fracs: {top_fracs})")
        niche_labels[cluster] = niche
    
    # Assign niche label to each cell based on its Leiden cluster
    adata.obs["spatial_niche"] = adata.obs["leiden"].map(niche_labels).astype("category")
    return adata

from liana.method import cellphonedb
import numpy as np

def get_LR_neibs_adata(adata, lig, rec):
    if lig.startswith("Q"):
        lig_q = "T Cells"
        lig_mask = (adata.obs['ct'] == "T Cells") & (adata.obs["spatial_bin"] == lig)
    else:
        lig_q = lig
        lig_mask = adata.obs['ct'] == lig

    if rec.startswith("Q"):
        rec_q = "T Cells"
        rec_mask = (adata.obs['ct'] == "T Cells") & (adata.obs["spatial_bin"] == rec)
    else:
        rec_q = rec
        rec_mask = adata.obs['ct'] == rec
    
    adata_sub = ad.concat([adata[lig_mask].copy(),adata[rec_mask].copy()])
    
    neibs_of_rec = get_query_tn_target_all(adata_sub, 'ct', lig_q, rec_q)
    print(neibs_of_rec.obs['TN_status:'+rec_q].value_counts())
    neibs_of_rec = neibs_of_rec[neibs_of_rec.obs['TN_status:'+rec_q] == "TN of "+rec_q].copy()
    

    neibs_of_lig = get_query_tn_target_all(adata_sub, 'ct', rec_q, lig_q)
    print(neibs_of_lig.obs['TN_status:'+lig_q].value_counts())
    neibs_of_lig = neibs_of_lig[neibs_of_lig.obs['TN_status:'+lig_q] == "TN of "+lig_q].copy()

    adata_neibs = ad.concat([neibs_of_lig,neibs_of_rec])
    return(adata_neibs)

def get_LR_df(adata_neibs, tcell_bin):
    adata_resp = adata_neibs[adata_neibs.obs['Response'] == 'Responder'].copy()
    adata_nonresp = adata_neibs[adata_neibs.obs['Response'] == 'Non_Responder'].copy()

    res_cpdb_resp = cellphonedb(adata_resp,groupby="ct",resource_name="consensus",expr_prop=0.1,key_added="cpdb_resp")
    res_cpdb_nonresp = cellphonedb(adata_nonresp,groupby="ct",resource_name="consensus",expr_prop=0.1, key_added="cpdb_nonresp")

    interactions_resp = adata_resp.uns['cpdb_resp'].copy()
    interactions_nonresp = adata_nonresp.uns['cpdb_nonresp'].copy()

    pval_cutoff = 0.05
    prop_cutoff = 0.2
    interactions_resp_filt = interactions_resp[
        (interactions_resp["cellphone_pvals"] < pval_cutoff) &
        (interactions_resp["ligand_props"] > prop_cutoff) &
        (interactions_resp["receptor_props"] > prop_cutoff)
    ].copy()

    interactions_nonresp_filt = interactions_nonresp[
        (interactions_nonresp["cellphone_pvals"] < pval_cutoff) &
        (interactions_nonresp["ligand_props"] > prop_cutoff) &
        (interactions_nonresp["receptor_props"] > prop_cutoff)
    ].copy()

    merged_lr = interactions_resp_filt.merge(
        interactions_nonresp_filt,
        on=["ligand", "receptor", "source", "target"],
        suffixes=("_Resp", "_NonResp"),
        how="outer"
    )

    merged_lr["lr_means_Resp"] = merged_lr["lr_means_Resp"].fillna(0)
    merged_lr["lr_means_NonResp"] = merged_lr["lr_means_NonResp"].fillna(0)

    merged_lr["cellphone_pvals_Resp"] = merged_lr["cellphone_pvals_Resp"].fillna(1)
    merged_lr["cellphone_pvals_NonResp"] = merged_lr["cellphone_pvals_NonResp"].fillna(1)

    merged_lr = merged_lr[
        (merged_lr["cellphone_pvals_Resp"] < pval_cutoff) |
        (merged_lr["cellphone_pvals_NonResp"] < pval_cutoff)
    ].copy()

    merged_lr["score_diff"] = (
        merged_lr["lr_means_Resp"] - merged_lr["lr_means_NonResp"]
    )

    merged_lr["log2_fc"] = np.log2(
        (merged_lr["lr_means_Resp"] + 1e-3) /
        (merged_lr["lr_means_NonResp"] + 1e-3)
    )

    merged_lr = merged_lr.sort_values(
        by="log2_fc",
        key=lambda x: np.abs(x),
        ascending=False
    )

    merged_lr["interaction_type"] = "Shared"

    merged_lr.loc[
        (merged_lr["cellphone_pvals_Resp"] < 0.05) &
        (merged_lr["cellphone_pvals_NonResp"] >= 0.05),
        "interaction_type"
    ] = "Responder-specific"

    merged_lr.loc[
        (merged_lr["cellphone_pvals_NonResp"] < 0.05) &
        (merged_lr["cellphone_pvals_Resp"] >= 0.05),
        "interaction_type"
    ] = "Non_Responder-specific"

    merged_lr = merged_lr[merged_lr['source'] != merged_lr['target']].copy()
    # Example: add Q-bin info to source and target
    merged_lr['source'] = np.where(
        merged_lr['source'] == "T Cells",
        merged_lr['source'] + f" ({tcell_bin})",
        merged_lr['source']
    )
    merged_lr['target'] = np.where(
        merged_lr['target'] == "T Cells",
        merged_lr['target'] + f" ({tcell_bin})",
        merged_lr['target']
    )

    merged_lr['pair_genes'] = merged_lr['ligand'] + '-' + merged_lr['receptor']
    merged_lr['pair_ct'] = merged_lr['source'] + '-' + merged_lr['target']
    merged_lr["pair_ct_genes"] = merged_lr['pair_ct'] + " | " + merged_lr['pair_genes']
    return(merged_lr)





def analyze_proportions(adata, split_obs, ct_obs, test_stat):
    # Step 1: Prepare the data
    df = adata.obs.copy()
    proportions = df.groupby(['TumorType', split_obs, ct_obs, 'Response']).size().reset_index(name='count')
    total_counts = proportions.groupby([split_obs, 'TumorType'])['count'].sum().reset_index(name='total_count')
    proportions = proportions.merge(total_counts, on=[split_obs, 'TumorType'])
    proportions['proportion'] = proportions['count'] / proportions['total_count']

    stat_results = []

    # Unique cell types
    cell_types = proportions[ct_obs].unique()

    # Step 4: Perform the paired Wilcoxon test for each cell type
    for cell_type in cell_types:
        # Filter for the current cell type
        cell_data = proportions[proportions[ct_obs] == cell_type]
        pivoted_data = cell_data.pivot(index=split_obs, columns='TumorType', values='proportion').reset_index()
        
        pre_values = pivoted_data['Pre-NIVO Primary']
        post_values = pivoted_data['Post-NIVO Primary']
        
        # Perform the Wilcoxon test
        if test_stat == 'wilcoxon':
            from scipy.stats import wilcoxon
            stat, p_value = wilcoxon(pre_values, post_values, method = 'exact')
        elif test_stat == 'ttest':
            from scipy.stats import ttest_rel
            # Perform the paired t-test 
            stat, p_value = ttest_rel(pre_values, post_values)
        elif test_stat == 'sign':
            from statsmodels.stats.descriptivestats import sign_test

            diffs = post_values - pre_values
            stat, p_value = sign_test(diffs)
            # Perform the sign test
            #stat, p_value = sign_test(pre_values, post_values)
        elif test_stat == 'permutation':
            from scipy.stats import permutation_test
            import numpy as np

            def statistic(x, y):
                return np.median(x - y)

            res = permutation_test((post_values, pre_values), statistic, vectorized=False,
                                permutation_type='pairings', alternative='two-sided', n_resamples=10000)
            p_value = res.pvalue
            stat = res.statistic
            print(f"Permutation test p-value: {res.pvalue:.4f}")
        else:
            raise ValueError(f"Unknown test statistic: {test_stat}")
        #ans = ans[ans['ct_1'].isin(['T Cells'])].copy()
        # Map treatment to numeric
        #ans['tt_num'] = ans['TumorType'].map({'Pre-NIVO Primary': 1, 'Post-NIVO Primary': 0})
        #model = smf.mixedlm("proportion ~ tt_num", ans, groups=ans["Patient"]).fit()
        #print(model.summary())
        #stat, p_value = ttest_rel(pre_values, post_values)
        
        # Append results to list
        stat_results.append({
            ct_obs: cell_type,
            'statistic': stat,
            'p_value': p_value
        })

    # Step 5: Create a DataFrame from the Wilcoxon results
    results_df = pd.DataFrame(stat_results)
    final_proportions = proportions.merge(results_df, on=ct_obs, how='left')
    #print(pivoted_data.head())
    #pivoted_data = pivoted_data.merge(results_df, on=ct_obs, how='left')
    # Save the results to a CSV file
    #final_proportions.to_csv(res_dir + split_obs + "_" + ct_obs + "_proportions_stat_results.csv", index=False)
    
    return final_proportions

def analyze_proportions(adata, obs1, obs2, ct_obs, resp_df):
    # Step 1: Prepare the data
    df = adata.obs.copy()
    proportions = df.groupby([obs1, obs2, ct_obs]).size().reset_index(name='count')
    total_counts = proportions.groupby([obs1, obs2])['count'].sum().reset_index(name='total_count')
    proportions = proportions.merge(total_counts, on=[obs1, obs2])
    proportions['proportion'] = proportions['count'] / proportions['total_count']
    proportions = proportions.merge(resp_df, on=[obs1], how='left')
    return proportions

def run_degs_by_celltype_and_response(adata,ct_obs="ct_1",resp_obs="Response",tt_obs="TumorType"):
    results = []
    tumor_groups_to_compare = adata.obs[tt_obs].unique().tolist()
    for ct in adata.obs[ct_obs].unique():
        for resp in adata.obs[resp_obs].unique():
            print(f"🔍 Processing: Cell type = {ct}, Response = {resp}")
            # Subset
            subset = adata[(adata.obs[ct_obs] == ct) &(adata.obs[resp_obs] == resp)].copy()
            subset = subset[subset.obs[tt_obs].isin(tumor_groups_to_compare)].copy()
            if len(subset.obs[tt_obs].unique()) < 2:
                print(f"⚠️ Skipping {ct} - {resp}: less than 2 tumor types.")
                continue
            
            sc.tl.rank_genes_groups(subset,groupby=tt_obs,method='wilcoxon',key_added="rank_genes")
            for gr in tumor_groups_to_compare:
                df = sc.get.rank_genes_groups_df(subset,group=gr,key="rank_genes")
                df['resp_tt'] = resp+ "::" + gr
                df[ct_obs] = ct
                df[resp_obs] = resp
                df[tt_obs] = gr  
                results.append(df)

    all_deg = pd.concat(results, ignore_index=True)
    return all_deg

def agg_by_patient_pseudo_bulk(adata, obs1, obs2, obs3, obs4):
    adata.obs['group'] = (
        adata.obs[obs1].astype(str) + '::' +
        adata.obs[obs2].astype(str) + '::' +
        adata.obs[obs3].astype(str) + '::' +
        adata.obs[obs4].astype(str)
    )

    expr_df = pd.DataFrame(adata.X.toarray(), index=adata.obs.index, columns=adata.var_names)
    expr_df['group'] = adata.obs['group'].values
    agg_expr = expr_df.groupby('group').mean()

    group_meta = agg_expr.index.to_series().str.split('::', expand=True)
    group_meta.columns = [obs1, obs2, obs3, obs4]
    group_meta = group_meta.reset_index().rename(columns={'index': 'group'})
    agg_expr = agg_expr.T
    return agg_expr,group_meta

def gene_distance_correlation(adata, genes=None, sortby = "Spearman_rho", distance_col="Distance_to_Tumor"):
    import pandas as pd
    from scipy.stats import spearmanr, pearsonr

    """
    Calculate Spearman and Pearson correlation between gene expression and distance to tumor.

    Parameters
    ----------
    adata : AnnData
        Single-cell AnnData object
    genes : list, optional
        List of genes to calculate correlations for. If None, use all genes in adata.var_names
    distance_col : str
        Column in adata.obs that contains distance to tumor

    Returns
    -------
    pd.DataFrame
        DataFrame with columns: Gene, Spearman_rho, Spearman_pval, Pearson_r, Pearson_pval
    """
    
    if genes is None:
        genes = adata.var_names.tolist()
    
    # Extract expression matrix
    expr = pd.DataFrame(
        adata[:, genes].X.toarray(),
        index=adata.obs.index,
        columns=genes
    )
    
    # Extract distance vector
    distance = adata.obs[distance_col].values
    
    results = []
    for gene in genes:
        y = expr[gene].values
        # Spearman correlation
        rho, rho_p = spearmanr(distance, y)
        # Pearson correlation
        r, r_p = pearsonr(distance, y)
        
        results.append({
            "Gene": gene,
            "Spearman_rho": rho,
            "Spearman_pval": rho_p,
            "Pearson_r": r,
            "Pearson_pval": r_p
        })
        
    return pd.DataFrame(results).sort_values(sortby, ascending=True)

def plot_signature_heatmap_matplotlib(scores_matrix, figsize=(12,8), cmap='RdBu_r', center=0):

    import matplotlib.pyplot as plt
    import numpy as np
    """
    Plot a heatmap of signature scores per cluster using matplotlib tiles.
    
    Parameters
    ----------
    scores_matrix : pd.DataFrame
        Rows = signatures, Columns = clusters, Values = scores
    figsize : tuple
        Width, height
    cmap : str
        Color map
    center : float
        Value to center the colormap
    """
    data = scores_matrix.values
    signatures = scores_matrix.index
    clusters = scores_matrix.columns

    # Create figure and axes
    plt.figure(figsize=figsize)
    im = plt.imshow(data, aspect='auto', cmap=cmap)

    # Center the colormap if needed
    if center is not None:
        from matplotlib import colors
        im.set_norm(colors.TwoSlopeNorm(vmin=data.min(), vcenter=center, vmax=data.max()))

    # Add colorbar
    cbar = plt.colorbar(im)
    cbar.set_label('Score', fontsize=12)

    # Set ticks
    plt.xticks(ticks=np.arange(len(clusters)), labels=clusters, rotation=45, fontsize=12)
    plt.yticks(ticks=np.arange(len(signatures)), labels=signatures, fontsize=12)

    # Add grid lines
    plt.grid(False)
    
    # Optional: annotate cells with numbers
    for i in range(len(signatures)):
        for j in range(len(clusters)):
            plt.text(j, i, f"{data[i, j]:.2f}", ha='center', va='center', fontsize=9, color='black')

    plt.xlabel('Cluster', fontsize=14)
    plt.ylabel('Signature', fontsize=14)
    plt.title('T cell Signature Scores per Cluster', fontsize=16, pad=20)
    plt.tight_layout()
    plt.show()

def plot_gene_distance_heatmap(
    adata,
    genes_ordered,
    gene_groups=None,      
    x_col = "spatial_bin",
    tit = "heatmap",
    ph=0,pw=0
):
    expr = pd.DataFrame(
        adata[:, genes_ordered].X.toarray(),
        index=adata.obs.index,
        columns=genes_ordered
    )
    expr[x_col] = adata.obs[x_col].values
    binned_expr = expr.groupby(x_col)[genes_ordered].mean()
    binned_expr_mm = (binned_expr - binned_expr.min()) / (binned_expr.max() - binned_expr.min())
    #norm = np.linalg.norm(binned_expr)
    #binned_expr_mm = binned_expr / norm if norm != 0 else binned_expr

    #binned_expr_mm = (binned_expr - binned_expr.mean()) / binned_expr.std()
    #binned_expr_mm = binned_expr
    if gene_groups is not None:
        group_palette = sns.color_palette("tab10", n_colors=len(gene_groups))
        group_color_map = dict(zip(gene_groups.keys(), group_palette))
        ytick_colors = []
        # For row group annotation
        row_group_positions = {}
        start = 0
        for group, g_genes in gene_groups.items():
            positions = [i for i, g in enumerate(binned_expr.columns) if g in g_genes]
            if positions:
                row_group_positions[group] = (min(positions), max(positions))
                ytick_colors.extend([group_color_map[group]]*len(positions))
        # If some genes not in groups
        if len(ytick_colors) < len(binned_expr.columns):
            ytick_colors.extend(["black"]*(len(binned_expr.columns)-len(ytick_colors)))
    else:
        ytick_colors = ["black"]*binned_expr.shape[1]

    # Column colors for response
    #col_meta = expr.groupby(x_col).agg(lambda x: x.mode()[0])
    col_meta = expr.groupby(x_col).agg(lambda x: x.mode().iloc[0] if not x.mode().empty else np.nan)
    #xtick_colors = col_meta.map(response_palette).values
    xtick_colors = ["black"]*binned_expr.shape[0]

    # -------------------
    # 3️⃣ Plot heatmap
    # -------------------
    if ph+pw == 0:
        pw = max(10, len(col_meta)*0.2)
        ph = len(genes_ordered)*0.3

    plt.figure(figsize=(pw, ph))
    ax = sns.heatmap(
        binned_expr_mm.T,
        cmap="viridis",
        cbar_kws={"label": "Mean Expression", "orientation": "horizontal"}
    )

    # Move colorbar to bottom
    cbar = ax.collections[0].colorbar
    cbar.ax.set_position([0.1, 0.5, 0.8, 0.03])
    #cbar.ax.set_position([0.1, 0.05, 0.8, 0.03])
    cbar.ax.xaxis.set_label_position('top')

    # Force all yticks (genes) to be shown
    ax.set_yticks(range(len(binned_expr.columns)))
    ax.set_yticklabels(binned_expr.columns, fontsize=8)

    # Force all xticks (samples/conditions) to be shown
    ax.set_xticks(range(len(binned_expr.index)))
    ax.set_xticklabels(binned_expr.index, fontsize=8, rotation=90)


    # Color ytick labels
    for tick_label, color in zip(ax.get_yticklabels(), ytick_colors):
        tick_label.set_color(color)
        tick_label.set_fontweight('bold')

    # Color xtick labels
    for tick_label, color in zip(ax.get_xticklabels(), xtick_colors):
        tick_label.set_color(color)
        tick_label.set_fontweight('bold')
        tick_label.set_rotation(90)

    ax.set_xlabel(f"Distance bins", fontsize=12, labelpad = 10)
    #ax.set_ylabel("Genes", fontsize=12)
    plt.title(tit, pad=20)
    plt.tight_layout()
    return plt, ax

def calc_tsub_proportion(adata, prop_obs, gropobs, meta_cols=None):
    grp_list = adata.obs[gropobs].unique().tolist()
    res = {}
    meta_data = {}
    
    for g in grp_list:
        df_subset = adata[adata.obs[gropobs] == g].copy()
        counts = df_subset.obs[prop_obs].value_counts(normalize=True)
        res[g] = counts
        if meta_cols:
            meta_data[g] = df_subset.obs[meta_cols].iloc[0].to_dict()
    res_df = pd.DataFrame(res).T.fillna(0)
    if meta_cols:
        meta_df = pd.DataFrame(meta_data).T
        res_df = pd.concat([meta_df, res_df], axis=1)
    return res_df

def calc_tsub_raw(adata, prop_obs, gropobs, meta_cols=None):
    grp_list = adata.obs[gropobs].unique().tolist()
    res = {}
    meta_data = {}
    
    for g in grp_list:
        df_subset = adata[adata.obs[gropobs] == g].copy()
        counts = df_subset.obs[prop_obs].value_counts()
        res[g] = counts
        if meta_cols:
            meta_data[g] = df_subset.obs[meta_cols].iloc[0].to_dict()
    res_df = pd.DataFrame(res).T.fillna(0)
    res_df["Total_Cells"] = res_df.sum(axis=1)
    if meta_cols:
        meta_df = pd.DataFrame(meta_data).T
        res_df = pd.concat([meta_df, res_df], axis=1)
    return res_df

def stratify_spatial_bins(adata, bins, labels, distance_col="Distance_to_Tumor", tn_col="TN_status:Tumor",
                       new_col="spatial_bin", target_ct = "Tumor", TN_flag = True):
    adata.obs[new_col] = pd.cut(
        adata.obs[distance_col],
        bins=bins,
        labels=labels,
        right=False
    ).astype(str)
    if (TN_flag == True):
        adata.obs.loc[adata.obs[tn_col] == "TN of " + target_ct, new_col] = "Q0: TN"
    return adata



from sklearn.decomposition import TruncatedSVD
from scipy.cluster.hierarchy import linkage, fcluster

def assign_spatial_niches(adata,pure_threshold=0.9,n_clusters=10,svd_components=5):
    neighbor_df = pd.DataFrame(
        adata.X.toarray() if hasattr(adata.X, "toarray") else adata.X,
        index=adata.obs_names,
        columns=adata.var_names
    )
    neighbor_df = neighbor_df.div(neighbor_df.sum(axis=1).replace(0,1), axis=0).fillna(0)
    
    # Layer 1: pure niches
    dominant_type = neighbor_df.idxmax(axis=1)
    dominant_fraction = neighbor_df.max(axis=1)
    
    spatial_labels = np.where(
        dominant_fraction >= pure_threshold,
        dominant_type + "_rich_niche",
        "mixed"
    )
    
    # Layer 2: hierarchical clustering on mixed cells
    mixed_mask = spatial_labels == "mixed"
    X_mixed = neighbor_df.loc[mixed_mask]
    
    if X_mixed.shape[0] > 0:
        n_comp = min(svd_components, X_mixed.shape[1]-1)
        if n_comp > 0:
            svd = TruncatedSVD(n_components=n_comp, random_state=0)
            X_reduced = svd.fit_transform(X_mixed)
        else:
            X_reduced = X_mixed.values.copy()
        
        Z = linkage(X_reduced, method="ward")
        cluster_labels = fcluster(Z, n_clusters, criterion="maxclust")
        
        spatial_labels[mixed_mask] = [f"Mixed_niche_{i}" for i in cluster_labels]
    
    adata.obs['spatial_niche'] = pd.Series(spatial_labels, index=adata.obs_names).astype("category")
    return adata

def relabel_mixed_niches(adata,dominant_thresh=0.6,max_top_neighbors=3):
    # Compute mean neighbor fractions per cluster
    adata.X = adata.X / adata.X.sum(axis=1, keepdims=True)
    cluster_means = adata.to_df().groupby(adata.obs['spatial_niche']).mean()
    # cluster_means: rows = cells, columns = neighbors → actually columns = features
    # For iteration, columns = clusters, rows = neighbors if you transpose:
    cluster_means = cluster_means.T

    import seaborn as sns
    import matplotlib.pyplot as plt
    plt.figure(figsize=(15,6))
    sns.heatmap(cluster_means, cmap="viridis", annot=True)
    plt.ylabel("Niche")
    plt.xlabel("Neighbor Cell Types")
    plt.show()
    new_labels = {}
    for cluster in cluster_means.columns:
        if cluster.lower().startswith("mixed_niche"):
            row = cluster_means[cluster].sort_values(ascending=False)
            top_neighbors = row.index[:max_top_neighbors]
            top_fracs = row.values[:max_top_neighbors]
            top1_frac = top_fracs[0]
            top2_frac = top_fracs[1] if len(top_fracs) > 1 else 0.0
            
            # Rule 1: dominant threshold
            if top1_frac >= dominant_thresh:
                new_label = f"{top_neighbors[0]}_rich_niche"
            # Rule 2: top gap threshold
            elif (top1_frac - top2_frac) >= 0.8 * top1_frac:
                new_label = f"{top_neighbors[0]}_rich_niche"
            else:
                # Rule 3: cumulative sum
                cum_sum = 0
                chosen_neighbors = []
                for n, f in row.items():
                    cum_sum += f
                    chosen_neighbors.append(n)
                    if cum_sum >= dominant_thresh or len(chosen_neighbors) == max_top_neighbors:
                        break
                new_label = "_".join(chosen_neighbors) + "_mixed"
            print(f"{cluster} → {new_label} (top fracs: {top_fracs})")
            new_labels[cluster] = new_label
        else:
            new_labels[cluster] = cluster
    # Apply new labels
    adata.obs['spatial_niche'] = adata.obs['spatial_niche'].replace(new_labels).astype('category')
    
    return adata

import pandas as pd
import numpy as np
from scipy.stats import chi2_contingency, fisher_exact
def calc_chi_pval(count1, total1, count2, total2):
    """
    Calculate Chi-squared or Fisher exact test p-value for 2x2 table.
    
    Parameters:
    count1, total1 : int : count and total for group 1
    count2, total2 : int : count and total for group 2
    
    Returns:
    dict : {'p_value': float, 'method': str, 'statistic': float or None}
    """
    # Build 2x2 table
    table = np.array([
        [count1, total1 - count1],
        [count2, total2 - count2]
    ])
    
    # Check expected counts
    chi2_test = chi2_contingency(table, correction=False)
    if (chi2_test[3] < 5).any():
        # Use Fisher if any expected count < 5
        stat, pval = fisher_exact(table)
        method = "Fisher"
        statistic = None
    else:
        # Use Chi-squared otherwise
        chi2_stat, pval, _, _ = chi2_contingency(table, correction=False)
        method = "Chi-squared"
        statistic = chi2_stat
        
    return {'p_value': pval, 'method': method, 'statistic': statistic}

from statsmodels.stats.proportion import proportions_ztest

def calc_prop_test(count1, total1, count2, total2):
    if any([x is None or np.isnan(x) for x in [count1, total1, count2, total2]]) or total1 == 0 or total2 == 0:
        return {'statistic': np.nan, 'p_value': np.nan, 'method': 'proportion_ztest'}
    
    counts = np.array([count1, count2])
    nobs   = np.array([total1, total2])

    z_stat, pval = proportions_ztest(counts, nobs)

    return {'p_value': pval, 'method': 'proportion_ztest', 'statistic': z_stat}


import pandas as pd
import numpy as np
from scipy.stats import chi2_contingency, fisher_exact

def calc_tcell_type_pvals(df_counts, id_vars=None, group_col='TumorType', group1='Pre', group2='Post'):
    df_long = df_counts.melt(
        id_vars=[col for col in df_counts.columns if not col.startswith("CD8_")],
        value_vars=[col for col in df_counts.columns if col.startswith("CD8_")],
        var_name="Subtype",
        value_name="Count"
    )
    results = []
    for keys, grp in df_long.groupby(id_vars):
        # Counts & totals
        count1  = grp.loc[grp[group_col]==group1, 'Count'].values[0] if group1 in grp[group_col].values else np.nan
        total1  = grp.loc[grp[group_col]==group1, 'Total_Cells'].values[0] if group1 in grp[group_col].values else np.nan
        prop1   = count1 / total1 if total1 > 0 else np.nan

        count2  = grp.loc[grp[group_col]==group2, 'Count'].values[0] if group2 in grp[group_col].values else np.nan
        total2  = grp.loc[grp[group_col]==group2, 'Total_Cells'].values[0] if group2 in grp[group_col].values else np.nan
        prop2   = count2 / total2 if total2 > 0 else np.nan
        # p-value
        if not np.isnan(count1) and not np.isnan(count2):
            pval_res = calc_chi_pval(count1, total1, count2, total2)
        else:
            pval_res = {'p_value': np.nan, 'method': None, 'statistic': None}

        result_dict = {col: val for col, val in zip(id_vars, keys)}  # retain all strata values
        result_dict.update({
            f'Count_{group1}': count1,
            f'Count_{group2}': count2,
            f'Prop_{group1}': prop1,
            f'Prop_{group2}': prop2,
            'p_value': pval_res['p_value'],
            'test_method': pval_res['method'],
            'statistic': pval_res['statistic']
        })
        results.append(result_dict)

    return pd.DataFrame(results)

def calc_tcell_composition_change_pval(df_counts, id_vars=None, group_col='TumorType', group1='Pre', group2='Post' ):
    df = df_counts.copy()
    subtype_cols = [col for col in df.columns if col.startswith("CD8_")]
    df[subtype_cols] = df[subtype_cols].fillna(0)
    results = []
    for keys, g in df.groupby(id_vars):
        g = g.sort_values(group_col)
        if g[group_col].nunique() < 2:
            continue
        contingency = g[subtype_cols].to_numpy()
        if contingency.shape[0] != 2:
            continue
        if contingency.shape == (2, 2) and (contingency < 5).any():
            _, p = fisher_exact(contingency)
            test = "fisher"
        else:
            chi2, p, dof, exp = chi2_contingency(contingency)
            test = "chi2"

        result_dict = {col: val for col, val in zip(id_vars, keys)}  # retain all strata values
        result_dict.update({
            f'Total_{group1}': contingency[0].sum(),
            f'Count_{group2}': contingency[1].sum(),
            'p_value': p,
            'test_method': test
        })
        results.append(result_dict)

    results_df = pd.DataFrame(results)
    return results_df
