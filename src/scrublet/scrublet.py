from .helper_functions import *
from sklearn.decomposition import PCA, TruncatedSVD

def compute_doublet_scores(E, n_neighbors=50, sim_doublet_ratio=3, expected_doublet_rate = 0.1, use_approx_neighbors=True, total_counts_normalize=True, min_counts=3, min_cells=3, vscore_percentile=85, gene_filter=None, scaling_method = 'zscore', n_prin_comps=30, get_doublet_neighbor_parents = False):
    ''' Predict cell doublets

    Given a counts matrix `E`, calculates a doublet score between 0 and 1 by 
    simulating doublets, performing PCA, building a k-nearest-neighbor graph, 
    and finding the fraction of each observed transcriptome's neighbors that are
    simulated doublets. This 
    
    Required inputs:
    - E: 2-D matrix with shape (n_cells, n_genes)
        scipy.sparse matrix or numpy array containing raw (unnormalized) 
        transcript counts. E[i,j] is the number of copies of transcript j 
        detected in cell i.

    Optional parameters:
    - n_neighbors: int (default = 50)
        Number of neighbors used to construct the kNN graph of observed
        transcriptomes and simulated doublets.
    - sim_doublet_ratio: float(default = 3)
        Number of doublets to simulate relative to the number of observed 
        transcriptomes.
    - expected_doublet_rate: float (default = 0.1)
        The estimated doublet rate for the experiment.
    - use_approx_neighbors: boolean (default = True)
        If true, use approximate nearest neighbors algorithm to contstruct the 
        kNN graph.
    - total_counts_normalize: boolean (default = True)
        If true, apply total counts normalization prior to PCA
    - gene_filter: list (default = None)
        For gene filtering prior to PCA. List of gene indices (columns of `E`) 
        to use in PCA.
    - min_counts and min_cells: float (default = 3 for both)
        For gene filtering prior to PCA. Keep genes with at least `min_counts` 
        counts in at least `min_cells` cells. Ignored if `gene_filter` is not
        `None`.
    - vscore_percentile: float (default = 85)
        For gene filtering prior to  PCA. Keep only highly variable genes, i.e.,
        those in the `vscore_percentile`-th percentile or above. Ignored if 
        `gene_filter` is not `None`.
    - scaling_method: str (default = "zscore")
        Method for gene-level scaling of transcript counts prior to PCA. Options
        are "zscore", "log", and "none".
    - n_prin_comps: int (default = 30)
        Number of principal components to use for embedding cells prior to
        building kNN graph.
    - get_doublet_neighbor_parents: boolean (default = False)
        If true, returns the list of parent cells used to generate each observed
        transcriptome's doublet neighbors.

    Returns: a dictionary with following items:
    - "doublet_scores_observed_cells": 1-D `array`
        List of doublet scores for each observed transcriptome
    - "doublet_scores_simulateed_doublets": 1-D `array`
        List of doublet scores for each synthetic doublet
    - "pca_observed_cells": 2-D `array`
        Principal component embedding of the observed transcriptomes
    - "pca_simulated_doublets":  
        Principal component embedding of the simulated doublets
    - "gene_filter":
        List of gene indices used for PCA.
    - "doublet_neighbor_parents": list of numpy arrays
        `None` if `get_doublet_neighbor_parents` is `False`. Otherwise, entry i
        is the list (1-D numpy array) of parent cells that generated the doublet 
        neighbors of cell i. Cells with no doublet neighbors have an empty list.
    '''

    # Initialize output: dictionary to store results 
    # and useful intermediate variables
    output = {}

    # Check that input is valid
    valid_scaling_methods = ['zscore', 'log', 'none']
    if not scaling_method in valid_scaling_methods:
        print('Select one of the following scaling methods:', valid_scaling_methods)
        return

    # Convert counts matrix to sparse format if necessary
    if not scipy.sparse.issparse(E):
        print('Converting to scipy.sparse.csc_matrix')
        E = scipy.sparse.csc_matrix(E)
    elif not scipy.sparse.isspmatrix_csc(E):
        print('Converting to scipy.sparse.csc_matrix')
        E = E.tocsc()

    # Simulate doublets 
    print('Simulating doublets')
    E_doub, parent_ix = simulate_doublets_from_counts(E, sim_doublet_ratio)

    # Total counts normalize observed cells and simulated doublets
    if total_counts_normalize:
        print('Total counts normalizing')
        E = tot_counts_norm(E)[0] 
        E_doub = tot_counts_norm(E_doub, target_mean=1e5)[0]

    # Filter genes (highly variable, expressed above minimum level)
    if gene_filter is None:
        print('Finding highly variable genes')
        gene_filter = filter_genes(E, min_vscore_pctl=vscore_percentile, min_counts=min_counts, min_cells=min_cells)
    else:
        gene_filter = np.array(gene_filter)

    # Total counts normalize observed cells to the same total as doublets
    if total_counts_normalize:
        E = tot_counts_norm(E, target_mean=1e5)[0]

    # Apply gene filter
    print('Filtering genes from {} to {}'.format(E.shape[1], len(gene_filter)))
    E = E[:, gene_filter]
    E_doub = E_doub[:, gene_filter]
    output['gene_filter'] = gene_filter

    # Rescale counts
    if scaling_method == 'log':
        print('Applying log normalization')
        E.data = np.log10(1 + E.data)
        E_doub.data = np.log10(1 + E_doub.data)
        # to do: option of TruncatedSVD to preserve sparsity
        pca = PCA(n_components = n_prin_comps)
        E_doub = E_doub.toarray()
        E = E.toarray()
    elif scaling_method == 'zscore':
        print('Applying z-score normalization')
        gene_means = E.mean(0)
        gene_stdevs = np.sqrt(sparse_var(E))
        E = sparse_zscore(E, gene_means, gene_stdevs)
        E_doub = sparse_zscore(E_doub, gene_means, gene_stdevs)
        pca = PCA(n_components = n_prin_comps)
    else:
        # to do: option of TruncatedSVD to preserve sparsity
        pca = PCA(n_components = n_prin_comps)
        E_doub = E_doub.toarray()
        E = E.toarray()

    # Fit PCA to observed cells, then apply same transformation to
    # simulated doublets.
    print('Running PCA')
    pca.fit(E)
    E_pca = pca.transform(E)
    E_doub_pca = pca.transform(E_doub)
    doub_labels = np.concatenate((np.zeros(E_pca.shape[0], dtype=int), 
                                  np.ones(E_doub_pca.shape[0], dtype=int)))

    output['pca_observed_cells'] = E_pca
    output['pca_simulated_cells'] = E_doub_pca

    # Calculate doublet scores using k-nearest-neighbor classifier
    print('Building kNN graph and calculating doublet scores')
    nn_outs = nearest_neighbor_classifier(np.vstack((E_pca, E_doub_pca)), 
                                          doub_labels, 
                                          k=n_neighbors, 
                                          use_approx_nn=use_approx_neighbors, 
                                          exp_doub_rate = expected_doublet_rate, 
                                          get_neighbor_parents = get_doublet_neighbor_parents,
                                          parent_cells = parent_ix
                                         ) 
    output['doublet_scores_observed_cells'] = nn_outs[0]
    output['doublet_scores_simulated_doublets'] = nn_outs[1]
    output['doublet_neighbor_parents'] = nn_outs[2]
    return output

#========================================================================================#

def plot_scrublet_results(coords, doublet_scores_obs, doublet_scores_sim, score_threshold, marker_size=5, order_points=False, scale_hist_obs='log', scale_hist_sim='linear', fig_size = (8,6)):

    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    called_doubs = doublet_scores_obs > score_threshold
    called_doubs_sim = doublet_scores_sim > score_threshold
    predictable_doub_frac = sum(called_doubs_sim) / float(len(called_doubs_sim))
    called_frac = sum(called_doubs) / float(len(called_doubs))

    print('{}/{} = {:.1f}% of cells are predicted doublets.'.format(sum(called_doubs), len(called_doubs), 
                                                              100 * called_frac))
    print('{:.1f}% of doublets are predicted to be detectable.'.format(100 * predictable_doub_frac))
    print('Predicted overall doublet rate = {:.1f}%'.format(100 * called_frac / predictable_doub_frac))    

    fig, axs = plt.subplots(2, 2, figsize = fig_size)

    ax = axs[0,0]
    ax.hist(doublet_scores_obs, np.linspace(0, 1, 50), color = 'gray', linewidth = 0, density=True)
    ax.set_yscale(scale_hist_obs)
    yl = ax.get_ylim()
    ax.set_ylim(yl)
    ax.plot([score_threshold, score_threshold], yl, c = 'black', linewidth = 1)
    ax.set_title('Observed cells')
    ax.set_xlabel('Doublet score')
    ax.set_ylabel('Prob. density')

    ax = axs[0,1]
    ax.hist(doublet_scores_sim, np.linspace(0, 1, 50), color = 'gray', linewidth = 0, density=True)
    ax.set_yscale(scale_hist_sim)
    yl = ax.get_ylim()
    ax.set_ylim(yl)
    ax.plot([score_threshold, score_threshold], yl, c = 'black', linewidth = 1)
    ax.set_title('Simulated doublets')
    ax.set_xlabel('Doublet score')
    ax.set_ylabel('Prob. density')

    x = coords[:,0]
    y = coords[:,1]
    xl = (x.min() - x.ptp() * .05, x.max() + x.ptp() * 0.05)
    yl = (y.min() - y.ptp() * .05, y.max() + y.ptp() * 0.05)
    
    if order_points:
        o = np.argsort(doublet_scores_obs)
    else:
        o = np.arange(len(doublet_scores_obs)) 
    
    ax = axs[1,0]
    pp = ax.scatter(x[o], y[o], s=marker_size, edgecolors='', c = doublet_scores_obs[o], cmap=darken_cmap(plt.cm.Reds, 0.9))
    ax.set_xlim(xl)
    ax.set_ylim(yl)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_title('Doublet score')

    ax = axs[1,1]
    ax.scatter(x[o], y[o], s=marker_size, edgecolors='', c=doublet_scores_obs[o] > score_threshold, cmap = custom_cmap([[.7,.7,.7], [0,0,0]]))
    ax.set_xlim(xl)
    ax.set_ylim(yl)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_title('Predicted doublets')
    singlet_marker = Line2D([], [], color=[.7,.7,.7], marker='o', markersize=5, label='Singlet', linewidth=0)
    doublet_marker = Line2D([], [], color=[.0,.0,.0], marker='o', markersize=5, label='Doublet', linewidth=0)
    ax.legend(handles = [singlet_marker, doublet_marker])

    fig.tight_layout()

    return fig, axs

#========================================================================================#

def simulate_doublets_from_counts(E, sim_doublet_ratio=1):
    '''
    Simulate doublets by summing the counts of random cell pairs.

    Inputs:
    E (numpy or scipy matrix of size (num_cells, num_genes)): counts matrix, ideally without total-counts normalization.
    sim_doublet_ratio (float): number of doublets to simulate, as a fraction of the number of cells in E.
                          A total of num_sim_doubs = int(sim_doublet_ratio * E[0]) doublets will be simulated.

    Returns:
    - Edoub (scipy sparse CSC matrix of size (num_cells+num_sim_doubs, num_genes)): counts matrix with the simulated doublet data appended to the original data matrix E.
    - doub_labels (array of size (num_cells+num_sim_doubs)): 0 if observed cell, 1 if simulated doublet
    - pair_ix (matrix of size(num_sim_doubs, 2)): each row gives the indices of the parent cells from E used to generate the simulated doublet
    '''

    if not scipy.sparse.issparse(E):
        E = scipy.sparse.csc_matrix(E)
    elif not scipy.sparse.isspmatrix_csc(E):
        E = E.tocsc()

    n_obs = E.shape[0]
    n_doub = int(n_obs * sim_doublet_ratio)
    pair_ix = np.random.randint(0, n_obs, size=(n_doub, 2))
    Edoub = E[pair_ix[:, 0],:] + E[pair_ix[:, 1],:]

    return Edoub, pair_ix

#========================================================================================#

def simulate_doublets_from_pca(PCdat, total_counts=None, sim_doublet_ratio=1):
    '''
    Simulate doublets by averaging PCA coordinates of random cell pairs.
    Average is weighted by total counts of each parent cell, if provided.

    Returns:
    PCdoub (matrix of size (num_cells+num_sim_doubs, num_pcs)): PCA matrix with the simulated doublet PCA coordinates appended to the original data matrix PCdat.
    doub_labels (array of size (num_cells+num_sim_doubs)): 0 if observed cell, 1 if simulated doublet
    pair_ix (matrix of size(num_sim_doubs, 2)): each row gives the indices of the parent cells used to generate the simulated doublet
    '''

    n_obs = PCdat.shape[0]
    n_doub = int(n_obs * sim_doublet_ratio)

    if total_counts is None:
        total_counts = np.ones(n_obs)

    pair_ix = np.random.randint(0, n_obs, size=(n_doub, 2))

    pair_tots = np.hstack((total_counts[pair_ix[:, 0]][:,None], total_counts[pair_ix[:, 1]][:,None]))
    pair_tots = np.array(pair_tots, dtype=float)
    pair_fracs = pair_tots / np.sum(pair_tots, axis=1)[:,None]

    PCdoub = PCdat[pair_ix[:, 0],:] * pair_fracs[:, 0][:,None] + PCdat[pair_ix[:, 1],:] * pair_fracs[:, 1][:,None]

    PCdoub = np.vstack((PCdat, PCdoub))
    doub_labels = np.concatenate((np.zeros(n_obs, dtype=int), np.ones(n_doub, dtype=int)))

    return PCdoub, doub_labels, pair_ix

#========================================================================================#

def nearest_neighbor_classifier(embedding, doub_labels, k=50, use_approx_nn=True, exp_doub_rate = 1.0, get_neighbor_parents = False, parent_cells = None):
    n_obs = sum(doub_labels == 0)
    n_sim = sum(doub_labels == 1)

    # Adjust k (number of nearest neighbors) based on the ratio of simulated to observed cells
    k_adj = int(round(k * (1+n_sim/float(n_obs))))

    # Find k_adj nearest neighbors
    neighbors = get_knn_graph(embedding, k=k_adj, dist_metric='euclidean', approx=use_approx_nn, return_edges = False)
    
    # Calculate doublet score based on ratio of simulated cell neighbors vs. observed cell neighbors
    doub_neigh_mask = doub_labels[neighbors] == 1
    n_sim_neigh = doub_neigh_mask.sum(1)
    n_obs_neigh = doub_neigh_mask.shape[1] - n_sim_neigh
    doub_score = n_sim_neigh / (n_sim_neigh + n_obs_neigh * n_sim / float(n_obs) / exp_doub_rate)

    # get parents of doublet neighbors, if requested
    neighbor_parents = None
    if get_neighbor_parents and parent_cells is not None:
        neighbors = neighbors - n_obs
        neighbor_parents = []
        for iCell in range(n_obs):
            this_doub_neigh = neighbors[iCell,:][neighbors[iCell,:] > -1]
            if len(this_doub_neigh) > 0:
                this_doub_neigh_parents = np.unique(parent_cells[this_doub_neigh,:].flatten())
                neighbor_parents.append(this_doub_neigh_parents)
            else:
                neighbor_parents.append([])
        neighbor_parents = np.array(neighbor_parents)

    
    return doub_score[doub_labels == 0], doub_score[doub_labels == 1], neighbor_parents

