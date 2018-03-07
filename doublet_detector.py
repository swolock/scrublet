from helper_functions import *

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


    Edoub = scipy.sparse.vstack((E, Edoub))
    doub_labels = np.concatenate((np.zeros(n_obs), np.ones(n_doub)))

    return Edoub, doub_labels, pair_ix

#========================================================================================#

def simulate_doublets_from_pca(PCdat, total_counts=[], sim_doublet_ratio=1):
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

    if len(total_counts) == 0:
        total_counts = np.ones(n_obs)

    pair_ix = np.random.randint(0, n_obs, size=(n_doub, 2))

    pair_tots = np.hstack((total_counts[pair_ix[:, 0]][:,None], total_counts[pair_ix[:, 1]][:,None]))
    pair_tots = np.array(pair_tots, dtype=float)
    pair_fracs = pair_tots / np.sum(pair_tots, axis=1)[:,None]

    PCdoub = PCdat[pair_ix[:, 0],:] * pair_fracs[:, 0][:,None] + PCdat[pair_ix[:, 1],:] * pair_fracs[:, 1][:,None]

    PCdoub = np.vstack((PCdat, PCdoub))
    doub_labels = np.concatenate((np.zeros(n_obs), np.ones(n_doub)))

    return PCdoub, doub_labels, pair_ix

#========================================================================================#

def calculate_doublet_scores(embedding, doub_labels, k=50, use_approx_nn=True, exp_doub_rate = 1.0):
    t00 = time.time()
    n_obs = sum(doub_labels == 0)
    n_sim = sum(doub_labels == 1)

    # Adjust k (number of nearest neighbors) based on the ratio of simulated to observed cells
    k_adj = int(round(k * (1+n_sim/float(n_obs))))

    # Find k_adj nearest neighbors
    neighbors = get_knn_graph(embedding, k=k_adj, dist_metric='euclidean', approx=use_approx_nn)[1]
    
    # Calculate doublet score based on ratio of simulated cell neighbors vs. observed cell neighbors
    n_sim_neigh = np.sum(doub_labels[neighbors] == 1, axis = 1)
    n_obs_neigh = np.sum(doub_labels[neighbors] == 0, axis = 1)
    
    doub_score = n_sim_neigh / (n_sim_neigh + n_obs_neigh * n_sim / float(n_obs) / exp_doub_rate)
    doub_score_obs = doub_score[doub_labels == 0]

    # return doublet scores for observed cells and simulated cells
    return doub_score[doub_labels == 0], doub_score[doub_labels == 1]


#========================================================================================#

def woublet(E=None, exp_doub_rate = 0.1, sim_doublet_ratio=3, k=50, use_approx_nn=False, precomputed_pca=None, total_counts=None, total_counts_normalize=True, norm_exclude_abundant_gene_frac=1, min_counts=3, min_cells=5, vscore_percentile=85, gene_filter=None, num_pc=50, sparse_pca=False):
    
    # Check that input is valid
    if E is None and precomputed_pca is None:
        print 'Please supply a counts matrix (E) or PCA coordinates (precomputed_pca)'
        return
    
    # Convert counts matrix to sparse format if necessary
    if not E is None:
        if not scipy.sparse.issparse(E):
            E = scipy.sparse.csc_matrix(E)
        elif not scipy.sparse.isspmatrix_csc(E):
            E = E.tocsc()

    # Get total counts per cell
    if total_counts is None:
        if E is None:
            total_counts = np.ones(precomputed_pca.shape[0])
        else:
            total_counts = E.sum(1).A.squeeze()

    # Run PCA if not provided
    if precomputed_pca is None:
        PCdat = preprocess_and_pca(E, total_counts_normalize=total_counts_normalize, norm_exclude_abundant_gene_frac=norm_exclude_abundant_gene_frac, min_counts=min_counts, min_cells=min_cells, vscore_percentile=vscore_percentile, gene_filter=gene_filter, num_pc=num_pc, sparse_pca=sparse_pca)
    else:
        PCdat = precomputed_pca

    # Simulate doublets
    print 'Simulating doublets'
    PCdat, doub_labels, parent_ix = simulate_doublets_from_pca(PCdat, total_counts=total_counts, sim_doublet_ratio=sim_doublet_ratio)

    # Calculate doublet scores using k-nearest-neighbor classifier
    print 'Running KNN classifier'
    doub_score_obs, doub_score_sim = calculate_doublet_scores(PCdat, doub_labels, k=k, use_approx_nn=use_approx_nn, exp_doub_rate = exp_doub_rate)
    return doub_score_obs, doub_score_sim

