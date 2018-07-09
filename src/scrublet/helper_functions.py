import numpy as np
import scipy
import scipy.stats
import scipy.sparse 
from sklearn.decomposition import PCA,TruncatedSVD
from sklearn.neighbors import NearestNeighbors
import time

########## LOADING DATA
def load_genes(filename, delimiter='\t', column=0, skip_rows=0):
    gene_list = []
    gene_dict = {}

    with open(filename) as f:
        for iL in range(skip_rows):
            f.readline()
        for l in f:
            gene = l.strip('\n').split(delimiter)[column]
            if gene in gene_dict:
                gene_dict[gene] += 1
                gene_list.append(gene + '__' + str(gene_dict[gene]))
                if gene_dict[gene] == 2:
                    i = gene_list.index(gene)
                    gene_list[i] = gene + '__1'
            else: 
                gene_dict[gene] = 1
                gene_list.append(gene)
    return gene_list


def make_genes_unique(orig_gene_list):
    gene_list = []
    gene_dict = {}

    for gene in orig_gene_list:
        if gene in gene_dict:
            gene_dict[gene] += 1
            gene_list.append(gene + '__' + str(gene_dict[gene]))
            if gene_dict[gene] == 2:
                i = gene_list.index(gene)
                gene_list[i] = gene + '__1'
        else:
           gene_dict[gene] = 1
           gene_list.append(gene)
    return gene_list

########## USEFUL SPARSE FUNCTIONS

def sparse_var(E, axis=0):
    ''' variance across the specified axis '''

    mean_gene = E.mean(axis=axis).A.squeeze()
    tmp = E.copy()
    tmp.data **= 2
    return tmp.mean(axis=axis).A.squeeze() - mean_gene ** 2

def sparse_multiply(E, a):
    ''' multiply each row of E by a scalar '''

    nrow = E.shape[0]
    w = scipy.sparse.lil_matrix((nrow, nrow))
    w.setdiag(a)
    return w * E

def sparse_zscore(E):
    ''' z-score normalize each column of E '''

    mean_gene = E.mean(0)
    stdev_gene = np.sqrt(sparse_var(E))
    return sparse_multiply((E - mean_gene).T, 1/stdev_gene).T

########## GENE FILTERING

def runningquantile(x, y, p, nBins):

    ind = np.argsort(x)
    x = x[ind]
    y = y[ind]

    dx = (x[-1] - x[0]) / nBins
    xOut = np.linspace(x[0]+dx/2, x[-1]-dx/2, nBins)

    yOut = np.zeros(xOut.shape)

    for i in range(len(xOut)):
        ind = np.nonzero((x >= xOut[i]-dx/2) & (x < xOut[i]+dx/2))[0]
        if len(ind) > 0:
            yOut[i] = np.percentile(y[ind], p)
        else:
            if i > 0:
                yOut[i] = yOut[i-1]
            else:
                yOut[i] = np.nan

    return xOut, yOut


def get_vscores(E, min_mean=0, nBins=50, fit_percentile=0.1, error_wt=1):
    '''
    Calculate v-score (above-Poisson noise statistic) for genes in the input counts matrix
    Return v-scores and other stats
    '''

    ncell = E.shape[0]

    mu_gene = E.mean(axis=0).A.squeeze()
    gene_ix = np.nonzero(mu_gene > min_mean)[0]
    mu_gene = mu_gene[gene_ix]

    tmp = E[:,gene_ix]
    tmp.data **= 2
    var_gene = tmp.mean(axis=0).A.squeeze() - mu_gene ** 2
    del tmp
    FF_gene = var_gene / mu_gene

    data_x = np.log(mu_gene)
    data_y = np.log(FF_gene / mu_gene)

    x, y = runningquantile(data_x, data_y, fit_percentile, nBins)
    x = x[~np.isnan(y)]
    y = y[~np.isnan(y)]

    gLog = lambda input: np.log(input[1] * np.exp(-input[0]) + input[2])
    h,b = np.histogram(np.log(FF_gene[mu_gene>0]), bins=200)
    b = b[:-1] + np.diff(b)/2
    max_ix = np.argmax(h)
    c = np.max((np.exp(b[max_ix]), 1))
    errFun = lambda b2: np.sum(abs(gLog([x,c,b2])-y) ** error_wt)
    b0 = 0.1
    b = scipy.optimize.fmin(func = errFun, x0=[b0], disp=False)
    a = c / (1+b) - 1


    v_scores = FF_gene / ((1+a)*(1+b) + b * mu_gene);
    CV_eff = np.sqrt((1+a)*(1+b) - 1);
    CV_input = np.sqrt(b);

    return v_scores, CV_eff, CV_input, gene_ix, mu_gene, FF_gene, a, b

def filter_genes(E, base_ix = [], min_vscore_pctl = 85, min_counts = 3, min_cells = 3, show_vscore_plot = False, sample_name = ''):
    ''' 
    Filter genes by expression level and variability
    Return list of filtered gene indices
    '''

    if len(base_ix) == 0:
        base_ix = np.arange(E.shape[0])

    Vscores, CV_eff, CV_input, gene_ix, mu_gene, FF_gene, a, b = get_vscores(E[base_ix, :])
    ix2 = Vscores>0
    Vscores = Vscores[ix2]
    gene_ix = gene_ix[ix2]
    mu_gene = mu_gene[ix2]
    FF_gene = FF_gene[ix2]
    min_vscore = np.percentile(Vscores, min_vscore_pctl)
    ix = (((E[:,gene_ix] >= min_counts).sum(0).A.squeeze() >= min_cells) & (Vscores >= min_vscore))
    
    if show_vscore_plot:
        import matplotlib.pyplot as plt
        x_min = 0.5*np.min(mu_gene)
        x_max = 2*np.max(mu_gene)
        xTh = x_min * np.exp(np.log(x_max/x_min)*np.linspace(0,1,100))
        yTh = (1 + a)*(1+b) + b * xTh
        plt.figure(figsize=(8, 6));
        plt.scatter(np.log10(mu_gene), np.log10(FF_gene), c = [.8,.8,.8], alpha = 0.3, edgecolors='');
        plt.scatter(np.log10(mu_gene)[ix], np.log10(FF_gene)[ix], c = [0,0,0], alpha = 0.3, edgecolors='');
        plt.plot(np.log10(xTh),np.log10(yTh));
        plt.title(sample_name)
        plt.xlabel('log10(mean)');
        plt.ylabel('log10(Fano factor)');
        plt.show()

    return gene_ix[ix]

########## CELL NORMALIZATION

def tot_counts_norm(E, exclude_dominant_frac = 1, included = [], target_mean = 0):
    ''' 
    Cell-level total counts normalization of input counts matrix, excluding overly abundant genes if desired.
    Return normalized counts, average total counts, and (if exclude_dominant_frac < 1) list of genes used to calculate total counts 
    '''

    E = E.tocsc()
    ncell = E.shape[0]
    if len(included) == 0:
        if exclude_dominant_frac == 1:
            tots_use = E.sum(axis=1)
        else:
            tots = E.sum(axis=1)
            wtmp = scipy.sparse.lil_matrix((ncell, ncell))
            wtmp.setdiag(1. / tots)
            included = np.asarray(~(((wtmp * E) > exclude_dominant_frac).sum(axis=0) > 0))[0,:]
            tots_use = E[:,included].sum(axis = 1)
            print('Excluded %i genes from normalization' %(np.sum(~included)))
    else:
        tots_use = E[:,included].sum(axis = 1)

    if target_mean == 0:
        target_mean = np.mean(tots_use)

    w = scipy.sparse.lil_matrix((ncell, ncell))
    w.setdiag(float(target_mean) / tots_use)
    Enorm = w * E

    return Enorm.tocsc(), target_mean, included

########## DIMENSIONALITY REDUCTION

def get_pca(E, base_ix=[], numpc=50, keep_sparse=False, normalize=True):
    '''
    Run PCA on the counts matrix E, gene-level normalizing if desired
    Return PCA coordinates
    '''
    # If keep_sparse is True, gene-level normalization maintains sparsity
    #     (no centering) and TruncatedSVD is used instead of normal PCA.

    if len(base_ix) == 0:
        base_ix = np.arange(E.shape[0])

    if keep_sparse:
        if normalize:
            zstd = np.sqrt(sparse_var(E[base_ix,:]))
            Z = sparse_multiply(E.T, 1 / zstd).T
        else:
            Z = E
        pca = TruncatedSVD(n_components=numpc)

    else:
        if normalize:
            zmean = E[base_ix,:].mean(0)
            zstd = np.sqrt(sparse_var(E[base_ix,:]))
            Z = sparse_multiply((E - zmean).T, 1/zstd).T
        else:
            Z = E
        pca = PCA(n_components=numpc)

    pca.fit(Z[base_ix,:])
    return pca.transform(Z)


def preprocess_and_pca(E, total_counts_normalize=True, norm_exclude_abundant_gene_frac=1, min_counts=3, min_cells=5, min_vscore_pctl=85, gene_filter=None, num_pc=50, sparse_pca=False, show_vscore_plot=False):
    '''
    Total counts normalize, filter genes, run PCA
    Return PCA coordinates and filtered gene indices
    '''

    if total_counts_normalize:
        print('Total count normalizing')
        E = tot_counts_norm(E, exclude_dominant_frac = norm_exclude_abundant_gene_frac)[0]

    if gene_filter is None:
        print('Finding highly variable genes')
        gene_filter = filter_genes(E, min_vscore_pctl=min_vscore_pctl, min_counts=min_counts, min_cells=min_cells, show_vscore_plot=show_vscore_plot)

    print('Using %i genes for PCA' %len(gene_filter))
    PCdat = get_pca(E[:,gene_filter], numpc=num_pc, keep_sparse=sparse_pca)

    return PCdat, gene_filter

########## GRAPH CONSTRUCTION

def get_knn_graph(X, k=5, dist_metric='euclidean', approx=False, return_edges=True):
    '''
    Build k-nearest-neighbor graph
    Return edge list and nearest neighbor matrix
    '''

    t0 = time.time()
    if approx:
        try:
            from annoy import AnnoyIndex
        except:
            approx = False
            print('Could not find library "annoy" for approx. nearest neighbor search')
    if approx:
        print('Using approximate nearest neighbor search')

        if dist_metric == 'cosine':
            dist_metric = 'angular'
        npc = X.shape[1]
        ncell = X.shape[0]
        annoy_index = AnnoyIndex(npc, metric=dist_metric)

        for i in range(ncell):
            annoy_index.add_item(i, list(X[i,:]))
        annoy_index.build(10) # 10 trees

        knn = []
        for iCell in range(ncell):
            knn.append(annoy_index.get_nns_by_item(iCell, k + 1)[1:])
        knn = np.array(knn, dtype=int)

    else:
        print('Using sklearn NearestNeighbors')

        if dist_metric == 'cosine':
            nbrs = NearestNeighbors(n_neighbors=k, metric=dist_metric, algorithm='brute').fit(X)
        else:
            nbrs = NearestNeighbors(n_neighbors=k, metric=dist_metric).fit(X)
        knn = nbrs.kneighbors(return_distance=False)

    if return_edges:
        links = set([])
        for i in range(knn.shape[0]):
            for j in knn[i,:]:
                links.add(tuple(sorted((i,j))))

        t_elapse = time.time() - t0
        print('kNN graph built in %.3f sec' %(t_elapse))

        return links, knn
    return knn

def build_adj_mat(edges, n_nodes):
    A = scipy.sparse.lil_matrix((n_nodes, n_nodes))
    for e in edges:
        i, j = e
        A[i,j] = 1
        A[j,i] = 1
    return A.tocsc()

########## FORCE LAYOUT

def get_force_layout(links, n_cells, n_iter=100, edgeWeightInfluence=1, barnesHutTheta=2, scalingRatio=1, gravity=0.05, jitterTolerance=1, verbose=False):
    from fa2 import ForceAtlas2
    import networkx as nx

    G = nx.Graph()
    G.add_nodes_from(range(n_cells))
    G.add_edges_from(list(links))

    forceatlas2 = ForceAtlas2(
                  # Behavior alternatives
                  outboundAttractionDistribution=False,  # Dissuade hubs
                  linLogMode=False,  # NOT IMPLEMENTED
                  adjustSizes=False,  # Prevent overlap (NOT IMPLEMENTED)
                  edgeWeightInfluence=edgeWeightInfluence,

                  # Performance
                  jitterTolerance=jitterTolerance,  # Tolerance
                  barnesHutOptimize=True,
                  barnesHutTheta=barnesHutTheta,
                  multiThreaded=False,  # NOT IMPLEMENTED

                  # Tuning
                  scalingRatio=scalingRatio,
                  strongGravityMode=False,
                  gravity=gravity,
                  # Log
                  verbose=verbose)

    positions = forceatlas2.forceatlas2_networkx_layout(G, pos=None, iterations=n_iter)
    positions = np.array([positions[i] for i in sorted(positions.keys())])
    return positions

########## CLUSTERING

def get_spectral_clusters(A, k):
    from sklearn.cluster import SpectralClustering
    spec = SpectralClustering(n_clusters=k, random_state = 0, affinity = 'precomputed', assign_labels = 'discretize')
    return spec.fit_predict(A)


def get_louvain_clusters(nodes, edges):
    import networkx as nx
    import community
    
    G = nx.Graph()
    G.add_nodes_from(nodes)
    G.add_edges_from(edges)
    
    return np.array(community.best_partition(G).values())

########## GENE ENRICHMENT

def rank_enriched_genes(E, gene_list, cell_mask, min_counts=3, min_cells=3):
    gix = (E[cell_mask,:]>=min_counts).sum(0).A.squeeze() >= min_cells
    print('%i cells in group' %(sum(cell_mask)))
    print('Considering %i genes' %(sum(gix)))
    
    gene_list = gene_list[gix]
    
    z = sparse_zscore(E[:,gix])
    scores = z[cell_mask,:].mean(0).A.squeeze()
    o = np.argsort(-scores)
    
    return gene_list[o], scores[o]
    

########## PLOTTING STUFF

def darken_cmap(cmap, scale_factor):
    cdat = np.zeros((cmap.N, 4))
    for ii in range(cdat.shape[0]):
        curcol = cmap(ii)
        cdat[ii,0] = curcol[0] * scale_factor
        cdat[ii,1] = curcol[1] * scale_factor
        cdat[ii,2] = curcol[2] * scale_factor
        cdat[ii,3] = 1
    cmap = cmap.from_list(cmap.N, cdat)
    return cmap

def custom_cmap(rgb_list):
    import matplotlib.pyplot as plt
    rgb_list = np.array(rgb_list)
    cmap = plt.cm.Reds
    cmap = cmap.from_list(rgb_list.shape[0],rgb_list)
    return cmap

def plot_groups(x, y, groups, lim_buffer = 50, saving = False, fig_dir = './', fig_name = 'fig', res = 300, close_after = False, title_size = 12, point_size = 3, ncol = 5):
    import matplotlib.pyplot as plt

    n_col = int(ncol)
    ngroup = len(np.unique(groups))
    nrow = int(np.ceil(ngroup / float(ncol)))
    fig = plt.figure(figsize = (14, 3 * nrow))
    for ii, c in enumerate(np.unique(groups)):
        ax = plt.subplot(nrow, ncol, ii+1)
        ix = groups == c

        ax.scatter(x[~ix], y[~ix], s = point_size, c = [.8,.8,.8], edgecolors = '')
        ax.scatter(x[ix], y[ix], s = point_size, c = [0,0,0], edgecolors = '')
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_xlim([min(x) - lim_buffer, max(x) + lim_buffer])
        ax.set_ylim([min(y) - lim_buffer, max(y) + lim_buffer])

        ax.set_title(str(c), fontsize = title_size)

    fig.tight_layout()

    if saving:
        if not os.path.exists(fig_dir):
            os.makedirs(fig_dir)
        plt.savefig(fig_dir + '/' + fig_name + '.png', dpi=res)

    if close_after:
        plt.close()
