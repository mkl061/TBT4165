import networkx as nx
import numpy as np
import copy

def s_core(G, s, genes, corrmat, pmat):
    """ Returns the s-core of network G.

    Parameters
    ----------
    G : NetworkX graph
        A graph object of the yeast PIN. Can either be a previously generated s-core or
        the full yeast PIN.
    s : float
        Strength cutoff for node deletion.
    genes: list
        List of genes (same order as rows and columns
        of corrmat and pmat).
    corrmat:
        Correlation matrix
    pmat:
        P-value matrix

    Returns
    -------
    s_core : NetworkX graph
        s-core of network G.
    """
    s_core = copy.deepcopy(G)
    n_genes = len(genes)
    p = 1

    # Add correlations as edge attributes if not already implemented
    if not nx.get_edge_attributes(s_core, "strength"):
        nx.set_edge_attributes(s_core, 0, "strength")    # default for all edges
        for i in range(n_genes):
            for j in range(i + 1, n_genes):
                gene_i = genes[i]
                gene_j = genes[j]

                # Check if edge exist, add correlation as edge strength if rho > 0 and p-val < 1
                if s_core.has_edge(gene_i, gene_j) and corrmat[i, j] > 0 and pmat[i, j] < p:
                    attr = {(gene_i, gene_j): {"strength": corrmat[i, j]}}
                    nx.set_edge_attributes(s_core, attr)

    n_old = s_core.number_of_nodes()    # current number of nodes
    n = 0   # number of nodes after node removal, set to 0 for initial iteration

    # Prune network
    while n_old != n:
        nodes = s_core.nodes()
        n_old = s_core.number_of_nodes()

        # Remove node if sum of strengths are smaller than cutoff s
        del_nodes = []
        for node in nodes:
            edges = list(s_core.edges(node))
            strengths = sum([s_core[edge[0]][edge[1]]["strength"] for edge in edges])
            if strengths < s:
                del_nodes.append(node)
        s_core.remove_nodes_from(del_nodes)
        n = s_core.number_of_nodes()     # number of nodes after node removal

    return s_core
