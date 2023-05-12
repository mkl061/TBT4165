import networkx as nx
import pandas as pd
import numpy as np
from scipy import special
import cobra
from builtins import object
import matplotlib.pyplot as plt
from palettable.colorbrewer import get_map
from itertools import combinations
import copy


"""
Network Analysis
"""

def s_core(G, s, genes, corrmat, pmat):
    """ Returns the s-core of network G.

    Parameters
    ----------
    G : NetworkX graph
        A graph object of the yeast PIN. Can either be a previously generated s-core or
        the full yeast PIN.
    s : float
        Strength cutoff for node deletion.
    genes : list
        List of genes (same order as rows and columns
        of corrmat and pmat).
    corrmat : Numpy array
        Correlation matrix
    pmat : Numpy array
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
    n = 0   # number of nodes after node removal

    # Prune network
    while n_old != n:
        nodes = s_core.nodes()
        n_old = s_core.number_of_nodes()

        # Remove node if sum of strengths is smaller than cutoff s
        del_nodes = []
        for node in nodes:
            edges = list(s_core.edges(node))
            strengths = sum([s_core[edge[0]][edge[1]]["strength"] for edge in edges])
            if strengths < s:
                del_nodes.append(node)
        s_core.remove_nodes_from(del_nodes)
        n = s_core.number_of_nodes()

    return s_core


def faster_corrcoef(matrix):
    """ Returns the Pearson correlation and p-value matrix
    between all pairs of rows in matrix.

    Adapted from the following Stack Overflow post: 
    https://stackoverflow.com/questions/24432101/correlation-coefficients-and-p-values-for-all-pairs-of-rows-of-a-matrix

    Parameters
    ----------
    matrix : Numpy array
        Matrix of numerical gene expression values from normalized_expression.tsv. Each
        row and column corresponds to a gene and an experimental condition, respectively.

    Returns
    -------
    corrmat : Numpy array
        Correlation matrix
    pmat : Numpy array
        P-value matrix
    """

    # Calculate row-wise Pearson correlations
    corrmat = np.corrcoef(matrix)

    # Calculate associated p-values
    rf = corrmat[np.triu_indices(corrmat.shape[0], 1)]
    df = matrix.shape[1] - 2
    ts = rf * rf * (df / (1 - rf * rf))
    pf = special.betainc(0.5 * df, 0.5, df / (df + ts))
    pmat = np.zeros(shape=corrmat.shape)
    pmat[np.triu_indices(pmat.shape[0], 1)] = pf
    pmat[np.tril_indices(pmat.shape[0], -1)] = pmat.T[np.tril_indices(pmat.shape[0], -1)]
    pmat[np.diag_indices(pmat.shape[0])] = np.ones(pmat.shape[0])

    return corrmat, pmat


"""
Metabolic Analysis
"""


"""
The following contains an implementation of Phenotype Phase Plane Analysis (PhPP).
(Edwards et al. 2001, Characterizing the metabolic phenotype: A phenotype phase plane analysis)
Adapted from code made by Kai Zhuang from the python package framed, as well as cobrapy v0.5.10

November 2021, Vetle Simensen
"""

class PhenotypePhasePlane(object):
    def __init__(self, rxn_x, rxn_y, rxn_x_range, rxn_y_range):
        self.rxn_x = rxn_x
        self.rxn_y = rxn_y

        # converting reaction ranges to numpy array and storing it inside self
        self.x_range = rxn_x_range
        self.y_range = rxn_y_range

        # find length of reaction ranges
        len_x = len(self.x_range)
        len_y = len(self.y_range)

        # creating empty arrays for storing analysis results
        self.f_objective = np.zeros((len_x, len_y))
        self.shadow_price_x = np.zeros((len_x, len_y))
        self.shadow_price_y = np.zeros((len_x, len_y))
        self.segments = np.zeros(self.f_objective.shape, dtype=np.int32)
        self.phases = []

    def plot_PhPP(self, theme='Dark2', new_figure=True, show_plot=True):
        """ Plots a 3D plot of the PhPP.

        Parameters
        ----------
        theme : str
            Color theme used for distinguishing the different phases.
        new_figure : bool
            If set to True, a new matplotlib figure will be created.
        show_plot : bool
            If set to True, current figure will be shown
        """

        if new_figure:
            fig = plt.figure()
            ax = fig.add_subplot(111, projection="3d")

        f = self.f_objective
        x = self.x_range
        y = self.y_range

        if 'EX_' in self.rxn_x:
            x = x * -1

        if 'EX_' in self.rxn_y:
            y = y * -1

        xgrid, ygrid = np.meshgrid(x, y)
        xgrid = xgrid.transpose()
        ygrid = ygrid.transpose()
        xgrid_scale, ygrid_scale = (1, 1)

        # Define theme colors
        colors = np.empty(self.f_objective.shape, dtype=np.dtype((str, 7)))
        n_segments = self.segments.max()
        color_list = get_map(theme, 'Qualitative', n_segments).hex_colors
        for i in range(n_segments):
            colors[self.segments == (i + 1)] = color_list[i]

        # Make surface plots, and add wireframe and axis labels
        ax.plot_surface(xgrid, ygrid, f, facecolors=colors, rstride=1, cstride=1, linewidth=0, antialiased=False)
        ax.plot_wireframe(xgrid, ygrid, f, color="black", rstride=xgrid_scale, cstride=ygrid_scale)
        ax.set_xlabel(self.rxn_x)
        ax.set_ylabel(self.rxn_y)
        ax.set_zlabel("Growth rate")
        ax.view_init(elev=30, azim=-135)
        fig.set_tight_layout(True)

        if show_plot:
            plt.show()

    def plot_shadow_price_x(self, new_figure=True, show_plot=True):
        """ Plots the shadow price associated with reaction x

        Parameters
        ----------
        new_figure : bool
            If set to True, a new matplotlib figure will be created.
        show_plot : bool
            If set to True, current figure will be shown
        """

        if new_figure:
            fig = plt.figure()
            ax = fig.add_subplot()
        sp_x = self.shadow_price_x
        x = self.x_range
        y = self.y_range

        if 'EX_' in self.rxn_x:
            x = x * -1

        if 'EX_' in self.rxn_y:
            y = y * -1

        plt.pcolormesh(x, y, np.transpose(sp_x))
        plt.colorbar()
        ax.set_xlabel(self.rxn_x)
        ax.set_ylabel(self.rxn_y)

        if show_plot:
            plt.show()

    def plot_shadow_price_y(self, new_figure=True, show_plot=True):
        """ Plots the shadow price associated with reaction y

        Parameters
        ----------
        new_figure : bool
            If set to True, a new matplotlib figure will be created.
        show_plot : bool
            If set to True, current figure will be shown
        """

        if new_figure:
            fig = plt.figure()
            ax = fig.add_subplot()
        sp_y = self.shadow_price_y
        x = self.x_range
        y = self.y_range

        if 'EX_' in self.rxn_x:
            x = x * -1

        if 'EX_' in self.rxn_y:
            y = y * -1

        plt.pcolormesh(x, y, np.transpose(sp_y))
        plt.colorbar()
        ax.set_xlabel(self.rxn_x)
        ax.set_ylabel(self.rxn_y)

        if show_plot:
            plt.show()

    def segment(self, threshold=0.001):
        """ This method attempts to segment the data and identify the various phases
        of the phenotype phase plane.

        Parameters
        ----------
        threshold : float
            Threshold for categorizing a new point as belonging to a different phase.
        """

        self.segments *= 0
        # each entry in phases will consist of the following tuple
        # ((x, y), shadow_price1, shadow_price2)
        self.phases = []
        segment_id = 0

        while self.segments.min() == 0:
            segment_id += 1
            # i and j are indices for a current point which has not been
            # assigned a segment yet
            i, j = np.unravel_index(self.segments.argmin(), self.segments.shape)
            # update the segment id for any point with a similar shadow price
            # to the current point
            d1 = abs(self.shadow_price_x - self.shadow_price_x[i, j])
            d2 = abs(self.shadow_price_y - self.shadow_price_y[i, j])
            self.segments[(d1 < threshold) * (d2 < threshold)] += segment_id
            # add the current point as one of the phases
            self.phases.append(((self.x_range[i], self.y_range[j]), self.shadow_price_x[i, j], self.shadow_price_y[i, j]))

def PhPP(model, rxn_x, rxn_y, rxn_x_range, rxn_y_range):
    """ Phenotype Phase Plane Analysis - analyzes the changes in the objective function and the shadow prices.

    Parameters
    ----------
    model : cobra model
        The constraint-based metabolic model.
    rxn_x : str
        ID of reaction to be plotted along the x-axis.
    rxn_y : str
        ID of reaction to be plotted along the y-axis.
    rxn_x_range : list or array 
        The range of reaction x
    rxn_y_range : list or array 
        The range of reaction y

    Returns
    -------
    phaseplane : PhenotypePhasePlane
        PhenotypePhasePlane object.
    """

    # save initial flux bounds (avoids issues when lb > ub)
    old_lb_x = model.reactions.get_by_id(rxn_x).lower_bound
    old_ub_x = model.reactions.get_by_id(rxn_x).upper_bound
    old_lb_y = model.reactions.get_by_id(rxn_y).lower_bound
    old_ub_y = model.reactions.get_by_id(rxn_y).upper_bound

    # find metabolite ids corresponding to reactions x and y
    met_x = list(model.reactions.get_by_id(rxn_x).metabolites)[0].id
    met_y = list(model.reactions.get_by_id(rxn_y).metabolites)[0].id

    # create a PhenotypePhasePlane instance for storing results
    phase_plane = PhenotypePhasePlane(rxn_x, rxn_y, rxn_x_range, rxn_y_range)

    for i, v_x in enumerate(rxn_x_range):
        model.reactions.get_by_id(rxn_x).lower_bound = v_x
        model.reactions.get_by_id(rxn_x).upper_bound = v_x

        for j, v_y in enumerate(rxn_y_range):
            model.reactions.get_by_id(rxn_y).lower_bound = v_y
            model.reactions.get_by_id(rxn_y).upper_bound = v_y

            try:
                solution = model.optimize()
                if solution.status == 'optimal':
                    phase_plane.f_objective[i, j] = solution.objective_value
                    phase_plane.shadow_price_x[i, j] = solution.shadow_prices[met_x]
                    phase_plane.shadow_price_y[i, j] = solution.shadow_prices[met_y]
            except UserWarning:
                pass

        model.reactions.get_by_id(rxn_y).upper_bound = 0
        model.reactions.get_by_id(rxn_y).lower_bound = 0

    phase_plane.segment()

    # set bounds back to initial values
    model.reactions.get_by_id(rxn_x).upper_bound = old_ub_x
    model.reactions.get_by_id(rxn_x).lower_bound = old_lb_x
    model.reactions.get_by_id(rxn_y).upper_bound = old_ub_y
    model.reactions.get_by_id(rxn_y).lower_bound = old_lb_y

    return phase_plane


"""
The following is a brute-force implementation of OptKnock as the OptKnock functionality found
in the Cameo package is buggy and rather annoying.
"""

def phonyOptKnock(model, target, reactions=None, knockouts=1):
    """ Phony implementation of OptKnock to identify mutants 
    with growth-coupled target metabolite production.

    Parameters
    ----------
    model : cobra model object
            Model object.
    
    target : string
             Reaction id of target reaction to be improved through model knockouts.

    reactions : cobra reaction object(s) 
                Reactions to be considered for knockout. By default
                all model reactions.

    knockouts : int
                Number of knockouts to consider, single (1),
                double (2), triple (3), etc.

    Returns
    ----------
    results : pandas DataFrame object
              Results for all proposed mutants containing the columns:

              reactions: knocked-out reaction(s)

              target: flux of target reaction when optimizing the model objective function

              biomass: corresponding optimal objective value of the model (i.e. growth if biomass is the objective)

              fva_min: minimal flux through the target reaction at optimal growth

              fva_max: maximal flux through the target reaction at optimal growth
    """

    del_model = copy.deepcopy(model)

    # By default all model reactions
    if reactions == None:
        reactions = del_model.reactions

    # Remove biomass, exchange, and reactions without associated genes
    reactions.remove(list(cobra.util.solver.linear_reaction_coefficients(del_model).keys())[0])
    reactions = [rxn for rxn in reactions if rxn not in cobra.medium.find_boundary_types(del_model, "exchange")]
    reactions = [rxn for rxn in reactions if not "s0001" in rxn.gene_reaction_rule]

    # All unique knockout combinations
    if knockouts > 1:
        candidates = [i for i in combinations(reactions, knockouts)]
    else:
        candidates = list(reactions)
    
    result = pd.DataFrame(index=range(len(candidates)),columns=["reactions", "target", "biomass", "fva_min", "fva_max"])

    sol_wt = del_model.optimize()

    # Constrain and simulate fluxes
    for i, cand in enumerate(candidates):
        skip = False
        if knockouts == 1:
            cand = [cand]
        else:
            cand = list(cand)

        old_lower_bounds = np.zeros(knockouts)
        old_upper_bounds = np.zeros(knockouts)
        for j, c in enumerate(cand):
            old_lower_bounds[j] = cand[j].lower_bound
            old_upper_bounds[j] = cand[j].upper_bound
            c.lower_bound = 0.0
            c.upper_bound = 0.0

        # Simulate fluxes
        sol = del_model.optimize()

        # Skip if growth is below 10% of wild type
        if sol.objective_value < 0.1 * sol_wt.objective_value:
            skip = True

        # Simulate flux variability, skip if infeasible
        try:
            sol_fva = cobra.flux_analysis.flux_variability_analysis(del_model, target)
        except:
            skip = True

        if not skip:
            # Add results to dataframe
            result["reactions"][i] = [c.id for c in cand]
            result["target"][i] = sol.fluxes[target]
            result["biomass"][i] = sol.objective_value
            result["fva_min"][i] = sol_fva["minimum"][0]
            if sol_fva["maximum"][0] < 1e-10:
                result["fva_max"][i] = 0
            else:
                result["fva_max"][i] = sol_fva["maximum"][0]
                

        # Reset previous flux bounds
        for j, rxn in enumerate(cand):
            rxn.upper_bound = old_upper_bounds[j]
            rxn.lower_bound = old_lower_bounds[j]
    
    # Post-process results - remove nan rows and sort
    result = result.dropna()
    result = result.sort_values(by=['fva_min'], ascending=[False])
    return result
