#############################################
#   This is a module created to hold        #
#   all custom functions used in multiple   #
#   project throughout the TBT4165 course   #
#############################################


## Imports:
import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
from pyvis.network import Network


## Functions:

def connected_ER(N, E):
    """
    N : int, Number of nodes
    E : int, Number of edges 
    """
    g = nx.gnm_random_graph(N, E)
    while nx.is_connected(g) == False:
        g = nx.gnm_random_graph(N, E)
    return g


def connected_BA(N, m):
    """
    N : int, Number of nodes
    m : int, Number of edges to attach from a new node to existing nodes
    """
    g = nx.barabasi_albert_graph(n=N, m=m)
    while nx.is_connected(g) == False:
        g = nx.barabasi_albert_graph(n=N, m=m)
    return g


def get_largest_component(g, is_undirected=True):
    if is_undirected:
        l_comp = g.subgraph(
            sorted(
                nx.connected_components(
                    nx.to_undirected(g.copy())
                    ),
                reverse=True
            )[0]
        ).copy()
    else:
        l_comp = g.subgraph(
            sorted(
                nx.connected_components(g),
                reverse=True
            )[0]
        ).copy()
    return l_comp


def graph_info(g):
    print(f"\nNumber of nodes: {g.number_of_nodes()}")
    print(f"Number of edges: {g.number_of_edges()}")
    try:
        print(f"Is connected:    {nx.is_connected(g)}")
    except:
        print("""\
The graph is not undirected. Therefore .is_connected() does not work.
        """)


def show_html(graph, name="nx", show=True, size="small"):
    """
    Generate and display the graph through pyvis
    """
    
    if size == "small":
        nt = Network("500px", "500px")
    else:
        nt = Network("1080px", "1920px")
        
    nt.from_nx(graph)
    if show == True:
        nt.show(f"{name}.html")



def nice_plot(title, xlab, ylab, show=True):
    """
    Function to generate title and axis labels to a plot.
    The argument "show" is by default True. 
    """
    plt.title(title)
    plt.xlabel(xlab)
    plt.ylabel(ylab)
    if (show == True):
        plt.show()


def plot_degree_distribution(graph):
    """
    Function that returns a scatter plot of the degree
    distribution of a graph. 
    """
    N = graph.number_of_nodes()  # Total number of nodes
    y = np.array(nx.degree_histogram(graph))/N  # All occurrences of "n" degree, divided by total number of nodes
    x = np.arange(len(y))[y != 0]  # x values
    y = y[y != 0]  # Remove the values = 0 from the array
    
    return [x, y]
    #return plt.scatter(x,y)