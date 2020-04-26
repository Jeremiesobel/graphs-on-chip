import numpy as np
import pandas
import os
import io
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import networkx as nx
import scipy as sp
import scipy.spatial as sptl
from scipy.spatial import KDTree
from scipy.spatial import Voronoi
from scipy.spatial import ConvexHull

def generate_geometric_graph(spheroid:dict, 
                            zRatio:float, 
                            dCells:int):
  
    """
  
    From an input dict generates networkx geometric graph.

    Returns:

    - networkx graph
  
    """
    cells = spheroid['cells']
    # Generate a dict of positions
    pos = {int(i): (int(cells[i]['x']), int(cells[i]['y']), int(cells[i]['z']/zRatio)) for i in cells.keys()}
    
    # Create random 3D network
    G = nx.random_geometric_graph(len(cells), dCells, pos=pos)

    if 'intensity' in spheroid.keys():

        for ind in list(G.nodes):

            assert 'color' in cells[ind].keys()

            G.add_node(ind, color = cells[ind]['color'])

    return G