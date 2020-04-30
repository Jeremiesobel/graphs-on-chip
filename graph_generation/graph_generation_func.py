import numpy as np
import pandas
import os
import io
import networkx as nx
from scipy.spatial import Delaunay

from collections import defaultdict


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
    for edge in list(G.edges()):

      G[edge[0]][edge[1]]['length'] = np.sqrt((cells[edge[0]]['x'] - cells[edge[1]]['x'])**2 + (cells[edge[0]]['y'] - cells[edge[1]]['y'])**2 + (cells[edge[0]]['z'] -cells[edge[1]]['z'])**2/zRatio**2)
 

    return G


def prep_points(cells:dict):

    return [[cells[cell]['x'], cells[cell]['y'], cells[cell]['z']] for cell in cells.keys()]


def find_neighbors(tess):
    neighbors = defaultdict(set)

    for simplex in tess.simplices:
        for idx in simplex:

            other = set(simplex)
            other.remove(idx)
            neighbors[idx] = neighbors[idx].union(other)

    return neighbors

def generate_voronoi_graph(spheroid:dict,
                            zRatio:float,
                            dCells:float = 60):

    """

    From an input dict generates voronoi graph.

    Returns:

    - networkx graph

    """

    cells = spheroid['cells']
    cells_pos = prep_points(cells)
    tri = Delaunay(cells_pos)

    neighbors = find_neighbors(tri)

    G=nx.Graph()
    neighbors = dict(neighbors)

    for key in neighbors:
        for node in neighbors[key]:

            G.add_edge(key, node)

    pos = {int(i): (cells[i]['x'], cells[i]['y'], cells[i]['z']) for i in cells.keys()}
    nx.set_node_attributes(G, pos, 'pos')

    if 'intensity' in spheroid.keys():

        for ind in list(G.nodes):

            assert 'color' in cells[ind].keys()

            cells_ind = list(cells.keys())[ind]

            G.add_node(ind, color = cells[cells_ind]['color'])
        for edge in list(G.edges()):

            G[edge[0]][edge[1]]['length'] = np.sqrt((cells[edge[0]]['x'] - cells[edge[1]]['x'])**2 + (cells[edge[0]]['y'] - cells[edge[1]]['y'])**2 + (cells[edge[0]]['z'] -ells[edge[1]]['z'])**2/zRatio**2)


    return trim_graph_voronoi(G, zRatio, dCells)

def trim_graph_voronoi(G, zRatio, dCells):

    pos =nx.get_node_attributes(G,'pos')
    edges = G.edges()
    to_remove = []

    for e in edges:

        i, j = e
        dx2 = (pos[i][0]-pos[j][0])**2
        dy2 = (pos[i][1]-pos[j][1])**2
        dz2 = (pos[i][2]-pos[j][2])**2

        if  np.sqrt(dx2+dy2+dz2/zRatio**2) > dCells:

            to_remove.append((i,j))

    for e in to_remove:

        i, j = e
        G.remove_edge(i,j)

    return G
