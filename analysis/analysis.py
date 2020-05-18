import numpy as np
import pandas
import networkx as nx
import math

def density_and_mean_degree(G):

  nb_nodes = len(G)

  nb_edges = G.number_of_edges()

  nb_edges_max = (nb_nodes**2 + nb_nodes)/2

  return nb_edges/nb_edges_max, 2*nb_edges/nb_nodes

def Energy_cells(G, Egg=1, Egr=1.5, Err=1):

  E=0
  colors = nx.get_node_attributes(G, 'color')
  neighbours = nx.get_node_attributes(G, 'neighbour')

  for edge in list(G.edges()):

    if colors[edge[0]] == 'g' and colors[edge[1]]  == 'g':
      E = E + Egg

    elif (colors[edge[0]]  == 'g' and colors[edge[1]]  == 'r') or (colors[edge[0]]  == 'r' and colors[edge[1]]  == 'g'):
      E = E + Egr

    elif colors[edge[0]]  == 'r' and colors[edge[1]]  == 'r':
      E = E + Err

  return E

def Energy_cells_distance(G, Egg=1, Egr=1.5, Err=1):
  E=0
  colors = nx.get_node_attributes(G, 'color')
  neighbours = nx.get_node_attributes(G, 'neighbour')
  lengths = nx.get_edge_attributes(G, 'length')

  for edge in list(G.edges()):

    if colors[edge[0]] == 'g' and colors[edge[1]]  == 'g':
      E = E + Egg*((dCells)**2-lengths[edge]**2)/4

    elif (colors[edge[0]]  == 'g' and colors[edge[1]]  == 'r') or (colors[edge[0]]  == 'r' and colors[edge[1]]  == 'g'):
      E = E + Egr*((dCells)**2-lengths[edge]**2)/4

    elif colors[edge[0]]  == 'r' and colors[edge[1]]  == 'r':
      E = E + Err*((dCells)**2-lengths[edge]**2)/4

  return E

def mean_lengths(G):
  total_length = 0
  total_length_gg = 0
  total_length_rr = 0
  total_length_gr = 0
  nb_gglinks = 0
  nb_rrlinks = 0
  nb_grlinks = 0
  lengths = nx.get_edge_attributes(G, 'length')
  colors = nx.get_node_attributes(G, 'color')


  for edge in list(G.edges()):
    total_length = total_length + lengths[edge]
    if colors[edge[0]] == 'g' and colors[edge[1]]  == 'g':
      total_length_gg = total_length_gg + lengths[edge]
      nb_gglinks = nb_gglinks +1
    elif (colors[edge[0]]  == 'g' and colors[edge[1]]  == 'r') or (colors[edge[0]]  == 'r' and colors[edge[1]]  == 'g'):
     total_length_gr = total_length_gr + lengths[edge]
     nb_grlinks = nb_grlinks +1
    elif colors[edge[0]] == 'r' and colors[edge[1]]  == 'r':
      total_length_rr = total_length_rr + lengths[edge]
      nb_rrlinks = nb_rrlinks +1


  return total_length/len(G.edges), total_length_gg/nb_gglinks, total_length_rr/nb_rrlinks, total_length_gr/nb_grlinks
