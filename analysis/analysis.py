import numpy as np
import pandas
import networkx
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

def surface_of_facet(facet, nodes):
    """
    From the facet between two cells, calculate the surface
    """
    hull = ConvexHull(points)
    center = points.mean(axis=0)
    vol = hull.volume
    dsph = (vol/4.19)**(1/3)
    if -1 in facet:
        facet = list(set(facet)).remove(-1)
    elif not(facet):
         None
    else:
        m = len(common_facet)
        if m > 2:
          facet_pos = np.zeros((m, 2))
          #creating the base of the flat to project the points
          v1 = nodes[facet[1]] - nodes[facet[0]]
          v2 = nodes[facet[2]] - nodes[facet[0]]
          normal = np.cross(v1, v2)
          v3 = np.cross(v1, normal)
          for i in range(m):
            dst = distance.euclidean(nodes[common_facet[i]], center)
            if dst<dsph:

                common_facet_pos[i][0] = np.dot(nodes[common_facet[i]], v1)
                common_facet_pos[i][1] = np.dot(nodes[common_facet[i]], v3)

            Surface_of_contact = ConvexHull(common_facet_pos).area
            return Surface_of_contact

def get_volume_and_surfaces(G):
    """
    From a graph, adding to each edge the surface as an attribute
    """
  npoints = G.number_of_nodes()
  pos = nx.get_node_attributes(G, 'pos')

  points = np.zeros((npoints, 3))

  for ind in range(npoints):
    points[ind][0] = pos[ind][0]
    points[ind][1] = pos[ind][1]
    points[ind][2] = pos[ind][2]
  vor = Voronoi(points)
  nodes = vor.vertices
  regions = vor.regions
  hull = ConvexHull(points)
  center = points.mean(axis=0)
  vol = hull.volume
  dsph = (3*vol/(4*pi))**(1/3)
  if len(regions) == len(G)+1:
    l = 0
    for ind in vor.point_region:
      G.add_node(l, ridge_vertices = regions[ind])
      n = len(regions[ind])
      L = np.zeros((n, 3))
      regions[ind] = np.asarray(regions[ind])
      if (np.any(regions[ind] < 0) and n > 4) or (np.all(regions[ind] > 0) and n>3):
        for i in range(n):
          if regions[ind][i] != -1:
            L[i][0] = nodes[regions[ind][i]][0]
            L[i][1] = nodes[regions[ind][i]][1]
            L[i][2] = nodes[regions[ind][i]][2]
        hull = ConvexHull(L)
        """else:
          L[i] = (nodes[ind] - center)/distance.euclidean(nodes[ind], center)*dsph"""


      if n !=0:
        vol_cell = hull_cell.volume
        surface_cell = hull_cell.area
        RF = 36*pi*vol_cell**2/surface_cell**3
        G.add_node(l, pos_ridge_vertices = L)
        G.add_node(l, hull = hull)
        G.add_node(l, equations = hull.equations)
        G.add_node(l, volume = vol_cell)
        G.add_node(l, area = surface_cell)
        G.add_node(l, Roundness_factor = RF)
      l = l+1
    ridge_vertices = nx.get_node_attributes(G,'ridge_vertices')
    pos_ridge_vertices = nx.get_node_attributes(G,'pos_ridge_vertices')
    equations = nx.get_node_attributes(G,'equations')
    hull = nx.get_node_attributes(G, 'hull')
    for edge in list(G.edges()):

      common_facet = list(set(ridge_vertices[edge[0]]).intersection(ridge_vertices[edge[1]]))
      G[edge[0]][edge[1]]['Area'] = surface_of_facet(common_facet, nodes)
