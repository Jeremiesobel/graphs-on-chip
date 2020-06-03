import numpy as np
import pandas
import networkx as nx
import math

def density_and_mean_degree(G):

"""
From a networkx graph, return:
- the density, ie the number of edges over the maximum number of edges given the number of nb_nodes
- the mean degree, ie the number of neighbours of a node on average
"""
  nb_nodes = len(G)

  nb_edges = G.number_of_edges()

  nb_edges_max = (nb_nodes**2 + nb_nodes)/2

  return nb_edges/nb_edges_max, 2*nb_edges/nb_nodes

def Energy_cells(G, Egg=1, Egr=1.5, Err=1):
"""
From a networkx graph, return an energy calculated from coefficients, and independantly from the surface of contact
"""
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
    """
    From a networkx graph, return an energy calculated from coefficients, taking into account the distance between cells
    """
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

 def Energy_cells_voronoi(G, Egg=1, Egr=1.5, Err=1, include_color:bool = False, V0 = (dCells/2)**3, S0 = (dCells/2)**2, Ka=1, Ks=1):
"""
From a networkx graph, return an energy calculated from coefficients, using the same Hamiltonian as Lisa Manning, and with cross terms taking into account the surface of contact between cells
"""
  E=0
  colors = nx.get_node_attributes(G, 'color')
  edge_surfaces  = nx.get_edge_attributes(G, 'Area')
  cell_area  = nx.get_node_attributes(G, 'area')
  cell_volume  = nx.get_node_attributes(G, 'volume')
  for ind in list(G.nodes):
    if ind in cell_area:

      E = E + Ks*(cell_area[ind] - S0) + Ka*(cell_volume[ind] - V0)

  if include_color:

      assert If is not None

      for edge in list(G.edges()):
        if edge in edge_surfaces:

          if colors[edge[0]] == 'g' and colors[edge[1]]  == 'g':
            E = E + Egg*edge_surfaces[edge]

          elif (colors[edge[0]]  == 'g' and colors[edge[1]]  == 'r') or (colors[edge[0]]  == 'r' and colors[edge[1]]  == 'g'):
            E = E + Egr*edge_surfaces[edge]

          elif colors[edge[0]]  == 'r' and colors[edge[1]]  == 'r':
            E = E + Err*edge_surfaces[edge]

  return E

def mean_lengths(G):
    """
    From a networkx graph, return the average distance between cell core, depending on the color of each cells
    """
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




def surface_of_facet(facet, nodes, dsph):
    """
    From the facet between two cells, calculate the surface
    facet: an array of the indexes of the vertexes composing the polygone between the two cells_ind
    nodes: an array of the positions of each vertex of the voronoi vertices
    dsph: the radius of the spheroid considering that it is a sphere
    """

    m = len(facet)
    if m > 2:
        facet_pos = np.zeros((m, 2))
        #creating the base of the flat to project the points
        v1 = (nodes[facet[1]] - nodes[facet[0]])/distance.euclidean(nodes[facet[1]], nodes[facet[0]])
        v2 = (nodes[facet[2]] - nodes[facet[0]])
        normal = np.cross(v1, v2)
        v3 = np.cross(v1, normal)
        v3 = v3/np.linalg.norm(v3)
        for i in range(m):

            facet_pos[i][0] = np.dot(nodes[facet[i]], v1)
            facet_pos[i][1] = np.dot(nodes[facet[i]], v3)

        Surface_of_contact = ConvexHull(facet_pos).area
    return Surface_of_contact

def get_volume_and_surfaces2(G):

    """
    From a networkx graph,
    - adding to each edge the surface as an attribute
    - adding to each node the volumen area and roundness factor as an attribute
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
    dsph = (3*vol/(4*math.pi))**(1/3)  #the radius of the spheroid considering that it is a sphere

    for i in range(len(nodes)):
        dst = distance.euclidean(nodes[i], center)
        if dst > dsph:
            nodes[i] = (nodes[i] - center)/distance.euclidean(nodes[i], center)*dsph

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
          hull_cell = ConvexHull(L)
          """else:
            L[i] = (nodes[ind] - center)/distance.euclidean(nodes[ind], center)*dsph"""



          vol_cell = hull_cell.volume
          surface_cell = hull_cell.area
          RF = 36*math.pi*vol_cell**2/surface_cell**3
          G.add_node(l, pos_ridge_vertices = L)
          G.add_node(l, hull = hull)
          G.add_node(l, equations = hull.equations)
          G.add_node(l, volume = vol_cell)
          G.add_node(l, area = surface_cell)
          G.add_node(l, Roundness_factor = RF)
        l = l+1
      ridge_vertices = nx.get_node_attributes(G,'ridge_vertices')
      pos_ridge_vertices = nx.get_node_attributes(G,'pos_ridge_vertices')
      hull = nx.get_node_attributes(G, 'hull')
      for edge in list(G.edges()):

        #pos_common_facet = list(set(pos_ridge_vertices[edge[0]]).intersection(pos_ridge_vertices[edge[1]]))
        common_facet = list(set(ridge_vertices[edge[0]]).intersection(ridge_vertices[edge[1]]))
        if -1 in common_facet:
            common_facet = list(set(common_facet)).remove(-1)
        elif not(common_facet):
            None
        else:

            G[edge[0]][edge[1]]['Area'] = surface_of_facet(common_facet, nodes, dsph)
      return G

def attribute_layer(G):
    """
    From a graph:
    - construct the convex hull to determine the outter layer: layer = 0 as an attribute
    - construct the convex hull of the spheroid without the outter cells to determine the first inner layer: layer = 1
    - and so on until there is not enough cells to construct the convex hull
    - attribute the last layer to the remaining cells
    """
  npoints = G.number_of_nodes()
  pos = nx.get_node_attributes(G, 'pos')

  points = np.zeros((npoints, 3))
  layer = 0
  for ind in range(npoints):
    points[ind][0] = pos[ind][0]
    points[ind][1] = pos[ind][1]
    points[ind][2] = pos[ind][2]
  pts = points.tolist()
  vor = Voronoi(points)
  hull = ConvexHull(points)
  L = []
  while len(points) > 3 :
    hull = ConvexHull(points)
    points = points.tolist()
    k = len(points)
    for l  in range(1, k+1):
      if k-l in hull.vertices:
        n = pts.index(points.pop(k-l))
        G.add_node(n, layer = layer)
        L.append(n)
    layer = layer +1
    points = np.asarray(points)
  if len(G) - len(L) >0:
    for l in range(len(G)):
      if l not in L:
        G.add_node(l, layer = layer)
        L.append(l)
  return G
