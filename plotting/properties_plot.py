import pandas
import seaborn as sns
import networkx as nx
import matplotlib.pyplot as plt

def properties_per_layer(G):

layer = nx.get_node_attributes(G, 'layer')
volume = nx.get_node_attributes(G, 'volume')
area = nx.get_node_attributes(G, 'area')
RF = nx.get_node_attributes(G, 'Roundness_factor')
volume_list = [0 for i in range(G.number_of_nodes())]
area_list = [0 for i in range(G.number_of_nodes())]
RF_list = [0 for i in range(G.number_of_nodes())]
degree_list = [0 for i in range(G.graph['nb_of_layer'] -1)]
nbcells_per_layer = [0 for i in range(G.graph['nb_of_layer'] -1)]
layer_list = [0 for i in range(G.number_of_nodes())]
degree_list = [0 for i in range(G.number_of_nodes())]
nodes_list = [i for i in range(G.number_of_nodes())]
layer_list2 = list(range(G.graph['nb_of_layer'] -1))
for ind in nodes_list:
  if ind in volume:
    volume_list[ind] = volume[ind]
    area_list[ind] = area[ind]
    RF_list[ind] = RF[ind]
    layer_list[ind] = layer[ind]
    degree_list[ind] = G.degree(ind)
    nbcells_per_layer[layer[ind]] = nbcells_per_layer[layer[ind]] + 1
spheroid_pd = pandas.DataFrame({"Nodes": nodes_list, "Layer": layer_list, "Volume" : volume_list, "Degree" : degree_list, "Area" : area_list, "Roundness Factor" : RF_list })
#degree_pd = pandas.DataFrame({"Nodes": nodes_list, "Layer": layer_list, "Degree" : degree_list})
plt.figure()
ax = sns.violinplot(x = "Layer", y = "Volume", data = spheroid_pd, cut=0, palette = "Set2", scale="count", inner="points")
plt.figure()
ax = sns.violinplot(x = "Layer", y = "Area", data = spheroid_pd, cut=0, palette="Set2", scale="count", inner="points")
plt.figure()
ax = sns.violinplot(x = "Layer", y = "Degree", data = spheroid_pd, cut=0, palette="Set2", scale="count", inner="points")
plt.figure()
ax = sns.violinplot(x = "Layer", y = "Roundness Factor", data = spheroid_pd, cut=0, palette="Set2", scale="count", inner="points")
