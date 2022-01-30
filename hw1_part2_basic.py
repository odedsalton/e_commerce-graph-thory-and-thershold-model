import numpy as np
import pandas as pd
import networkx as nx
import random
import csv
from scipy.special import comb
import matplotlib.pyplot as plt
random.seed(0)

class ModelData:
    """The class reads 5 files as specified in the init function, it creates basic containers for the data.
    See the get functions for the possible options, it also creates and stores a unique index for each user and movie
    """

    def __init__(self, dataset):
        self.init_graph = pd.read_csv(dataset)
        self.nodes = self.init_graph[['source', 'target']]
def removing_highest_weight(links_dataset , number):
    G = nx.from_pandas_edgelist(links_dataset.init_graph, 'target', 'source', edge_attr=True, create_using=nx.Graph())
    weight_dict = dict()
    for edge in G.edges():
        node1 = edge[0]
        node2 = edge[1]
        weight_dict[node1,node2] = G[edge[0]][edge[1]]['weight']
    weight_orderd_dict = dict(sorted(weight_dict.items(), key=lambda item: item[1], reverse=True))
    iter = 0
    list_of_removed = list()
    for key1,key2 in weight_orderd_dict.keys():
        if iter == number:
            break
        iter = iter + 1
        G.remove_edge(key1,key2)
        list_of_removed.append([key1,key2])
    return list_of_removed