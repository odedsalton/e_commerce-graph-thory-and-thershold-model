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
        """Expects data set file with index column (train and test) """
        self.init_graph = pd.read_csv(dataset)
        self.nodes = self.init_graph[['source', 'target']]
def competitive_part(links_dataset , number):

    return [[5,9],[7,6]]
def write_file_competition (list_of_removed):
    filename = "hw1_competition.csv"
    with open(filename, 'w') as csv_file:
        writer = csv.writer(csv_file, lineterminator='\n')
        fieldnames2 = ["source", "target"]
        writer.writerow(fieldnames2)
        for edge in list_of_removed:
            writer.writerow([edge[0], edge[1]])

