import math
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import networkx as nx
import random
from scipy.special import comb
random.seed(0)

""" The return values in all functions are only for demonstrating the format, and should be overwritten """

class ModelData:
    def __init__(self, dataset):
        self.init_graph = pd.read_csv(dataset)
        self.nodes = self.init_graph[['source', 'target']]

def graphA(links_data):
    # Implement for Q1
    G = nx.from_pandas_edgelist(links_data.init_graph, 'source', 'target', edge_attr=True, create_using=nx.Graph())
    return G

def calc_best_power_law(G):
    # Implement for Q1
    out_degree = dict(G.degree())
    count_number_nodes = {}
    elements_degree = list(out_degree.values())
    m = max(elements_degree)
    n = 1
    for i in range(m):
        for j in range(len(elements_degree)):
            if elements_degree[j]==i :
                count_number_nodes.update({i: n})
                n = n + 1
        n = 1
    X = list(count_number_nodes.keys())
    Y = list(count_number_nodes.values())
    sum1 = len(Y)
    for i in range(sum1):
        Y[i] = math.log(Y[i],10)
        if X[i] == 0:
            X.remove(X[i])
            Y.remove(Y[i])
        else:
            X[i] = math.log(X[i],10)
    print(len(Y))
    plt.scatter(X,Y,label ='log(count)')
    beta, alpha = np.polyfit(X, Y,1)
    # p = np.poly1d((beta,alpha))
    return alpha , beta
def plot_histogram(G, alpha, beta):
    p = np.poly1d((beta, alpha))
    plt.plot([0,2.7],[p(0),p(2.7)],color= 'red', label='power law={0:.2f}{1:2f}*log.(out-degree)'.format(alpha,beta))
    plt.xlabel('log(out-degree)[base 10]')
    plt.ylabel('log(number_of_nodes)[base10]')
    plt.show()
    return
def G_features(G): #########return dict of dicts ########
    # Implement for Q2
    SLCC = max(nx.connected_components(G), key=len)
    list_w = list(SLCC)
    list_n = list(G.nodes())
    for node in list_n:
        if node not in list_w:
            G.remove_node(node)
    n = len(list_w)
    a = {}
    x = 0
    for i in list_w:
        x = x + 1
        dist = nx.shortest_path_length(G, i, target=None, weight=None)
        if ((sum(dist.values())) != 0):
            a.update({i: (n - 1) / (sum(dist.values()))})
        else:
            a.update({i: 0})
    b = betweenness_centrality(G)
    nested_dict = {'a': a, 'b': b}
    return nested_dict

def calc_flow_val(betweenness, S, P, sigma, s):
    delta = dict.fromkeys(S, 0)
    while S:
        w = S.pop()
        flow_val = (1.0 + delta[w]) / sigma[w]
        for vv in P[w]:
            delta[vv] += sigma[vv] * flow_val
        if w != s:
            betweenness[w] += delta[w]

    return betweenness
def bfs(G, s):  # bfs for all shortes path
    S = []
    P = {}
    for v in G:
        P[v] = []
    sigma = dict.fromkeys(G, 0.0)
    D = {}
    sigma[s] = 1.0
    D[s] = 0
    Q = [s]
    x = 0
    while Q:  # while Q is not empty
        x = x + 1
        vi = Q.pop(0)
        S.append(vi)
        Dv = D[vi]
        sigmav = sigma[vi]
        for w in G[vi]:
            if w not in D:
                Q.append(w)
                D[w] = Dv + 1
            if D[w] == Dv + 1:  # relax
                sigma[w] += sigmav
                P[w].append(vi)  # predecessors
    return S, P, sigma
def betweenness_centrality(G):
    betweenness = dict.fromkeys(G, 0.0)
    for s in G.nodes:
        S, P, sigma = bfs(G, s)
        betweenness = calc_flow_val(betweenness, S, P, sigma, s)

    n = G.number_of_nodes()
    nCr = math.factorial(n - 1) / (math.factorial(2) * math.factorial(n - 3))
    multiplier = (1 / nCr)

    for bet in betweenness:
        betweenness[bet] = betweenness[bet] / 2 * multiplier

    return betweenness
def create_undirected_weighted_graph(links_data , users_data, question): #########return Graph ########
    # Implement for Q3
    G = nx.from_pandas_edgelist(links_data.init_graph, 'target', 'source', edge_attr=True, create_using=nx.Graph())
    nodes = list(users_data['node'])
    is_infected = list(users_data[question])
    nodes_dict = dict()
    for i in range(len(nodes)):
        nodes_dict[nodes[i]] = is_infected[i]
    nx.set_node_attributes(G,nodes_dict,'is_infected')
    return G
def run_k_iterations(WG, k , Threshold):
    # Implement for Q3
    S = dict()
    dict_nodes_value = dict()
    iter_list_nodes = list()
    for node in WG.nodes():
        if WG.nodes[node]['is_infected'] == 'YES':
            iter_list_nodes.append(node)
    S[0] = iter_list_nodes
    for i in range(k):
        iter_list_nodes = list()
        list_n = list()
        for node in WG.nodes():
            if WG.nodes[node]['is_infected'] == 'YES':
                list_n.append(node)
            else:
                dict_nodes_value[node] = 0
        for node in list_n:
            for node2 in WG.neighbors(node):
                if WG.nodes[node2]['is_infected'] == 'YES':
                    continue
                else:
                    dict_nodes_value[node2] = dict_nodes_value[node2] + WG[node][node2]['weight']
                    if dict_nodes_value[node2] > Threshold:
                        WG.nodes[node2]['is_infected'] = 'YES'
                        iter_list_nodes.append(node2)
        S[i + 1] = iter_list_nodes
    return S
# return a dictionary of with the branching factors R1,...Rk
def calculate_branching_factor(S,k):
    # Implement for Q3
    branching_factor = dict()
    for i in range(k - 1):
        list_l = S[i]
        list_m = S[i+1]
        branching_factor[i+1] = len(list_m)/len(list_l)
    return branching_factor
# return a dictionary of the h nodes with the highest degree, sorted by decreasing degree
# {index : node}
def find_maximal_h_deg(WG,h):
    # implement for Q4
    dict_deg = dict(WG.degree())
    x = sorted(dict_deg.items(), key=lambda item: item[1], reverse=True)
    sum = len(list(WG.nodes()))
    while sum > h:
        x.pop()
        sum = sum -1
    nodes_dict = dict(x)
    return nodes_dict
# return a dictionary of all nodes with their clustering coefficient
# {node : CC}
def calculate_clustering_coefficient(WG,nodes_dict):
    # implement for Q4
    nodes_list = list(nodes_dict.keys())
    nodes_dict_clustering = dict()
    for node in nodes_list:
        number_triplrts = 0
        node_negihbors_list = list(WG.neighbors(node))
        if (len(node_negihbors_list)<2):
            nodes_dict_clustering[node] = 0
            continue
        else:
            combinations = comb(len(node_negihbors_list), 2)
        for node2 in node_negihbors_list:
            node2_negihbors_list = WG.neighbors(node2)
            for node3 in node2_negihbors_list:
                if node3 in node_negihbors_list:
                    number_triplrts = number_triplrts+1
                    continue
        nodes_dict_clustering[node] = (number_triplrts / combinations)/2
    return nodes_dict_clustering
def infected_nodes_after_k_iterations(WG, k, Threshold):
    # implement for Q4
    final_infectd_nodes = 0
    dict_nodes_value = dict()
    for i in range(k):
        iter_list_nodes = list()
        list_n = list()
        for node in WG.nodes():
            if WG.nodes[node]['is_infected'] == 'YES':
                list_n.append(node)
            else:
                dict_nodes_value[node] = 0
        for node in list_n:
            for node2 in WG.neighbors(node):
                if WG.nodes[node2]['is_infected'] == 'YES':
                    continue
                else:
                    dict_nodes_value[node2] = dict_nodes_value[node2] + WG[node][node2]['weight']
                    if dict_nodes_value[node2] > Threshold:
                        final_infectd_nodes = final_infectd_nodes + 1
                        WG.nodes[node2]['is_infected'] = 'YES'
                        iter_list_nodes.append(node2)
    return final_infectd_nodes
# return the first [number] nodes in the list
def slice_dict(dict_nodes,number):
    # implement for Q4
    nodes_Tot_list = list(dict_nodes.keys())
    nodes_list = list()
    for i in range(number):
        nodes_list.append(nodes_Tot_list[i])
    return nodes_list
# remove all nodes in [nodes_list] from the graph WG, with their edges, and return the new graph
def nodes_removal(WG,nodes_list):
    # implement for Q4
    for node in nodes_list:
        WG.remove_node(node)
    return WG
# plot the graph according to Q4 , add the histogram to the pdf and run the program without it
def graphB(number_nodes_1,number_nodes__2,number_nodes_3):
    # implement for Q4
    plt.plot(number_nodes_1.keys(), number_nodes_1.values(), 'r--', number_nodes__2.keys(), number_nodes__2.values(), 'bs', number_nodes_3.keys(), number_nodes_3.values() , 'g^')
    plt.ylabel('user infected')
    plt.xlabel('user in_quarntine')
    plt.show()
    return