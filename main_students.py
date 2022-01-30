import hw1_part1 as st
import pandas as pd
import matplotlib.pyplot as plt
def main():
    ###########################################################################
    #######################            Part A           #######################
    ###########################################################################
    """Read the data files with the class ModelData into the a data object"""
    file_name = 'result_part_A.csv'
    links_data = st.ModelData('links_dataset.csv')
    G = st.graphA(links_data)
    alpha,beta = st.calc_best_power_law(G)
    print("alpha is {}".format(alpha),"beta is {}".format(beta))
    st.plot_histogram(G, alpha, beta)
    features_dict = st.G_features(G)
    student_result_df = pd.DataFrame(features_dict)
    student_result_df.to_csv(file_name)
    print(features_dict)
    ###########################################################################
    #######################            Part B           #######################
    # ###########################################################################
    ####################### PartB - Q3 ################
    users_data = pd.read_csv('infection_information_set.csv')
    file_name1 = 'result_part_B1.csv'
    file_name2 = 'result_part_B2.csv'
    iterations = 6
    Threshold = 0.41
    WG3 = st.create_undirected_weighted_graph(links_data,users_data,'Q3')
    S = st.run_k_iterations(WG3, iterations, Threshold)
    branching_fac= st.calculate_branching_factor(S, iterations)
    print("Q3 - the dictionary of nodes per iteration is", S )
    print("Q3 - Branching_fac is {}".format(branching_fac))
    branching_fac[0] = 0
    new_row = branching_fac
    for key,value in S.items():
        S[key] = len(value)
    print(S)
    student_result_df = pd.DataFrame(S,index=[0])
    student_result_df = student_result_df.append(new_row, ignore_index=True)
    student_result_df.to_csv(file_name1)
    student_result_df = pd.DataFrame(branching_fac, index=[0])
    student_result_df.to_csv(file_name2)
    # ##################### partB - Q4 ##################

    file_name3 = 'result_part_B3.csv'
    file_name4 = 'result_part_B4.csv'
    file_name5 = 'result_part_B5.csv'
    iterations = 5
    Threshold = 0.3
    h = 600
    WG4 = st.create_undirected_weighted_graph(links_data, users_data,'Q4')
    nodes_dict = st.find_maximal_h_deg(WG4,h)
    nodes_dict_CC = st.calculate_clustering_coefficient(WG4,nodes_dict)
    ##################### find clustered_dsc iterations############
    nodes_dict_sorted_dsc_CC = dict(sorted(nodes_dict_CC.items(), key=lambda item: item[1], reverse=True))
    iter = 1
    infected_number_nodes_sorted_dsc_CC = dict()
    while iter < 31 :
        nodes_list = st.slice_dict(nodes_dict_sorted_dsc_CC,(iter*20))
        WG4 = st.nodes_removal(WG4,nodes_list)
        infected_number_nodes_sorted_dsc_CC[(iter*20)] = st.infected_nodes_after_k_iterations(WG4,iterations,Threshold)
        WG4 = st.create_undirected_weighted_graph(links_data, users_data, 'Q4')
        iter = iter + 1
    student_result_df = pd.DataFrame(infected_number_nodes_sorted_dsc_CC, index=[0])
    student_result_df.to_csv(file_name3)
    ##################### find clustered_asc iterations############
    nodes_dict_sorted_asc_CC = dict(sorted(nodes_dict_sorted_dsc_CC.items(), key=lambda item: item[1], reverse=False))
    iter = 1
    infected_number_nodes_sorted_asc_CC = dict()
    while iter < 31:
        nodes_list = st.slice_dict(nodes_dict_sorted_asc_CC,(iter*20))
        WG4 = st.nodes_removal(WG4,nodes_list)
        infected_number_nodes_sorted_asc_CC[(iter*20)] = st.infected_nodes_after_k_iterations(WG4,iterations,Threshold)
        WG4 = st.create_undirected_weighted_graph(links_data, users_data, 'Q4')
        iter = iter + 1
    print("the dictionary of iteration asc CC:",infected_number_nodes_sorted_asc_CC)
    student_result_df = pd.DataFrame(infected_number_nodes_sorted_asc_CC, index=[0])
    student_result_df.to_csv(file_name4)
    # #################### find random sliced iterations#######
    iter = 1
    dict_random_nodes = dict(sorted(nodes_dict_sorted_dsc_CC.items(), key=lambda item: item[0],  reverse=False))
    infected_number_nodes_random = dict()
    while iter < 31 :
        nodes_list = st.slice_dict(dict_random_nodes,(iter*20))
        WG4 = st.nodes_removal(WG4,nodes_list)
        infected_number_nodes_random[(iter*20)] = st.infected_nodes_after_k_iterations(WG4,iterations,Threshold)
        WG4 = st.create_undirected_weighted_graph(links_data, users_data, 'Q4')
        iter = iter + 1
    print("the dictionary of iteration random:",infected_number_nodes_random)
    student_result_df = pd.DataFrame(infected_number_nodes_random, index=[0])
    student_result_df.to_csv(file_name5)
    # ##################plot a graph ###########################
    ##### uncomment and run it only on your computer and add the plot to the pdf , don't run at the VM ####
    st.graphB(infected_number_nodes_random,infected_number_nodes_sorted_asc_CC,infected_number_nodes_sorted_dsc_CC)
 ################# good luck    ###########################
if __name__ == '__main__':
    main()
