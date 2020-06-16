import numpy as np
import networkx as nx
import scipy.sparse as sp


def my_data_load(graph_path, feature_path):
    miRNA_nodes = []
    mRNA_nodes = []

    with open(feature_path, "r") as ff:
        for i in range(88214):
            line = ff.readline()
            if i == 0:
                continue
            line_split = line.strip().split(",")
            if line_split[2] not in miRNA_nodes and line_split[2] != "":
                miRNA_nodes.append(line_split[2])
            if line_split[17] not in mRNA_nodes and line_split[17] != "":
                mRNA_nodes.append(line_split[17])


    cur = 0
    mRNA_graph = nx.Graph()
    weight_edges = []
    with open(graph_path, "r") as fg:
        for line in fg.readlines():
            line_split = line.split()
            if cur == 0:
                cur += 1
                continue
            if len(line_split) > 1:
                # GLYMA_05G019000    Glyma.04G131800
                nodeA = line_split[0].replace("lyma.", "LYMA_")
                nodeB = line_split[1].replace("lyma.", "LYMA_")
                lls = float(line_split[2])
                if nodeA in mRNA_nodes and nodeB in mRNA_nodes:
                    weight_edges.append((nodeA, nodeB, lls))

    mRNA_graph.add_weighted_edges_from(weight_edges)
    # 移除度<16的节点
    remove_nodes = []
    for node in mRNA_graph.nodes:
        if mRNA_graph.degree(node) < 16 :
            remove_nodes.append(node)
    mRNA_graph.remove_nodes_from(remove_nodes)

    nx.write_adjlist(mRNA_graph, "/data/mRNA graph adj.dat")

    adj = nx.adjacency_matrix(mRNA_graph)

    print("adj OK")


    feature_matrix = np.zeros((len(list(mRNA_graph.nodes)), len(miRNA_nodes)), dtype=np.int)
    with open(feature_path, "r") as ff:
        for i in range(88214):  #88214
            line = ff.readline()
            if i == 0:
                continue
            line_split = line.strip().split(",")
            if line_split[2] != "" and line_split[17] != "":
                if line_split[17] in list(mRNA_graph.nodes):
                    row_index = list(mRNA_graph.nodes).index(line_split[17])
                    col_index = miRNA_nodes.index(line_split[2])
                    feature_matrix[row_index][col_index] = 1

    np.savetxt("/data/miRNA-4X4 matrix.txt", np.transpose(feature_matrix), fmt="%1d")

    features = sp.csr_matrix(feature_matrix)

    print("feature OK")

    # 输出miRNA_nodes与mRNA_nodes的name
    help_print(miRNA_nodes, "miRNA nodes")
    help_print(list(mRNA_graph.nodes), "4X4 nodes")

    print(features.shape)
    print(mRNA_graph.number_of_nodes())

    return adj, features


def help_print(one_list, name):
    file = "/data/" + name + ".name"

    with open(file, "w") as fw:
        for elem in one_list:
            fw.write(elem + "\n")


if __name__ == '__main__':
    adj, features = my_data_load("/data/SoyNet.txt",
                 "/data/final gma data dup filled.csv")
