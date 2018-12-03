import sys
import os
import pdb
import math

first_arg = sys.argv[1]


def component(graph, node):
    components = {}
    count = 0
    while node:
        start = node[0]
        count += 1
        q = [start]
        path = [start]
        while q:
            v = q.pop(0)
            node.remove(v)
            for u in graph[v]:
                if u not in path:
                    path.append(u)
                    q.append(u)
                graph[u].remove(v)
        components[count] = path
    return count, components


def cocitation_opsahl(graph, lst_nodes, opsahl_result, ds):
    common_nodes = {}
    constructed_graph = {}
    for i in range(0, len(lst_nodes)):
        # pdb.set_trace()
        print("{0}/{1}".format(i, len(lst_nodes)))
        for j in range(i + 1, len(lst_nodes)):
            common_nodes[(lst_nodes[i], lst_nodes[j])] = list(set(graph[lst_nodes[i]]).intersection(graph[lst_nodes[j]]))
            # if len(common_nodes[(lst_nodes[i], lst_nodes[j])])==2 and lst_nodes[i] != 'C0024117' and lst_nodes[j] != 'C0024117' and lst_nodes[i] !='C0017567' and lst_nodes[j] !='C0017567':
            #     pdb.set_trace()
            sumi = 0
            for node in common_nodes[(lst_nodes[i], lst_nodes[j])]:
                # pdb.set_trace()
                if (lst_nodes[i], node) in ds:
                    w1=float(ds[(lst_nodes[i], node)])
                else:
                    w1=float(ds[(node, lst_nodes[i])])
                if (lst_nodes[j], node) in ds:
                    w2=float(ds[(lst_nodes[j], node)])
                else:
                    w2=float(ds[(node, lst_nodes[j])])
                # pdb.set_trace()
                sumi += (w1+w2)/opsahl_result[node]
            constructed_graph[(lst_nodes[i], lst_nodes[j])] = sumi
    # pdb.set_trace()
    return constructed_graph


def read_models(data):
    file = []
    graph = {}
    dataset = {}
    with open(data, 'r') as f:
        for line in f:
            src, des, weight = line.strip().split('\t')
            if src not in graph:
                graph[src] = [des]
            else:
                graph[src].append(des)
            if des not in graph:
                graph[des] = [src]
            else:
                graph[des].append(src)
            node = graph.keys()
            key = (src, des)
            if key not in dataset:
                dataset[key] = weight
            file.append(float(weight))
    return dataset, graph, file, node


def degree_of_node(file):
    degree = {}
    for node, neighbour in file.items():
        degree[node] = len(neighbour)
    return degree


def Opsahl_method(alpha, weights, graph):
    weight_degree = degree_of_node_weighted(graph, weights)
    degree = degree_of_node(graph)
    opsahl_result = {}
    for node in graph:
        opsahl_result[node] = float((math.pow(degree[node], 1 - float(alpha)) * math.pow(weight_degree[node], float(alpha))))
    return opsahl_result


def degree_of_node_weighted(file, weights):
    weight_degree = {}
    for node, neighbours in file.items():
        sumi = 0
        for neighbour in neighbours:
            if (node, neighbour) in weights:
                sumi += float(weights[(node, neighbour)])
            else:
                sumi += float(weights[(neighbour, node)])
        weight_degree[node] = sumi
    return weight_degree


def ds_cleaner(dataset=first_arg):
    # data="sample.txt"
    diseases = {}
    genes = {}
    lst_diseas = []
    lst_genes = []
    network = {}
    with open(dataset, 'r') as f:
        for line in f:
            # g=line.strip().split('\t')[0]
            # d= line.strip().split('\t')[2]
            # w= line.strip().split('\t')[4]
            g = line.strip().split('\t')[0]
            d = line.strip().split('\t')[1]
            w = line.strip().split('\t')[2]
            if (g, d) not in network and (d, g) not in network:
                network[(g, d)] = w
            if d not in diseases:
                diseases[d] = [g]
                lst_diseas.append(d)
            else:
                diseases[d].append(g)
            if g not in genes:
                genes[g] = [d]
                lst_genes.append(g)
            else:
                genes[g].append(d)


            # ----------finding the giant component in the bipartite network-------------#
            # graph=dict(genes.items()+diseases.items())
            # # pdb.set_trace()
            # count, components=component(graph, graph.keys())
            # max_component=0
            # max_index=-1
            # for key, value in components.items():
            # 	components[key]=sorted(value)
            # 	if len(value)>max_component:
            # 		max_component=len(value)
            # 		max_index=key

            # pdb.set_trace()

            # with open("wbipartite.txt",  "w") as myfile:
            # 	for key, item in network.items():
            # 		g,d=key
            # 		if g in components[max_index] and d in components[max_index]:
            # 			myfile.write(g+"\t"+d+"\t"+item+"\n")

            # -----------end of finding giant component ----------------#

    # ----------------wrting a componnet in txt file-----------------#
    # nodes = ['C4014321', 'C1721005', '3860', 'C4011926', '3851']
    # partial_graph = {}
    # for node in nodes:
    #     if node in diseases:
    #         for neigh in diseases[node]:
    #             if (neigh, node) not in partial_graph:
    #                 if (node, neigh) in network:
    #                     partial_graph[(node, neigh)] = network[(node, neigh)]
    #                 else:
    #                     partial_graph[(node, neigh)] = network[(neigh, node)]
    #     else:
    #         for neigh in genes[node]:
    #             if (neigh, node) not in partial_graph:
    #                 if (node, neigh) in network:
    #                     partial_graph[(node, neigh)] = network[(node, neigh)]
    #                 else:
    #                     partial_graph[(node, neigh)] = network[(neigh, node)]
    # with open("partial_graph.txt",  "w") as myfile:
    # 	for key, item in partial_graph.items():
    # 		g,d=key
    # 		myfile.write(g+"\t"+d+"\t"+item+"\n")
    # ----------------end of wrting a componnet in txt file-----------------#


    # pdb.set_trace()
    alphas = [0, 0.5, 1, 1.5]
    alphas = [1]

    for alpha in alphas:
        # #test opsahl
        ds,graph,c,d = read_models(dataset)
        # opsahl=Opsahl_method(alpha, ds, graph)
        # print("alpha: {0}\n".format(alpha))
        # print(opsahl)
        # print("----------")
        # #end test opsahl

        opsahl_genes = Opsahl_method(alpha, network, genes)  # genes show nodes which are connected to each gene
        opsahl_diseases = Opsahl_method(alpha, network, diseases)  # network show weight of edges, diseases show nodes which are connected to each disease
        gene_network = cocitation_opsahl(genes, lst_genes, opsahl_diseases, ds)
        # disease_network = cocitation_opsahl(diseases, lst_diseas, opsahl_genes, ds)
        # pdb.set_trace()
        with open("gene_network_" + str(alpha) + ".txt", "w") as myfile:
            for key, value in gene_network.items():
                if value > 0:
                    src, des = key
                    myfile.write(src + "\t" + des + "\t" + str(value) + "\n")

        # with open("disease_network_" + str(alpha) + ".txt", "w") as myfile:
        #     for key, value in disease_network.items():
        #         if value > 0:
        #             src, des = key
        #             myfile.write(src + "\t" + des + "\t" + str(value) + "\n")


if __name__ == '__main__':
    ds_cleaner()
