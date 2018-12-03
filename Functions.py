
def degree_function(graph):
    degree = {}
    for dise, gene in graph.items():
        degree[dise] = len(gene)
    degree_distribution = {}
    for d, deg in degree.items():
        count = 0
        if deg in degree_distribution:
            continue
        for k, v in degree.items():
            if deg == v:
                count += 1
        degree_distribution[deg] = count
    return degree, degree_distribution


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

import pdb
def read_bipartite(file_name):
    genes = {}
    diseas = {}
    lst_diseas = []
    lst_genes = []
    with open(file_name, 'r') as f:
        for line in f:
            (g, d) = line.strip().split('\t')[0:2]
            if d not in diseas:
                diseas[d] = [g]
                lst_diseas.append(d)
            else:
                diseas[d].append(g)
            if g not in genes:
                genes[g] = [d]
                lst_genes.append(g)
            else:
                genes[g].append(d)
    return diseas, genes, lst_diseas, lst_genes

def read_models(data):
    file=[]
    graph={}
    dataset={}
    with open(data, 'r') as f:
        for line in f:
            src, des, weight= line.strip().split('\t')
            if src not in graph:
                graph[src] = [des]
            else:
                graph[src].append(des)
            if des not in graph:
                graph[des]=[src]
            else:
                graph[des].append(src)
            node=graph.keys()
            key=(src, des)
            if key not in dataset:
                dataset[key]=float(weight)
            file.append(float(weight))
    return dataset, graph, file, node

#to normalize values between 0 and 1
def normalization(file):
    normalized_file = []
    for i in range(0, len(file)):
        normalized_file.append((file[i] - min(file)) / (max(file) - min(file)))
    return normalized_file