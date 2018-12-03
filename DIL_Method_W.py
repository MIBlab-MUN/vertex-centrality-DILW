import sys
import pdb, math
import matplotlib.pyplot as plt

from creat_plot import creat_plot


from Functions import  degree_of_node_weighted
from Functions import read_models, normalization, degree_function

first_arg = sys.argv[1] #graph to compute dil
#  Uncomment to read node importance scores calculated by other measures in format of    NodeName	Score
#  in order to be compared with DIL-W method
#  second_arg_cc = sys.argv[2]

from scipy.stats import spearmanr

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

def DIL_Weighted_Opsahl(graph, ds, alpha):
    degree= Opsahl_method(alpha, ds, graph)
    # degree=degree_of_node_weighted(graph, ds)
    # sum_link = degree_of_node_weighted(graph, ds)
    # degree, dd = degree_function(graph)

    importance_node = {}
    importance_of_line = {}
    for key, item in graph.items():
        sigma = 0
        for node in item:
            if (node, key) not in importance_of_line:
                common=list(set(graph[key]).intersection(graph[node]))
                sum_key=0
                sum_node=0
                for vertex in common:
                    if (vertex, node) in ds:
                        sum_node+=ds[(vertex, node)]
                    else:
                        sum_node+=ds[(node, vertex)]
                    if (vertex, key) in ds:
                        sum_key+=ds[(vertex, key)]
                    else:
                        sum_key += ds[(key, vertex)]
                triangle=sum_key+sum_node
                try:
                    p_key=float((math.pow(len(common), 1 - float(alpha)) * math.pow(sum_key, float(alpha))))
                    p_node=float((math.pow(len(common), 1 - float(alpha)) * math.pow(sum_node, float(alpha))))
                except:
                    p_key = 0
                    p_node = 0
                u = (degree[key] - p_key) * (degree[node] - p_node)
                # landa = triangle / 2 + 1
                landa=float(p_key+p_node)/2+1
                Iij = float(u) / float(landa)
                # if (key, node)==('a', 'c'):
                #     pdb.set_trace()
                importance_of_line[(key, node)] = Iij
                sigma = sigma + (Iij * (float(degree[key]) / float(degree[key] + degree[node])))
            else:
                sigma = sigma + (importance_of_line[node, key] * (float(degree[key]) / float(degree[key] + degree[node])))

        importance_node[key] = degree[key] + sigma
        # importance_node[key] = len(graph[key]) + sigma
    # pdb.l()
    # print(importance_of_line)
    # pdb.set_trace()
    return importance_node


def DIL_Weighted(graph, ds):
    degree = degree_of_node_weighted(graph, ds)
    # degree, dd = degree_function(graph)

    importance_node = {}
    importance_of_line = {}
    for key, item in graph.items():
        sigma = 0
        for node in item:
            if (node, key) not in importance_of_line:
                common=list(set(graph[key]).intersection(graph[node]))
                sum_key=0
                sum_node=0
                for vertex in common:
                    if (vertex, node) in ds:
                        sum_node+=ds[(vertex, node)]
                    else:
                        sum_node+=ds[(node, vertex)]
                    if (vertex, key) in ds:
                        sum_key+=ds[(vertex, key)]
                    else:
                        sum_key += ds[(key, vertex)]
                triangle=sum_key+sum_node
                u = (degree[key] - sum_key) * (degree[node] - sum_node)

                # triangle = len(list(set(graph[key]).intersection(graph[node])))
                # u = (degree[key] - triangle - 1) * (degree[node] - triangle - 1)

                landa = triangle / 2 + 1
                Iij = float(u) / float(landa)
                # if (key, node)==('a', 'c'):
                #     pdb.set_trace()
                importance_of_line[(key, node)] = Iij
                # sigma = sigma + (Iij * (float(degree[key] - 1) / float(degree[key] + degree[node] - 2)))
                sigma = sigma + (Iij * (float(degree[key]) / float(degree[key] + degree[node])))
                # sigma = sigma + (Iij * (float(len(graph[key])) / float(len(graph[key]) + len(graph[node]))))

            else:
                # sigma = sigma + (importance_of_line[node, key] * (float(degree[key] - 1) / float(degree[key] + degree[node] - 2)))
                sigma = sigma + (importance_of_line[node, key] * (float(degree[key]) / float(degree[key] + degree[node])))
                # sigma = sigma + (importance_of_line[node, key] * (float(len(graph[key])) / float(len(graph[key]) + len(graph[node]))))

        importance_node[key] = degree[key] + sigma
        # importance_node[key] = len(graph[key]) + sigma
    # pdb.set_trace()
    return importance_node


if __name__ == '__main__':
    ds, graph, c, d= read_models(data=first_arg)
    # lst_obsahl_alpha=[0.5, 1, 1.5]
    lst_obsahl_alpha = [1]
    for obsahl_alpha in lst_obsahl_alpha:
        importance_node = DIL_Weighted_Opsahl(graph, ds, obsahl_alpha)
        importance_node_sorted = sorted(importance_node.items(), key=lambda kv: kv[1], reverse=True)
        
		try:
            type = (first_arg.split('.')[0]).split('_')[1]
			# If you want to compare DIL-W method with other meaures uncommpent lines below
            # alpha=(second_arg_cc.split('.')[1]).split('_')[2]
        except:
            type=""
		# If you want to compare DIL-W method with other meaures uncommpent lines below
        #     alpha=""
		

        with open("DIL_"+str(obsahl_alpha)+"_"+type+".txt", "w") as myfile:
            for i in range(0, len(importance_node_sorted)):
                key, value = importance_node_sorted[i]
                myfile.write(str(key) + "\t" + str(value) + "\n")

#For calculating DIL_W as normal formula uncomment following commands
    # importance_node= DIL_Weighted(graph, ds)
    # importance_node_sorted= sorted(importance_node.items(), key=lambda kv: kv[1], reverse=True)
    # try:
    #     type = (first_arg.split('.')[0]).split('_')[1]
    #     alpha = (first_arg.split('.')[0]).split('_')[3]
    # except:
    #     type=""
    #     alpha=""
    # with open("DIL_"+type+"_"+alpha+".txt", "w") as myfile:
    #     for i in range(0, len(importance_node_sorted)):
    #         key, value = importance_node_sorted[i]
    #         myfile.write(str(key) + "\t" + str(value) + "\n")
#---------------end ---------------#


# Comparison of DIL_W vs CC
#     ds_dc = second_arg_cc
#     CC = {}
#     with open(ds_dc, 'r') as f:
#         for line in f:
#             node=line.strip().split('\t')[0]
#             closeness = line.strip().split('\t')[1]
#             # node, closeness, NCC = line.strip().split('\t')
#             CC[node]=closeness
#     x_DIL = []
#     y_closeness = []
#     for key, item in CC.items():
#         x_DIL.append(importance_node[key])
#         y_closeness.append(float(item))
#
#     # x_DIL=normalization(x_DIL)
#     y_closeness=normalization(y_closeness)
#
#     x_DIL2=[]
#     y_closeness2=[]
#     for i in range (0, len(x_DIL)):
#         x_DIL2.append(x_DIL[i]/sum(x_DIL))
#         y_closeness2.append(y_closeness[i]/sum(y_closeness))
#
#     # pdb.set_trace()
#     print(spearmanr(x_DIL, y_closeness))
#
#     plt.plot(x_DIL, y_closeness, 'o', mfc='none', color='red')
#     plt.xlabel('EDIL-W centrality')
#     plt.ylabel('Closeness centrality')
#
#     plt.xscale('log')
#     # plt.yscale('log')
#
#     # plt.xlim(0.5, 1000)
#     # plt.xlim(-0.01, 1.01)
#     # plt.ylim(-0.01, 1.01)
#     # plt=creat_plot(x_DIL, y_closeness,'o', 'DIL_W Centrality', 'Closeness Centrality', 'Correlation between DIL and CC in '+type+' '+alpha)
#
#     plt.savefig("CC_vs_DIL_" + type + "_" + alpha + ".eps")
#     plt.show()
#     plt.clf()

#Comparison of DIL_W vs DC
    # # dc_w = degree_of_node_weighted(graph, ds)
    # # pdb.set_trace()
    # dc_w = Opsahl_method(obsahl_alpha, ds, graph)
    # x_DIL = []
    # y_DC = []
    # for key, item in dc_w.items():
    #     x_DIL.append(importance_node[key])
    #     y_DC.append(float(item))
    #
    # # x_DIL=normalization(x_DIL)
    # # y_DC=normalization(y_DC)
    #
    # x_DIL2=[]
    # y_DC2=[]
    # for i in range (0, len(x_DIL)):
    #     x_DIL2.append(x_DIL[i]/sum(x_DIL))
    #     y_DC2.append(y_DC[i]/sum(y_DC))
    #
    # plt.plot(x_DIL2, y_DC2, 'o', mfc='none', color='red')
    # plt.xlabel('EDIL-W centrality')
    # plt.ylabel('Degree centrality')
    # plt.xscale('log')
    # plt.yscale('log')
    # # plt.xlim(0.5, 1000)
    # # plt.xlim(-0.01, 1.01)
    # # plt.ylim(-0.01, 1.01)
    # # plt=creat_plot(x_DIL, y_closeness,'o', 'DIL_W Centrality', 'Closeness Centrality', 'Correlation between DIL and CC in '+type+' '+alpha)
    #
    # plt.savefig("DC_vs_DIL_" + type + "_" + alpha + ".eps")
    # plt.show()
    # plt.clf()


#Comparison of DIL_W vs BC
    # ds_dc = second_arg_cc
    # BC = {}
    # with open(ds_dc, 'r') as f:
    #     for line in f:
    #         node, betweenness = line.strip().split('\t')
    #         BC[node]=betweenness
	# 
    # x_DIL = []
    # y_betweenness = []
    # count_zero=0
    # for key, item in BC.items():
    #     x_DIL.append(importance_node[key])
    #     y_betweenness.append(float(item))
    #     if item =="0":
    #         count_zero+=1
	# 
    # print (count_zero)
    # # x_DIL=normalization(x_DIL)
    # # y_betweenness=normalization(y_betweenness)
	# 
    # x_DIL2=[]
    # y_betweenness2=[]
    # for i in range (0, len(x_DIL)):
    #     x_DIL2.append(x_DIL[i]/sum(x_DIL))
    #     y_betweenness2.append(y_betweenness[i]/sum(y_betweenness))
	# 
    # print(spearmanr(x_DIL, y_betweenness))
	# 
    # # pdb.set_trace()
    # plt.plot(x_DIL, y_betweenness, 'o', mfc='none', color='red')
    # plt.xlabel('EDIL-W centrality')
    # plt.ylabel('Betweenness centrality')
    # plt.xscale('log')
    # plt.yscale('log')
	# 
    # # plt.xlim(0.5, 1000)
    # # plt.ylim(-100, 10000)
	# 
    # plt.savefig("BC_vs_DIL_" + type + "_" + alpha + ".eps")
    # plt.show()
    # plt.clf()