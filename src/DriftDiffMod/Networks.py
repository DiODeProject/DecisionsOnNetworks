'''
Created on 19 Dec 2016

@author: Andreagiovanni Reina.
University of Sheffield, UK.
'''

import collections
import numpy as np
import matplotlib.pyplot as plt
import networkx as nx #@UnresolvedImport

# wrap a few graph generation functions so they have the same signature
# n : int  -- The number of nodes.
# m : int --  the number of random edges to add for each new node
# k : int --  Each node is joined with its ``k`` nearest neighbors in a ring topology.
# p : float -- The probability of adding a new edge for each edge.

def random_lobster(n, m, k, p):
    return nx.random_lobster(n, p, p / m)

def powerlaw_cluster(n, m, k, p, randomSeed):
    return nx.powerlaw_cluster_graph(n, m, p, randomSeed)

def erdos_renyi(n, m, k, p):
    return nx.erdos_renyi_graph(n, p)

def newman_watts_strogatz(n, m, k, p):
    return nx.newman_watts_strogatz_graph(n, k, p)

def plot_random_graph(n, m, k, p, generator):
    g = generator(n, m, k, p)
    nx.draw(g)
    plt.show()
    
    
################ TESTING A PLAYING ########################

def plotDensityHistogram(G):
    G = nx.gnp_random_graph(100, 0.02)
    #degree_sequence=sorted([d for n,d in G.degree()], reverse=True) # degree sequence
    tmp=G.degree()
    #degree_sequence=sorted(tmp.items(), key=operator.itemgetter(1))
    #degree_sequence=degree_sequence.reverse()
    degree_sequence=sorted(tmp, key=tmp.get, reverse=True)
    print ("Degree sequence", degree_sequence)
    degreeCount=collections.Counter(degree_sequence)
    deg, cnt = zip(*degreeCount.items())
    
    #fig, ax = plt.subplots()

    plt.bar(deg, cnt, width=0.80, color='b')
    
    plt.title("Degree Histogram")
    plt.ylabel("Count")
    plt.xlabel("Degree")

#     ax.set_xticks([d+0.4 for d in deg])
#     ax.set_xticklabels(deg)
#     
#     # draw graph in inset
#     plt.axes([0.4, 0.4, 0.5, 0.5])
#     Gcc=sorted(nx.connected_component_subgraphs(G), key = len, reverse=True)[0]
#     pos=nx.spring_layout(G)
#     plt.axis('off')
#     nx.draw_networkx_nodes(G, pos, node_size=20)
#     nx.draw_networkx_edges(G, pos, alpha=0.4)
    
    #plt.savefig("degree_histogram.png")
    plt.show()
    
def testNet(numOfNodes, numEdges, randomSeed):
    graph=powerlaw_cluster(numOfNodes, numEdges, 0, 0, randomSeed)
    print( "number of edges: " + str(graph.number_of_edges())  + " avg on " + str(graph.number_of_nodes()) + ": " + str(graph.number_of_edges()/(numOfNodes-2) ))
    print( "density: " + str(nx.density(graph)))
    
    conn = list()
    for node in nx.nodes(graph):
        e = len(list( nx.all_neighbors(graph,node) ))
        conn.append( e )
    
    print(conn)
    connCount = collections.Counter(conn)
    print(connCount)
    print( np.arange( len(connCount) ) )
    #plt.bar( np.arange( len(connCount) ), connCount.values())
    plt.bar( connCount.keys(), connCount.values())
    plt.show()
    #plotDensityHistogram(graph)

if __name__ == '__main__':
    numOfNodes = 30
    numEdges = 4
    randomSeed = 333
    testNet(numOfNodes, numEdges, randomSeed)
    
