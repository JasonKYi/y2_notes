import numpy as np
import matplotlib.pyplot as plt
import networkx as nx

# Creating the graph using a list of tuples
G = nx.Graph()
e = [(1, 2), (1, 5), (2, 3), (2, 5), (3, 4), (4, 5), (4, 6)]
G.add_edges_from(e)

# Drawing the graph
# plt.figure()
# nx.draw(G, with_labels = True)
# plt.show()

# Adjacency matric
A = nx.adjacency_matrix(G)
# We need to convert it to an proper array first as its currently a sparse array 
# Clustering coefficients as a dictionary
C = nx.clustering(G)
# Degree and degree distributions
D = nx.degree(G)
D_distrib = nx.degree_histogram(G)

# Generate random graph

H = nx.gnp_random_graph(1000, 0.05)
hist = nx.degree_histogram(H)
plt.figure()
plt.plot(hist, "bo")
plt.show()