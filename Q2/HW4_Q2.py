# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <headingcell level=1>

# Homework 4, Question 2

# <codecell>

import networkx as nx
from collections import defaultdict,deque
from random import choice,shuffle
from ggplot import *

# <markdowncell>

# Create a small test graph based on lecture slides.

# <codecell>

tg = nx.Graph()
tg.add_edges_from([(0,1),(0,2),(0,3),(0,4),(1,2),(1,5),(2,5),(3,6),(3,7),(4,7),(6,8),(5,8),(6,9),(7,9),(8,10),(9,10)])

# <markdowncell>

# Perform a BFS from a vertex s

# <codecell>

def BFS(g,v):
    # Book-keeping stuff
    distance = {}
    sigma = defaultdict(float)
    parents = defaultdict(set)
    queue = deque([v])
    done = set()
    
    distance[v] = 0
    sigma[v] = 1.0
    
    induced = nx.DiGraph()
    
    while len(queue) > 0:
        n = queue.popleft()
        done.add(n)
        # Add neighbors of n that are not done to the queue
        for nbr in g.neighbors_iter(n):
            if nbr not in done:
                if nbr in distance and distance[nbr] <= distance[n]: 
                    pass
                else: 
                    distance[nbr] = distance[n]+1
                    if n not in parents[nbr]:
                        parents[nbr].add(n)
                        sigma[nbr] += sigma[n]
                        induced.add_edge(n,nbr)
                    queue.append(nbr)
    return sigma,induced

# <markdowncell>

# Now that we have the induced tree, calculate dependency ($\delta_s$) of an edge $\{v,w\}$ recursively.

# <codecell>

def dependency(v,w,induced,sigma):
    if induced.has_edge(v,w):
        d = sigma[v]/sigma[w]
        if len(induced.neighbors(w)) > 0:
            return d*(1+sum([dependency(w,c,induced,sigma) for c in induced.neighbors(w)]))
        else:
            return d 
    else:
        return 0.0

# <codecell>

def dependencies(graph):
    dependencies = defaultdict(list)
    for src in graph.nodes_iter():
        sigma,induced = BFS(graph,src)
        for v,w in graph.edges_iter():
            dependencies[src].append((v,w,dependency(v,w,induced,sigma)))
    return dependencies

# <headingcell level=2>

# Algorithm 1: Exact Betweenness

# <codecell>

def exact_betweenness(graph):
    betweenness = defaultdict(float) 
    depends = dependencies(graph)
    for src in depends:
        for deps in depends[src]:
            v,w,d = deps
            betweenness[(v,w)] += d
    return dict(betweenness)

# <codecell>

exact_betweenness(tg)

# <codecell>

nx.edge_betweenness_centrality(tg,normalized=False,weight=1)

# <headingcell level=2>

# Algorithm 2: Approximate Betweennes

# <codecell>

def approx_betweenness(graph,c=5,s=10):
    nodes = graph.nodes()
    shuffle(nodes)
    N = float(len(nodes))
    max_D = c*N
    max_samples = int(N/s)
        
    betweenness = defaultdict(float)
    k = defaultdict(int)
    
    n_samples = 0
    for src in nodes[:max_samples]:
        sigma,induced = BFS(graph,src)
        change = False
        for v,w in graph.edges_iter():
            if betweenness[(v,w)] <= max_D:
                betweenness[(v,w)] += dependency(v,w,induced,sigma)
                k[(v,w)] += 1 
                change = True
        if not change:break
    
    return {edge:((N/k[edge])*betweenness[edge]) for edge in graph.edges_iter()}

# <codecell>

approx_betweenness(tg,s=1)

# <headingcell level=2>

# Run on preferential attachment graph

# <codecell>

g = nx.barabasi_albert_graph(1000,4)

# <codecell>

exact_b = exact_betweenness(g)

# <codecell>

approx_b = approx_betweenness(g)

# <codecell>

approx_exact_b = approx_betweenness(g,s=1)

# <codecell>

nx_b = nx.betweenness_centrality(g,normalized=False,weight=1)

# <codecell>

x = range(len(exact_b))

# <codecell>

y1 = sorted(exact_b.values(), reverse=True)

# <codecell>

y2 = sorted(approx_b.values(), reverse=True)

# <codecell>

y3 = sorted(approx_exact_b.values(), reverse=True)

# <codecell>

y4 = sorted(nx_b.values(), reverse=True)

# <codecell>

p=plt.plot(x,y1,label="Exact")
p=plt.plot(x,y2,'--',label="Approximate")
plt.yscale('log')
plt.legend()
plt.title("Comparison of exact and approximate betweenness centrality")
plt.xlabel("Centrality ranking")
plt.ylabel("Centrality value")

# <codecell>


