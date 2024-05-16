from sitesyncro import Model

from sitesyncro.utils.fnc_phase import (reduce_earlier_than)

from matplotlib import pyplot
import networkx as nx
import numpy as np
from sitesyncro import pygraphviz
import os

def to_agraph(N):
	
	directed = N.is_directed()
	strict = nx.number_of_selfloops(N) == 0 and not N.is_multigraph()
	A = pygraphviz.AGraph(name=N.name, strict=strict, directed=directed)

	# default graph attributes
	A.graph_attr.update(N.graph.get("graph", {}))
	A.node_attr.update(N.graph.get("node", {}))
	A.edge_attr.update(N.graph.get("edge", {}))

	A.graph_attr.update(
		(k, v) for k, v in N.graph.items() if k not in ("graph", "node", "edge")
	)

	# add nodes
	for n, nodedata in N.nodes(data=True):
		A.add_node(n)
		# Add node data
		a = A.get_node(n)
		a.attr.update({k: str(v) for k, v in nodedata.items()})

	# loop over edges
	if N.is_multigraph():
		for u, v, key, edgedata in N.edges(data=True, keys=True):
			str_edgedata = {k: str(v) for k, v in edgedata.items() if k != "key"}
			A.add_edge(u, v, key=str(key))
			# Add edge data
			a = A.get_edge(u, v)
			a.attr.update(str_edgedata)

	else:
		for u, v, edgedata in N.edges(data=True):
			str_edgedata = {k: str(v) for k, v in edgedata.items()}
			A.add_edge(u, v)
			# Add edge data
			a = A.get_edge(u, v)
			a.attr.update(str_edgedata)

	return A

def pygraphviz_layout(G, prog="dot", root=None, args=""):
	
	if root is not None:
		args += f"-Groot={root}"
	A = to_agraph(G)
	A.layout(prog=prog, args=args)
	node_pos = {}
	for n in G:
		node = pygraphviz.Node(A, n)
		try:
			xs = node.attr["pos"].split(",")
			node_pos[n] = tuple(float(x) for x in xs)
		except:
			print("no position for node", n)
			node_pos[n] = (0.0, 0.0)
	return node_pos

def get_graph(earlier_than, samples, dists, gap = 0.2):
	
	t_range = max(dists[s][0] for s in samples) - min(dists[s][0] for s in samples)
	G = nx.convert_matrix.from_numpy_array(earlier_than, create_using=nx.DiGraph)
	pos = pygraphviz_layout(G, prog='dot')  # {i: [x,y], ...}
	x_range = max(x for x, _ in pos.values()) - min(x for x, _ in pos.values())
	
	done = set()
	groups = []
	for c in nx.connected_components(G.to_undirected()):
		G_sub = G.subgraph(c)
		groups.append(list(G_sub.nodes()))
		done.update(groups[-1])
	for n in pos:
		if n not in done:
			groups.append([n])
	x_last = 0
	for group in groups:
		x_min = min(pos[n][0] for n in group)
		for n in group:
			pos[n] = ((pos[n][0] - x_min) + x_last, pos[n][1])
		x_last = max(pos[n][0] for n in group) + x_range * gap
	
	height = max(y for x, y in pos.values()) - min(y for x, y in pos.values())
	for i, s in enumerate(samples):
		pos[i] = (pos[i][0] * height / t_range, dists[s][0])
	
	labels = {i: s for i, s in enumerate(samples)}
	return G, pos, labels

if __name__ == '__main__':
	
	
	print("Loading data")
	model1 = Model(directory="stage_1")
	model2 = Model(directory="stage_2")
	model3 = Model(directory="stage_3")
	
	samples = list(model1.samples.keys())
	
	print("Getting ranges")
	dists0 = {}
	for name in samples:
		m = model1.samples[name].likelihood_mean
		r1, r2 = model1.samples[name].likelihood_range
		dists0[name] = [m, r1, r2]
	dists1 = {}
	for name in samples:
		m = model1.samples[name].posterior_mean
		r1, r2 = model1.samples[name].posterior_range
		dists1[name] = [m, r1, r2]
	dists2 = {}
	for name in samples:
		m = model2.samples[name].posterior_mean
		r1, r2 = model2.samples[name].posterior_range
		dists2[name] = [m, r1, r2]
	dists3 = {}
	for name in samples:
		m = model3.samples[name].posterior_mean
		r1, r2 = model3.samples[name].posterior_range
		dists3[name] = [m, r1, r2]
	
	print("Populating graph")
	earlier_than1, samples1 = model1.mphasing.create_earlier_than_matrix()
	assert(samples == samples1)
	earlier_than1 = reduce_earlier_than(earlier_than1)
	
	earlier_than2, samples2 = model2.mphasing.create_earlier_than_matrix()
	assert(samples == samples2)
	earlier_than2 = reduce_earlier_than(earlier_than2)

	earlier_than3, samples3 = model3.mphasing.create_earlier_than_matrix()
	assert(samples == samples3)
	earlier_than3 = reduce_earlier_than(earlier_than3)
	
	G1, pos1, labels1 = get_graph(earlier_than1, samples, dists1)
	
	nx.draw(G1, pos1, labels=labels1, with_labels=True, arrows=True, font_size=6)
	pyplot.show()
	