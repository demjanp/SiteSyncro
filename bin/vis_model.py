from sitesyncro import Model

from sitesyncro.utils.fnc_phase import (reduce_earlier_than)

from matplotlib import pyplot
import networkx as nx
from sitesyncro import pygraphviz
import os
import pickle
import copy


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


def get_graph(earlier_than, samples, gap = 0.05):
	
	G = nx.convert_matrix.from_numpy_array(earlier_than, create_using=nx.DiGraph)
	pos = pygraphviz_layout(G, prog='dot')  # {i: [x,y], ...}
	
	x_range = max(x for x, _ in pos.values()) - min(x for x, _ in pos.values())
	gap *= x_range
	
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
		x_last = max(pos[n][0] for n in group) + gap
	
	x_last = 0
	for n in sorted(pos.keys(), key = lambda n: pos[n][0]):
		if pos[n][0] < x_last + gap:
			pos[n] = (x_last + gap, pos[n][1])
			x_last = pos[n][0]
	
	groups = dict([(n, i) for i, group in enumerate(groups) for n in group])
	
	return G, pos, groups, gap


def update_pos_by_dists(pos, samples, dists, group_lines):
	
	pos = copy.deepcopy(pos)
	t_range = max(dists[s][0] for s in samples) - min(dists[s][0] for s in samples)
	height = max(y for x, y in pos.values()) - min(y for x, y in pos.values())
	for i, s in enumerate(samples):
		pos[i] = (pos[i][0] * height / t_range, dists[s][0])
	
	group_lines_ = [x * height / t_range for x in group_lines]
	
	return pos, group_lines_


def get_G(earlier_than):
	
	earlier_than = reduce_earlier_than(earlier_than)
	return nx.convert_matrix.from_numpy_array(earlier_than, create_using=nx.DiGraph)


def plot_graph(G, pos, node_color, edge_color, dists0, dists, samples, group_lines, name):
	pyplot.figure(figsize=(20, 10))
	x_values = [pos[i][0] for i in range(len(samples))]
	if dists:
		graph_width = max(x_values) - min(x_values)
		whisker_length = graph_width * 0.005
		for i, s in enumerate(samples):
			x = pos[i][0]
			
			m, r1, r2 = dists0[s]
			pyplot.errorbar(x, m, yerr=[[m-r2],[r1-m]], fmt='.', color='lightgrey', markersize=5, linewidth=2, zorder=1)
			pyplot.hlines([r1, r2], x - whisker_length, x + whisker_length, color='lightgrey', linewidth=2, zorder=1)
			
			m, r1, r2 = dists[s]
			pyplot.errorbar(x, m, yerr=[[m-r2],[r1-m]], fmt='.', color='k', markersize=5, linewidth=.5, zorder=2)
			pyplot.hlines([r1, r2], x - whisker_length, x + whisker_length, color='k', linewidth=.5, zorder=2)
	
	
	nodes = nx.draw_networkx_nodes(G, pos, node_color=node_color, node_size=25)
	nx.draw_networkx_edges(G, pos, edge_color=edge_color, arrows=True, width=.5)
	nodes.set_zorder(3)
	
	y_min = pyplot.gca().get_ylim()[0]
	y_max = pyplot.gca().get_ylim()[1]
	gap = 0.005 * (pyplot.gca().get_ylim()[1] - y_min)
	last_x = min(x_values)
	for i, x in enumerate(group_lines):
		pyplot.axvline(x, color='grey', linewidth=1, linestyle='--')
		pyplot.text((x + last_x) / 2, y_min + gap, "Gr. %d" % (i+1), verticalalignment='top', horizontalalignment='center', fontsize=8)
		last_x = x
	
	for i, s in enumerate(samples):
		pyplot.text(pos[i][0], y_max + gap, s.split("_")[0], rotation=90, verticalalignment='top', horizontalalignment='center', fontsize=8)
	
	pyplot.gca().invert_yaxis()
	pyplot.tight_layout()
	pyplot.savefig("%s.png" % name)
	pyplot.close()


fdata = "vis_model_data.pickle"

if __name__ == '__main__':
	
	if not os.path.isfile(fdata):
		print("Loading data")
		model0 = Model(directory="stage_0")
		model1 = Model(directory="stage_1")
		model2 = Model(directory="stage_2")
		model3 = Model(directory="stage_3")
		
		samples = list(model0.samples.keys())
		
		eaps = dict([(i, model0.samples[s].excavation_area_phase) for i, s in enumerate(samples)])
		groups = dict([(i, model1.samples[s].group) for i, s in enumerate(samples)])
		phases1 = dict([(i, model1.samples[s].phase) for i, s in enumerate(samples)])
		phases2 = dict([(i, model2.samples[s].phase) for i, s in enumerate(samples)])
		phases3 = dict([(i, model3.samples[s].phase) for i, s in enumerate(samples)])
		
		outliers = [samples.index(s) for s in model1.outliers]
		data = model2.clusters[model2.cluster_opt_n]
		clusters = {}
		for label in data:
			for s in data[label]:
				clusters[samples.index(s)] = label
		
		print("Getting ranges")
		dists0 = {}
		for name in samples:
			m = model0.samples[name].likelihood_mean
			r1, r2 = model0.samples[name].likelihood_range
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
		earlier_than0, samples0 = model0.mphasing.create_earlier_than_matrix()
		assert(samples == samples0)
		
		earlier_than1, samples1 = model1.mphasing.create_earlier_than_matrix()
		assert(samples == samples1)
		
		earlier_than2 = model1.mphasing.update_earlier_than_by_dating(earlier_than1, samples)
		
		earlier_than3, samples3 = model3.mphasing.create_earlier_than_matrix()
		assert(samples == samples3)
#		earlier_than3 = model2.mphasing.update_earlier_than_by_clustering(earlier_than2, samples)
		
		with open(fdata, "wb") as f:
			pickle.dump([
				samples, eaps, outliers, clusters,
				phases1, phases2, phases3,
				dists0, dists1, dists2, dists3, 
				earlier_than0, earlier_than1, earlier_than2, earlier_than3
			], f, pickle.HIGHEST_PROTOCOL, fix_imports = False)
	else:
		with open(fdata, "rb") as f:
			samples, eaps, outliers, clusters, phases1, phases2, phases3, dists0, dists1, dists2, dists3, earlier_than0, earlier_than1, earlier_than2, earlier_than3 = pickle.load(f, fix_imports = False)
	
	G0, pos, groups, gap = get_graph(earlier_than0, samples)
	G1 = get_G(earlier_than1)
	G2 = get_G(earlier_than2)
	G3 = get_G(earlier_than3)
	
	group_lines = []
	last_group = None
	for n in sorted(pos.keys(), key=lambda n: pos[n][0]):
		if groups[n] != last_group:
			group_lines.append(pos[n][0] - gap / 2)
			last_group = groups[n]
	group_lines.append(pos[n][0])
	group_lines = group_lines[1:]
	
	node_colors0 = []
	for i in G0.nodes():
		node_colors0.append('r' if i in outliers else 'k')
	
	edge_colors2 = []
	for i, j in G2.edges():
		if (not earlier_than1[i, j]) and (groups[i] != groups[j]):
			color = 'r'
		else:
			color = 'k'
		edge_colors2.append(color)
	
	edge_colors3 = []
	for i, j in G3.edges():
		if not earlier_than2[i,j]:
			color = 'k'
		else:
			color = 'None'
		edge_colors3.append(color)
	
	pos0, group_lines0 = update_pos_by_dists(pos, samples, dists0, group_lines)
	pos1, group_lines1 = update_pos_by_dists(pos, samples, dists1, group_lines)
	pos2, group_lines2 = update_pos_by_dists(pos, samples, dists2, group_lines)
	pos3, group_lines3 = update_pos_by_dists(pos, samples, dists3, group_lines)
	
	plot_graph(G0, pos, 'k', 'k', dists0,{}, samples, group_lines, 'g0a')
	plot_graph(G0, pos0, 'k', 'k', dists0, dists0, samples, group_lines0, 'g0b')
	plot_graph(G1, pos0, node_colors0, 'k', dists0, dists0, samples, group_lines0, 'g1a')
	plot_graph(G1, pos1, node_colors0, 'k', dists0, dists1, samples, group_lines1, 'g1b')
	plot_graph(G2, pos1, node_colors0, edge_colors2, dists0, dists1, samples, group_lines1, 'g2a')
	plot_graph(G2, pos2, node_colors0, edge_colors2, dists0, dists2, samples, group_lines2, 'g2b')
	plot_graph(G3, pos2, node_colors0, edge_colors3, dists0, dists2, samples, group_lines2, 'g3a')
	plot_graph(G3, pos3, node_colors0, edge_colors3, dists0, dists3, samples, group_lines3, 'g3b')
	
