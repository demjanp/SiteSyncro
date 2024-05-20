from sitesyncro import Model

from sitesyncro.utils.fnc_phase import (reduce_earlier_than)

from itertools import permutations
from collections import defaultdict
from matplotlib import pyplot
import networkx as nx
from sitesyncro import pygraphviz
import numpy as np
import math
import os
import pickle
import copy
from natsort import natsorted


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
	
	def _update_pos_by_groups(pos, groups, gap = 0.05):
		
		x_range = max(x for x, _ in pos.values()) - min(x for x, _ in pos.values())
		gap *= x_range
		
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
		
		return pos
		
	earlier_than = reduce_earlier_than(earlier_than)
	G = nx.convert_matrix.from_numpy_array(earlier_than, create_using=nx.DiGraph)
	pos = pygraphviz_layout(G, prog='dot')  # {i: [x,y], ...}
	
	done = set()
	groups = []
	for c in nx.connected_components(G.to_undirected()):
		G_sub = G.subgraph(c)
		groups.append(list(G_sub.nodes()))
		done.update(groups[-1])
	for n in pos:
		if n not in done:
			groups.append([n])
	
	groups = natsorted(groups, key = lambda group: samples[group[0]])
	pos = _update_pos_by_groups(pos, groups)
	
	groups = dict([(n, i) for i, group in enumerate(groups) for n in group])
	
	return G, pos, groups, gap


def get_G(earlier_than):
	
	earlier_than = reduce_earlier_than(earlier_than)
	return nx.convert_matrix.from_numpy_array(earlier_than, create_using=nx.DiGraph)


def plot_graph(title, G, pos, node_color, edge_color, likelihoods, posteriors, samples, outliers, time_groups, time_groups_title, name, t_step = 100):
	
	def _update_pos_by_dists(pos, samples, dists):
	
		pos = copy.deepcopy(pos)
		t_range = max(dists[s][0] for s in samples) - min(dists[s][0] for s in samples)
		height = max(y for x, y in pos.values()) - min(y for x, y in pos.values())
		for i, s in enumerate(samples):
			pos[i] = (pos[i][0] * height / t_range, dists[s][0])
			
		return pos
	
	if posteriors:
		pos = _update_pos_by_dists(pos, samples, posteriors)
	
	pyplot.figure(figsize=(20, 10))
	pyplot.title(title)
	x_values = [pos[i][0] for i in range(len(samples))]
	x_min = min(x_values)
	x_max = max(x_values)
	step = (x_max - x_min) / len(pos)
	if posteriors:
		graph_width = max(x_values) - min(x_values)
		whisker_length = graph_width * 0.005
		t_min, t_max = np.inf, -np.inf
		for i, s in enumerate(samples):
			x = pos[i][0]
			
			color = "k"
			if i in outliers:
				color = "r"
			m, r1, r2 = likelihoods[s]
			pyplot.errorbar(x, m, yerr=[[m-r2],[r1-m]], fmt='.', color=color, markersize=5, linewidth=2, zorder=1, alpha=0.25)
			pyplot.hlines([r1, r2], x - whisker_length, x + whisker_length, color=color, linewidth=2, zorder=1, alpha=0.25)
			t_min = min(t_min, r2)
			t_max = max(t_max, r1)
			
			m, r1, r2 = posteriors[s]
			pyplot.errorbar(x, m, yerr=[[m-r2],[r1-m]], fmt='.', color='k', markersize=5, linewidth=.5, zorder=2)
			pyplot.hlines([r1, r2], x - whisker_length, x + whisker_length, color='k', linewidth=.5, zorder=2)
			t_min = min(t_min, r2)
			t_max = max(t_max, r1)
		
		t0 = np.floor((1950 - t_max) / t_step) * t_step
		t1 = np.ceil((1950 - t_min) / t_step) * t_step
		yticks = np.arange(t0, t1 + t_step, t_step)
		yticks = list(zip(1950 - yticks, yticks.astype(int)))
		pyplot.ylabel("Year CE", labelpad=25)
		
	else:
		yticks = sorted(list(set([y for _, y in pos.values()])))
		yticks = list(zip(yticks, np.arange(len(yticks), dtype=int)[::-1] + 1))
		pyplot.ylabel("Stratigraphic Phase", rotation = 90, labelpad=10)
	
	nodes = nx.draw_networkx_nodes(G, pos, node_color=node_color, node_size=25)
	nx.draw_networkx_edges(G, pos, edge_color=edge_color, arrows=True, width=.5)
	nodes.set_zorder(3)
	
	y_min = pyplot.gca().get_ylim()[0]
	y_max = pyplot.gca().get_ylim()[1]
	gap = 0.005 * (pyplot.gca().get_ylim()[1] - y_min)
	if not time_groups:
		
		group_lines = []
		last_group = None
		for n in sorted(pos.keys(), key=lambda n: pos[n][0]):
			if groups[n] != last_group:
				group_lines.append(pos[n][0] - step / 2)
				last_group = groups[n]
		group_lines.append(pos[n][0])
		group_lines = group_lines[1:]
		
		last_x = x_min
		for i, x in enumerate(group_lines):
			prefix = "Group" if ((x - last_x) / step) > 1.5 else "Gr."
			pyplot.text((x + last_x) / 2, y_min + gap, "%s %d" % (prefix, i+1), verticalalignment='top', horizontalalignment='center', fontsize=10)
			if i < len(group_lines) - 1:
				pyplot.axvline(x, color='grey', linewidth=1, linestyle='--')
			last_x = x
	
	for i, s in enumerate(samples):
		pyplot.text(pos[i][0], y_max + gap, s.split("_")[0], rotation=90, verticalalignment='top', horizontalalignment='center', fontsize=8)
	pyplot.xlabel("Sample", labelpad=60)
	
	for y, label in yticks:
		pyplot.text(step / 2, y, str(label), verticalalignment='top', horizontalalignment='right', fontsize=8)
	
	if time_groups:
		time_group_lines = []
		last_group = None
		last_y = 0
		labels = []
		for n in sorted(pos.keys(), key=lambda n: pos[n][1]):
			if time_groups[n] != last_group:
				time_group_lines.append((pos[n][1] + last_y)/2)
				last_group = time_groups[n]
				labels.append(last_group)
			last_y = pos[n][1]
		time_group_lines.append((pos[n][1] + last_y)/2)
		time_group_lines = time_group_lines[1:]
		time_group_lines[-1] += gap
		
		for i, y in enumerate(time_group_lines):
			pyplot.axhline(y, color='grey', linewidth=1, linestyle='--')
		
		last_y = y_min
		for i, y in enumerate(time_group_lines):
			pyplot.text(x_max + 0.6*step, (last_y + y)/2, "%s %d" % (time_groups_title, labels[i]), verticalalignment='center', horizontalalignment='left', fontsize=10)
			last_y = y
		
	pyplot.xlim(x_min - step/2, x_max + step/2)
	pyplot.gca().invert_yaxis()
	pyplot.tight_layout()
	pyplot.savefig("%s.png" % name)
	pyplot.close()


fdata = "vis_model_data.pickle"

if __name__ == '__main__':
	
	if not os.path.isfile(fdata):
		print("Loading data")
		model0 = Model(directory="_non_reduced/stage_0")
		model1 = Model(directory="_non_reduced/stage_1")
		model2 = Model(directory="_non_reduced/stage_2")
		model3 = Model(directory="_non_reduced/stage_3")
		
		samples = list(model0.samples.keys())
		
		eaps = dict([(i, model0.samples[s].excavation_area_phase) for i, s in enumerate(samples)])
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
		# earlier_than2, samples2 = model2.mphasing.create_earlier_than_matrix()
		# assert(samples == samples2)
		
		earlier_than3, samples3 = model3.mphasing.create_earlier_than_matrix()
		assert(samples == samples3)
		
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
	G3 = get_G(np.zeros(earlier_than0.shape, dtype = bool))
	
	et_by_dating = np.zeros(earlier_than0.shape, dtype = bool)
	for i, j in zip(*np.where(earlier_than2)):
		if (not earlier_than1[i, j]) and (groups[i] != groups[j]):
			et_by_dating[i,j] = True
	G_by_dating = get_G(et_by_dating)
	G2 = get_G(earlier_than1)
	for i, j in G_by_dating.edges():
		G2.add_edge(i, j)
	
	node_colors0 = 'k'
	node_colors_outlier = [('r' if i in outliers else 'k') for i in range(len(samples))]
	
	edge_colors2 = []
	for i, j in G2.edges():
		if et_by_dating[i, j]:
			color = 'r'
		else:
			color = 'k'
		edge_colors2.append(color)
	
	plot_graph("Stratigraphic Phasing", G0, pos, 'k', 'k', {}, {}, samples, [], None, None, 'g0')
	plot_graph("Outlier Detection", G1, pos, node_colors_outlier, 'k', dists0, dists0, samples, outliers, None, None, 'g1a')
	plot_graph("Chronological Modeling 1", G1, pos, 'k', 'k', dists0, dists1, samples, outliers, None, None, 'g1b')
	plot_graph("Inter-Group Chronological Relations", G2, pos, 'k', edge_colors2, dists0, dists1, samples, outliers, None, None, 'g2a')
	plot_graph("Chronological Modeling 2", G2, pos, 'k', 'lightgrey', dists1, dists2, samples, outliers, phases2, 'Phase', 'g2b')
	plot_graph("Chronological Clustering", G3, pos, 'k', 'k', dists1, dists2, samples, outliers, clusters, 'Cluster', 'g3a')
	plot_graph("Chronological Modeling 3", G3, pos, 'k', 'k', dists2, dists3, samples, outliers, phases3, 'Phase', 'g3b')
	
