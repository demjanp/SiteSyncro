import matplotlib.pyplot as plt
import networkx as nx
import numpy as np

def check_circular_relationships(earlier_than, samples):
	
	G = nx.convert_matrix.from_numpy_array(earlier_than, create_using=nx.DiGraph)
	cycles = list(nx.simple_cycles(G))
	if cycles:
		print("Circular relationships detected:")
		for cycle in cycles:
			cycle_samples = [samples[i] for i in cycle]
			print(" -> ".join(cycle_samples))
		return False
	
	print("No circular relationships detected.")
	return True

def visualize_earlier_than(earlier_than, samples):
	
	G = nx.convert_matrix.from_numpy_array(earlier_than, create_using=nx.DiGraph)
	labels = {i: sample for i, sample in enumerate(samples)}
	pos = nx.spring_layout(G)
	nx.draw(G, pos, labels=labels, with_labels=True)
	plt.show()

def extend_earlier_than(earlier_than):
	# Create a directed graph from the earlier_than matrix
	G = nx.convert_matrix.from_numpy_array(earlier_than, create_using=nx.DiGraph)

	# Compute the transitive closure of the graph
	transitive_closure = nx.transitive_closure(G)

	# Convert the transitive closure graph back to a matrix
	extended_earlier_than = nx.convert_matrix.to_numpy_array(transitive_closure)

	return extended_earlier_than.astype(bool)

def find_groups(earlier_than):
	
	G = nx.convert_matrix.from_numpy_array(earlier_than, create_using = nx.Graph)
	groups = []
	for c in nx.connected_components(G):
		G_sub = G.subgraph(c)
		groups.append(list(G_sub.nodes))
	
	# groups = {group: [idx, ...], ...}; idx = index in earlier_than
	return dict(enumerate(sorted(groups, key = lambda group: len(group), reverse = True), start = 1))

def create_earlier_than_matrix(context_phase, earlier_than_rel, samples, context_samples, context_area):
	
	# Create a matrix of earlier-than relationships
	earlier_than = np.zeros((len(samples), len(samples)), dtype=bool)
	for s1 in earlier_than_rel:
		i = samples.index(s1)
		for s2 in earlier_than_rel[s1]:
			j = samples.index(s2)
			earlier_than[i][j] = True
	
	# Update earlier-than relationships based on phases
	for c1 in context_samples:
		a1 = context_area[c1]
		for c2 in context_samples:
			a2 = context_area[c2]
			if (a1 != a2):
				continue
			for sample1 in context_samples[c1]:
				i = samples.index(sample1)
				for sample2 in context_samples[c2]:
					j = samples.index(sample2)
					if (context_phase[c1] < context_phase[c2]):
						earlier_than[i][j] = True
					elif (context_phase[c2] < context_phase[c1]):
						earlier_than[j][i] = True
	
	# Check if earlier_than has circular relationships
	if not check_circular_relationships(earlier_than, samples):
		# Visualize earlier_than as a DAG
		visualize_earlier_than(earlier_than, samples)
		raise Exception("Circular relationships detected")
		
	# Extend the earlier_than matrix to include computed relations
	earlier_than = extend_earlier_than(earlier_than)
	
	# earlier_than: matrix[n_samples x n_samples] = [True/False, ...]; sample in row is earlier than sample in column based on stratigraphy
	
	return earlier_than

def update_earlier_than_matrix(earlier_than, samples, clusters, means):
	# Update earlier_than based on temporal clustering of the samples
	#
	# clusters = {label: [sample, ...], ...}
	# means = {label: mean, ...}
	
	# Sort clusters from oldest to youngest
	labels = sorted(clusters.keys(), key=lambda label: means[label], reverse=True)
	
	phases_clu = dict((label, idx + 1) for idx, label in enumerate(labels))
	# phases_clu = {label: phase, ...}; lower phase = earlier
	
	phases = {}
	for label in phases_clu:
		for sample in clusters[label]:
			phases[sample] = phases_clu[label]
	# phases = {sample: phase, ...}
	
	# Update earlier-than relationships based on phases derived from clustering
	for i, s1 in enumerate(samples):
		for j, s2 in enumerate(samples):
			if phases[s1] < phases[s2]:
				earlier_than[i][j] = True
			elif phases[s2] < phases[s1]:
				earlier_than[j][i] = True
	
	# Check if earlier_than has circular relationships
	if not check_circular_relationships(earlier_than, samples):
		# Visualize earlier_than as a DAG
		visualize_earlier_than(earlier_than, samples)
		raise Exception("Circular relationships detected")
	
	# Extend the earlier_than matrix to include computed relations
	earlier_than = extend_earlier_than(earlier_than)
	
	return earlier_than

def get_phases_gr(earlier_than):
	n_samples = earlier_than.shape[0]
	phasing = np.full(n_samples, np.nan)  # phasing[si] = phase; lower = earlier
	
	# assign phase to samples latest to earliest
	mask_todo = earlier_than.copy()
	phase = 0
	while mask_todo.any():
		latest = (mask_todo.any(axis=0) & ~mask_todo.any(axis=1))
		phasing[latest] = phase
		mask_todo[:, latest] = False
		phase += 1
	
	# assign phases to samples earliest to latest, if not already assigned
	mask_todo = earlier_than.copy()
	phase = n_samples
	while mask_todo.any():
		earliest = (mask_todo.any(axis=1) & ~mask_todo.any(axis=0))
		phasing[np.isnan(phasing) & earliest] = phase
		mask_todo[earliest] = False
		phase -= 1
	
	# minimize range of phases
	vals = np.unique(phasing[~np.isnan(phasing)])
	vals.sort()
	collect = phasing.copy()
	for val_new, val in enumerate(vals):
		collect[phasing == val] = val_new
	phasing = collect
	
	mask = (~np.isnan(phasing))
	if mask.any():
		phasing[mask] = phasing[mask].max() - phasing[mask]
	phasing[~mask] = -1
	return phasing

def get_groups_and_phases(earlier_than, samples):
	
	groups = find_groups(earlier_than)
	# groups = {group: [idx, ...], ...}; idx = index in earlier_than
	
	# Calculate phasing for each group
	phases = {}
	for gi in groups:
		earlier_than_gr = earlier_than[np.ix_(groups[gi], groups[gi])]
		samples_gr = [samples[i] for i in groups[gi]]
		phases_gr = get_phases_gr(earlier_than_gr)
		phases[gi] = dict([(samples_gr[i], int(phases_gr[i]) + 1) for i in range(len(groups[gi]))])
	groups = dict([(gi, [samples[i] for i in groups[gi]]) for gi in groups])
	
	# groups = {group: [sample, ...], ...}
	# phases = {group: {sample: phase, ...}, ...}
	return groups, phases

