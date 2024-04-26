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

def create_earlier_than_matrix(model):
	
	samples = sorted(list(model.samples.keys()))
	
	# Create a matrix of earlier-than relationships
	earlier_than = np.zeros((len(samples), len(samples)), dtype=bool)
	for i, s1 in enumerate(samples):
		for s2 in model.samples[s1].earlier_than:
			j = samples.index(s2)
			earlier_than[i][j] = True
	
	# Update earlier-than relationships based on phases
	for i, s1 in enumerate(samples):
		for j, s2 in enumerate(samples):
			if s2 == s1:
				continue
			if model.samples[s1].area != model.samples[s2].area:
				continue
			if model.samples[s1].area_excavation_phase < model.samples[s2].area_excavation_phase:
				earlier_than[i][j] = True
	
	# Check if earlier_than has circular relationships
	if not check_circular_relationships(earlier_than, samples):
		# Visualize earlier_than as a DAG
		visualize_earlier_than(earlier_than, samples)
		raise Exception("Circular relationships detected")
		
	# Extend the earlier_than matrix to include computed relations
	earlier_than = extend_earlier_than(earlier_than)
	
	# earlier_than: matrix[n_samples x n_samples] = [True/False, ...]; sample in row is earlier than sample in column based on stratigraphy
	
	return earlier_than, samples

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
	# returns groups_phases = {sample: [group, phase], ...}
	
	groups = find_groups(earlier_than)
	# groups = {group: [idx, ...], ...}; idx = index in earlier_than
	
	groups_phases = dict([(name, [None, None]) for name in samples])
	for gi in groups:
		for i in groups[gi]:
			groups_phases[samples[i]][0] = gi
	
	# Calculate phasing for each group
	for gi in groups:
		earlier_than_gr = earlier_than[np.ix_(groups[gi], groups[gi])]
		samples_gr = [samples[i] for i in groups[gi]]
		phases_gr = get_phases_gr(earlier_than_gr)
		for i in range(len(groups[gi])):
			groups_phases[samples_gr[i]][1] = int(phases_gr[i]) + 1
	
	return groups_phases

def update_earlier_than_by_clustering(model, earlier_than, samples):
	# Update earlier_than based on temporal clustering of the samples
	
	if model.clusters_opt_n is None:
		raise Exception("No clustering found")
	
	clusters = model.clusters[model.clusters_opt_n]
	
	# Sort clusters from oldest to youngest
	labels = sorted(clusters.keys(), key=lambda label: means[label], reverse=True)
	
	phases_clu = dict((label, idx + 1) for idx, label in enumerate(labels))
	# phases_clu = {label: phase, ...}; lower phase = earlier
	
	collect = {}
	for label in phases_clu:
		for sample in clusters[label]:
			collect[sample] = phases_clu[label]
	phases_clu = collect
	# phases_clu = {sample: phase, ...}
	
	# Update earlier-than relationships based on phases derived from clustering
	errors = []
	for i, s1 in enumerate(samples):
		for j, s2 in enumerate(samples):
			if phases_clu[s1] < phases_clu[s2]:
				earlier_than[i][j] = True
				if (model.samples[s1].group == model.samples[s2].group) and (model.samples[s1].phase > model.samples[s2].phase):
					errors.append([s1, s2, phases_clu[s1], phases_clu[s2], model.samples[s1].phase, model.samples[s2].phase])
	if errors:
		print("Warning, collisions detected between stratigraphic phasing and clustering:")
		for s1, s2, clu1, clu2, ph1, ph2 in errors:
			print("%s (Strat. phase %s, Clu. phase %s), %s (Strat. phase %s, Clu. phase %s)" % (s1, ph1, clu1, s2, ph2, clu2))
		print()
	
	# Check if earlier_than has circular relationships
	if not check_circular_relationships(earlier_than, samples):
		# Visualize earlier_than as a DAG
		visualize_earlier_than(earlier_than, samples)
		raise Exception("Circular relationships detected")
	
	# Extend the earlier_than matrix to include computed relations
	earlier_than = extend_earlier_than(earlier_than)
	
	return earlier_than

