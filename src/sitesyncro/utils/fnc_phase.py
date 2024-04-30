from sitesyncro.utils.fnc_stat import (calc_sum, samples_to_distributions)

from itertools import combinations
import matplotlib.pyplot as plt
import networkx as nx
from tqdm import tqdm
import numpy as np
import math

from typing import List, Dict

def check_circular_relationships(earlier_than: np.ndarray, samples: List[str]) -> bool:
	G = nx.convert_matrix.from_numpy_array(earlier_than, create_using=nx.DiGraph)
	cycles = list(nx.simple_cycles(G))
	if cycles:
		print("Circular relationships detected:")
		for cycle in cycles:
			cycle_samples = [samples[i] for i in cycle]
			print(" -> ".join(cycle_samples))
		return False
	return True

def visualize_earlier_than(earlier_than: np.ndarray, samples: List[str]):
	G = nx.convert_matrix.from_numpy_array(earlier_than, create_using=nx.DiGraph)
	labels = {i: sample for i, sample in enumerate(samples)}
	pos = nx.spring_layout(G)
	nx.draw(G, pos, labels=labels, with_labels=True)
	plt.show()

def extend_earlier_than(earlier_than: np.ndarray) -> np.ndarray:
	# Create a directed graph from the earlier_than matrix
	G = nx.convert_matrix.from_numpy_array(earlier_than, create_using=nx.DiGraph)

	# Compute the transitive closure of the graph
	transitive_closure = nx.transitive_closure(G)

	# Convert the transitive closure graph back to a matrix
	extended_earlier_than = nx.convert_matrix.to_numpy_array(transitive_closure)

	return extended_earlier_than.astype(bool)

def find_groups(earlier_than: np.ndarray) -> Dict[int, List[int]]:
	
	if earlier_than.sum():
		G = nx.convert_matrix.from_numpy_array(earlier_than, create_using = nx.Graph)
		groups = []
		for c in nx.connected_components(G):
			G_sub = G.subgraph(c)
			groups.append(list(G_sub.nodes))
	else:
		groups = [np.arange(earlier_than.shape[0], dtype = int)]
	
	# groups = {group: [idx, ...], ...}; idx = index in earlier_than
	return dict(enumerate(sorted(groups, key = lambda group: len(group), reverse = True), start = 1))

def create_earlier_than_matrix(model: object):
	samples = sorted(list(model.samples.keys()))
	
	# Create a matrix of earlier-than relationships
	earlier_than = np.zeros((len(samples), len(samples)), dtype=bool)
	for i, s1 in enumerate(samples):
		for s2 in model.samples[s1].earlier_than:
			j = samples.index(s2)
			earlier_than[i][j] = True
	
	# Update earlier-than relationships based on phases
	for i, s1 in enumerate(samples):
		if model.samples[s1].excavation_area_phase is None:
			continue
		for j, s2 in enumerate(samples):
			if s2 == s1:
				continue
			if model.samples[s2].excavation_area_phase is None:
				continue
			if model.samples[s1].area != model.samples[s2].area:
				continue
			if model.samples[s1].excavation_area_phase < model.samples[s2].excavation_area_phase:
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

def get_phases_gr(earlier_than: np.ndarray):
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

def get_groups_and_phases(earlier_than: np.ndarray, samples: List[str]) -> Dict[str, List[int]]:
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


def prob_earlier_than(dist1: np.ndarray, dist2: np.ndarray):
	# Calculate the probability that dist1 is earlier than dist2
	#
	# Returns p
	
	p = 1 - np.sum(np.cumsum(dist2[::-1])[::-1] * np.cumsum(dist1))
	return p

def update_earlier_than_by_dating(model: object, earlier_than: np.ndarray, samples: List[str]) -> np.ndarray:
	# Update earlier_than based on probability distributions of the samples
	
	distributions, samples_d, joined = samples_to_distributions(model.samples.values())
	# distributions = [[p, ...], ...]
	# samples = [combined_name, ...] ordered by distributions
	# joined = {combined_name: [sample name, ...], ...}
	
	sample_lookup = {}
	for combined_name in joined:
		for sample in joined[combined_name]:
			sample_lookup[combined_name] = sample
	for sample in samples_d:
		if sample not in sample_lookup:
			sample_lookup[sample] = sample
	
	# Update earlier_than based on probability distributions
	for i, s1 in enumerate(samples_d):
		for j, s2 in enumerate(samples_d):
			if i == j:
				continue
			if prob_earlier_than(distributions[i], distributions[j]) >= (1 - model.p_value):
				earlier_than[samples.index(sample_lookup[s1])][samples.index(sample_lookup[s2])] = True
	
	# Extend the earlier_than matrix to include computed relations
	earlier_than = extend_earlier_than(earlier_than)
	
	return earlier_than
	
def update_earlier_than_by_clustering(model: object, earlier_than: np.ndarray, samples: List[str]) -> np.ndarray:
	# Update earlier_than based on temporal clustering of the samples
	
	if model.cluster_opt_n is None:
		raise Exception("No clustering found")
	
	clusters = model.clusters[model.cluster_opt_n]
	means = model.cluster_means[model.cluster_opt_n]
	
	if len(clusters) < 2:
		return earlier_than
	
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
				if (model.samples[s1].group == model.samples[s2].group) and (model.samples[s1].phase is not None) and (model.samples[s2].phase is not None) and (model.samples[s1].phase > model.samples[s2].phase):
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

def find_dating_outliers(model: object):
	
	def _find_outliers_idxs(idxs: set, earlier_than: np.ndarray, samples: List[str], ranges: Dict[str, List[float]], check_only: bool = False):
		found = set()
		for i in idxs:
			for j in np.where(earlier_than[i])[0]:
				if ranges[samples[i]][0] < ranges[samples[j]][1]:
					if check_only:
						return True
					found.add(i)
					found.add(j)
		return list(found)
	
	def _pick_outliers(candidates: List[int], sample_idxs: List[int], earlier_than: np.ndarray, samples: List[str], ranges: Dict[str, List[float]]):
		removed = set(candidates)
		base = set(sample_idxs).difference(removed)
		# Try adding back individual outliers first
		addable = set()
		for i in removed:
			if not _find_outliers_idxs(base.union({i}), earlier_than, samples, ranges, check_only = True):
				addable.add(i)
		# Try adding back increasingly larger groups of outliers
		outliers = []
		added_n = -1
		with tqdm() as pbar:
			for n in range(1, len(addable) + 1):
				found = False
				n_combs = math.comb(len(addable), n)
				pbar.reset()
				pbar.total = n_combs
				pbar.refresh()
				pbar.set_description("Returning %d/%d" % (n, len(addable)))
				picked = []
				for added in combinations(addable, n):
					pbar.update(1)
					if not _find_outliers_idxs(base.union(added), earlier_than, samples, ranges, check_only = True):
						picked.append(added)
				for added in picked:
					found = True
					if n > added_n:
						added_n = n
						outliers = []
					outliers.append(list(removed.difference(added)))
				if not found:
					return outliers
		return outliers
	
	earlier_than, samples = create_earlier_than_matrix(model)
	if not earlier_than.sum():
		return [], []
	ranges = {}
	for name in model.samples:
		if model.samples[name].outlier:
			continue
		rng = model.samples[name].get_range()
		if None not in rng:
			ranges[name] = rng
	sample_idxs = [idx for idx, name in enumerate(samples) if name in ranges]
	if not sample_idxs:
		return [], []
	candidates = _find_outliers_idxs(sample_idxs, earlier_than, samples, ranges)
	if not candidates:
		return [], []
	print("Found %d candidates for outliers" % (len(candidates)))
	outliers = _pick_outliers(candidates, sample_idxs, earlier_than, samples, ranges)
	if outliers:
		candidates = set()
		for row in outliers:
			candidates.update(row)
		candidates = list(candidates)
		outliers = max(outliers, key = lambda row: sum([int(model.samples[samples[i]].redeposited) for i in row]))
	outliers = sorted([samples[i] for i in outliers])
	candidates = sorted([samples[i] for i in candidates])
	for name in model.outliers:
		if name not in candidates:
			candidates.append(name)
	print("Found %d outliers" % (len(outliers)))
	
	return outliers, candidates
	
