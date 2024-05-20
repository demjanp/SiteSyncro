from typing import List, Dict

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import copy

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


def visualize_earlier_than(earlier_than: np.ndarray, samples: List[str]) -> None:
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


def reduce_earlier_than(earlier_than: np.ndarray) -> np.ndarray:
	# Create a directed graph from the earlier_than matrix
	G = nx.convert_matrix.from_numpy_array(earlier_than, create_using=nx.DiGraph)
	
	# Compute the transitive reduction
	reduced = nx.transitive_reduction(G)
	
	# Convert the reduced graph back to a matrix
	reduced_earlier_than = nx.convert_matrix.to_numpy_array(reduced)
	
	return reduced_earlier_than.astype(bool)


def find_groups(earlier_than: np.ndarray) -> Dict[int, List[int]]:
	if earlier_than.sum():
		G = nx.convert_matrix.from_numpy_array(earlier_than, create_using=nx.Graph)
		groups = []
		for c in nx.connected_components(G):
			G_sub = G.subgraph(c)
			groups.append(list(G_sub.nodes))
	else:
		groups = [np.arange(earlier_than.shape[0], dtype=int)]
	
	# groups = {group: [idx, ...], ...}; idx = index in earlier_than
	return dict(enumerate(sorted(groups, key=lambda group: len(group), reverse=True), start=1))


def eap_to_int(eap: str) -> float:
	"""
	Converts an excavation area phase name to an interval of integers.
	Args:
		eap (str): "1" or "1a" or "1-2" or "1a-b" or "1a-2b", higher = earlier (older) phase

	Returns:
		[[int, int], [int, int]]: [[major from, minor from], [major to, minor to]]

	"""
	
	def _name_to_int(name: str) -> [int, int]:
		# convert excavation area phase name to two numbers [major, minor]
		# 1 => [1,0]
		# 1a => [1,1]
	
		if not name:
			return None
		
		if not name[0].isdigit():
			return [None, ord(name) - ord("a") + 1]
		
		i = 0
		while i < len(name) and name[i].isdigit():
			i += 1
		if i == len(name):
			return [int(name), None]
		return [int(name[:i]), ord(name[i:]) - ord("a") + 1]
	
	if not eap:
		return None
		
	if "-" in eap:
		eap = eap.split("-")
		if len(eap) != 2:
			return None
		eap = [eap[0].strip(), eap[1].strip()]
	else:
		eap = [eap.strip(), eap.strip()]
	eap = [_name_to_int(eap[0]), _name_to_int(eap[1])]
	if eap[1][0] is None:
		eap[1][0] = eap[0][0]
	
	return eap


def get_phases_gr(earlier_than: np.ndarray, ranges_gr: List) -> np.ndarray:
	
	def _calc_dice(phasing, ranges_gr):
		"""
		Calculate the avg. dice coefficient of overlap between the sample dating ranges and combined ranges of their assigned phases
		Args:
			phasing: List of integers. The phasing values for each sample.
			ranges_gr: List of lists of floats. Each list contains the maximum and minimum values of the combined ranges for each sample.

		Returns:
			Float. The average dice coefficient of all samples
		"""
		
		def dice_coefficient(r1, r2, s1, s2):
			L_r = r2 - r1
			L_s = s2 - s1
			intersection_length = max(0, min(r2, s2) - max(r1, s1))
			dice_coef = (2 * intersection_length) / (L_r + L_s)
			return dice_coef
		
		ranges_found = dict([(ph, [-np.inf, np.inf]) for ph in set(phasing) if ph > -1])
		for i, ph in enumerate(phasing):
			if ph == -1:
				continue
			ranges_found[ph][0] = max(ranges_found[ph][0], ranges_gr[i][0])
			ranges_found[ph][1] = min(ranges_found[ph][1], ranges_gr[i][1])
		d_sum = 0
		for i, ph in enumerate(phasing):
			if ph == -1:
				continue
			r2, r1 = ranges_gr[i]
			s2, s1 = ranges_found[ph]
			d_sum += dice_coefficient(r1, r2, s1, s2)
		return d_sum / len(phasing)
	
	def _optimize_phasing(phasing_ranges, ranges_gr):
		"""
		Optimize phasing to maximize the overlap of dating ranges of individual samples with combined dating ranges of their phases
		Args:
			phasing_ranges: List of lists of integers. Each list contains the minimum and maximum phasing values for each sample.
			ranges_gr: List of lists of floats. Each list contains the maximum and minimum values of the grouped ranges for each sample.

		Returns:
			List of integers. The optimized phasing values for each sample.
		"""
		
		phase_max = max(rng[1] for rng in phasing_ranges)
		idxs_all = list(range(len(phasing_ranges)))
		idxs_moving = [i for i in idxs_all if phasing_ranges[i][0] != phasing_ranges[i][1]]
		phasing_opt = [-1]*len(phasing_ranges)
		for idx in idxs_all:
			if idx not in idxs_moving:
				phasing_opt[idx] = phasing_ranges[idx][0]
		while -1 in phasing_opt:
			ds0 = []
			res0 = []
			for i0 in idxs_moving:
				phasing0 = np.array(phasing_opt, dtype = int)
				for ph0 in range(phasing_ranges[i0][0], phasing_ranges[i0][1]+1):
					phasing0[i0] = ph0
					while (phasing0 == -1).any():
						ds1 = []
						res1 = []
						for i1 in np.where(phasing0 == -1)[0]:
							phasing1 = phasing0.copy()
							for ph1 in range(phasing_ranges[i1][0], phasing_ranges[i1][1]+1):
								phasing1[i1] = ph1
								ds1.append(_calc_dice(phasing1, ranges_gr))
								res1.append([i1, ph1])
						i1, ph1 = res1[np.argmax(ds1)]
						phasing0[i1] = ph1
					ds0.append(_calc_dice(phasing0, ranges_gr))
					res0.append(phasing0.tolist())
			if ds0:
				phasing_opt = res0[np.argmax(ds0)]
		d_opt = _calc_dice(phasing_opt, ranges_gr)
		changed = True
		while changed:
			changed = False
			for idx in idxs_moving:
				phs = list(range(phasing_ranges[idx][0], phasing_ranges[idx][1]+1))
				phasing = copy.copy(phasing_opt)
				ds = []
				for ph in phs:
					phasing[idx] = ph
					ds.append(_calc_dice(phasing, ranges_gr))
				d_max = max(ds)
				if d_max > d_opt:
					phasing_opt[idx] = phs[np.argmax(ds)]
					d_opt = d_max
					changed = True
		phasing_multi = {}
		for idx in idxs_moving:
			ph_collect = [phasing_opt[idx]]
			phasing = copy.copy(phasing_opt)
			for ph in range(phasing_ranges[idx][0], phasing_ranges[idx][1]+1):
				if ph in ph_collect:
					continue
				phasing[idx] = ph
				d = _calc_dice(phasing, ranges_gr)
				if d>=d_opt:
					ph_collect.append(ph)
			if len(ph_collect) > 1:
				phasing_multi[idx] = ph_collect
		return phasing_opt, phasing_multi
	
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
	
	phasing_ranges = []
	idxs_later = [np.where(earlier_than[idx])[0] for idx in range(earlier_than.shape[0])]
	idxs_earlier = [np.where(earlier_than[:,idx])[0] for idx in range(earlier_than.shape[0])]
	for idx in range(len(phasing)):
		phase_max = int(phasing.max())
		phase_min = 0
		ph_later = phasing[idxs_later[idx]]
		ph_later = ph_later[~np.isnan(ph_later)]
		if ph_later.size:
			phase_max = int(ph_later.min()) - 1
		ph_earlier = phasing[idxs_earlier[idx]]
		ph_earlier = ph_earlier[~np.isnan(ph_earlier)]
		if ph_earlier.size:
			phase_min = int(ph_earlier.max()) + 1
		phase_max = max(phase_min, phase_max)
		phasing_ranges.append([phase_min, phase_max])
	
	phasing, phasing_multi = _optimize_phasing(phasing_ranges, ranges_gr)
	
	return phasing, phasing_multi


def get_groups_and_phases(earlier_than: np.ndarray, samples: List[str], ranges: List) -> Dict[str, List[int or None]]:
	"""
	Determines the groups and phases for each sample based on the "earlier than" matrix.

	Parameters:
	earlier_than: matrix[n_samples x n_samples] = [True/False, ...]; sample in row is earlier than sample in column based on stratigraphy
	samples: A list of sample names.
	ranges: A list of sample ranges ordered by samples
	
	Returns:
	(groups_phases, phasing_multi)
		groups_phases={sample: [group, phase], ...}
		phasing_multi={sample: [phase, ...], ...}
	"""
	
	groups = find_groups(earlier_than)
	# groups = {group: [idx, ...], ...}; idx = index in earlier_than
	
	groups_phases = dict([(name, [None, None]) for name in samples])
	for gi in groups:
		for i in groups[gi]:
			groups_phases[samples[i]][0] = gi
	
	# Calculate phasing for each group
	phasing_multi = {}
	for gi in groups:
		earlier_than_gr = earlier_than[np.ix_(groups[gi], groups[gi])]
		samples_gr = [samples[i] for i in groups[gi]]
		ranges_gr = [ranges[i] for i in groups[gi]]
		phases_gr, phasing_multi_gr = get_phases_gr(earlier_than_gr, ranges_gr)
		for i in range(len(groups[gi])):
			groups_phases[samples_gr[i]][1] = int(phases_gr[i]) + 1
		for i in phasing_multi_gr:
			phasing_multi[samples_gr[i]] = phasing_multi_gr[i]
	
	return groups_phases, phasing_multi

