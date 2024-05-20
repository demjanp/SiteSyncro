from sitesyncro.utils.fnc_mp import (process_mp)

from typing import List, Dict

import matplotlib.pyplot as plt
from itertools import product
from functools import reduce
from tqdm import tqdm
import networkx as nx
import numpy as np
import operator

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


def calc_rng(phasing, ranges_gr):
	"""
	Calculate the range of the grouped ranges for the given phasing values.
	Args:
		phasing: List of integers. The phasing values for each sample.
		ranges_gr: List of lists of floats. Each list contains the maximum and minimum values of the grouped ranges for each sample.

	Returns:
		Float. The range of the grouped ranges.
	"""
	
	ranges_found = dict([(ph, [-np.inf, np.inf]) for ph in set(phasing)])
	for i, ph in enumerate(phasing):
		ranges_found[ph][0] = max(ranges_found[ph][0], ranges_gr[i][0])
		ranges_found[ph][1] = min(ranges_found[ph][1], ranges_gr[i][1])
	rng = sum((ranges_found[ph][0] - ranges_found[ph][1]) for ph in ranges_found) / len(ranges_found)
	return rng


def worker_fnc(phasing: List, ranges_gr: List):
	
	rng = calc_rng(phasing, ranges_gr)
	return rng, phasing


def collect_fnc(data: List, rng_opt: List, phasing_opt: List, pbar: tqdm):
	pbar.update(len(data))
	for rng, phasing in data:
		if rng < rng_opt[0]:
			rng_opt[0] = rng
			phasing_opt[0] = phasing


def get_phases_gr(earlier_than: np.ndarray, ranges_gr: List, g_cnt, g_cmax) -> np.ndarray:
	
	def _reduce_phasing_sp(values, ranges_gr):
		"""
		Optimize phasing ranges to minimize the range of the grouped ranges.
		Args:
			phasing_ranges: List of lists of integers. Each list contains the minimum and maximum phasing values for each sample.
			ranges_gr: List of lists of floats. Each list contains the maximum and minimum values of the grouped ranges for each sample.

		Returns:
			List of integers. The optimized phasing values for each sample.
		"""
		
		phasing_opt = None
		rng_opt = np.inf
		for phasing in product(*values):
			rng = calc_rng(phasing, ranges_gr)
			if rng < rng_opt:
				rng_opt = rng
				phasing_opt = phasing
		phasing = np.array(list(phasing_opt))
		return phasing
	
	def _reduce_phasing(phasing_ranges, ranges_gr):
		"""
		Optimize phasing ranges to minimize the range of the grouped ranges.
		Args:
			phasing_ranges: List of lists of integers. Each list contains the minimum and maximum phasing values for each sample.
			ranges_gr: List of lists of floats. Each list contains the maximum and minimum values of the grouped ranges for each sample.

		Returns:
			List of integers. The optimized phasing values for each sample.
		"""
		
		phasing_opt = None
		rng_opt = np.inf
		values = [list(range(phase_min, phase_max+1)) for phase_min, phase_max in phasing_ranges]
		rng_opt = [np.inf]
		phasing_opt = [None]
		cnt = [0]
		cmax = reduce(operator.mul, (len(value) for value in values), 1)
		if cmax < 100000:
			return _reduce_phasing_sp(values, ranges_gr)
		with tqdm(total=cmax) as pbar:
			pbar.set_description("Optimizing phasing for group %d/%d" % (g_cnt, g_cmax))
			process_mp(worker_fnc, product(*values), [ranges_gr], collect_fnc = collect_fnc, collect_args = [rng_opt, phasing_opt, pbar], max_queue_size = 10000, batch_size = 10000)
		return np.array(list(phasing_opt[0]))
	
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
	
	idxs_later = [np.where(earlier_than[idx])[0] for idx in range(earlier_than.shape[0])]
	idxs_earlier = [np.where(earlier_than[:,idx])[0] for idx in range(earlier_than.shape[0])]
	
	phasing_ranges = []
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
	
	# Iterate over all possible combinations of phasing and find the one with the smallest combined dating range per phase
	phasing = _reduce_phasing(phasing_ranges, ranges_gr)
	
	return phasing


def get_groups_and_phases(earlier_than: np.ndarray, samples: List[str], ranges: List) -> Dict[str, List[int or None]]:
	"""
	Determines the groups and phases for each sample based on the "earlier than" matrix.

	Parameters:
	earlier_than: matrix[n_samples x n_samples] = [True/False, ...]; sample in row is earlier than sample in column based on stratigraphy
	samples: A list of sample names.
	ranges: A list of sample ranges ordered by samples

	Returns:
	groups_phases: {sample: [group, phase], ...}
	"""
	
	groups = find_groups(earlier_than)
	# groups = {group: [idx, ...], ...}; idx = index in earlier_than
	
	groups_phases = dict([(name, [None, None]) for name in samples])
	for gi in groups:
		for i in groups[gi]:
			groups_phases[samples[i]][0] = gi
	
	# Calculate phasing for each group
	cnt = 1
	cmax = len(groups)
	for gi in groups:
		earlier_than_gr = earlier_than[np.ix_(groups[gi], groups[gi])]
		samples_gr = [samples[i] for i in groups[gi]]
		ranges_gr = [ranges[i] for i in groups[gi]]
		phases_gr = get_phases_gr(earlier_than_gr, ranges_gr, cnt, cmax)
		for i in range(len(groups[gi])):
			groups_phases[samples_gr[i]][1] = int(phases_gr[i]) + 1
		cnt += 1
	
	return groups_phases

