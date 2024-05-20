import math
from itertools import combinations
from typing import List, Dict

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
from tqdm import tqdm

from sitesyncro.utils.fnc_stat import (samples_to_distributions, calc_range)
from sitesyncro.utils.fnc_phase import (check_circular_relationships, visualize_earlier_than, extend_earlier_than,
                                        eap_to_int, get_groups_and_phases)

class MPhasing(object):
	"""
	A class implementing phasing functionality
	
	:param model: The parent Model object
	:type model: Model
	"""
	def __init__(self, model: 'Model'):
		
		self.model = model
	
	def create_earlier_than_matrix(self) -> (np.ndarray, List[str]):
		"""
		Creates a matrix representing the "earlier than" relationships between samples.
		
		Returns:
		(earlier_than, samples)
			- earlier_than: matrix[n_samples x n_samples] = [True/False, ...]; sample in row is earlier than sample in column based on stratigraphy
			- samples: The list of sample names.
		"""
		
		def _is_earlier(eap1: [[int, int], [int, int]], eap2: [[int, int], [int, int]]):
			# Check if eap1 is earlier than eap2
			# eap = [[major from, minor from], [major to, minor to]]
			
			if eap1 is None:
				return False
			if eap2 is None:
				return False
			
			eap1 = eap_to_int(eap1)
			eap2 = eap_to_int(eap2)
			
			for val in [eap1, eap2]:
				if val is None:
					return False
				if val[0] is None:
					return False
			
			if (eap1[0][1] is None) or (eap2[1][1] is None):
				return (eap1[0][0] > eap2[1][0])
			return (eap1[0] > eap2[1])
		
		samples = sorted(list(self.model.samples.keys()))
		
		# Create a matrix of earlier-than relationships
		earlier_than = np.zeros((len(samples), len(samples)), dtype=bool)
		for i, s1 in enumerate(samples):
			for s2 in self.model.samples[s1].earlier_than:
				j = samples.index(s2)
				earlier_than[j][i] = True
		
		# Update earlier-than relationships based on excavation area phases
		for i, s1 in enumerate(samples):
			if self.model.samples[s1].excavation_area_phase is None:
				continue
			for j, s2 in enumerate(samples):
				if s2 == s1:
					continue
				if self.model.samples[s1].area != self.model.samples[s2].area:
					continue
				earlier_than[i][j] = _is_earlier(self.model.samples[s1].excavation_area_phase, self.model.samples[s2].excavation_area_phase)
		
		# Update earlier-than relationships based on groups and phases
		for i, s1 in enumerate(samples):
			if self.model.samples[s1].phase is None:
				continue
			for j, s2 in enumerate(samples):
				if s2 == s1:
					continue
				if self.model.samples[s2].group is None:
					continue
				if self.model.samples[s1].group != self.model.samples[s2].group:
					continue
				if self.model.samples[s1].phase < self.model.samples[s2].phase:
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
	
	def update_earlier_than_by_dating(self, earlier_than: np.ndarray, samples: List[str]) -> np.ndarray:
		"""
		Updates the "earlier than" matrix based on the probability distributions of the samples.
		
		A sample is considered earlier than another sample if the probability of its dating being earlier is greater than 1 - p_value.

		Parameters:
		earlier_than: matrix[n_samples x n_samples] = [True/False, ...]; sample in row is earlier than sample in column based on stratigraphy
		samples: The list of sample names.

		Returns:
		np.ndarray: The updated "earlier than" matrix.
		"""
		
		earlier_than = earlier_than.copy()
		distributions, names_dist, joined = samples_to_distributions(self.model.samples.values())
		
		# distributions = [[p, ...], ...]
		# names_dist = [combined_name, ...] ordered by distributions
		# joined = {combined_name: [sample name, ...], ...}
		
		def _get_combined(name, joined):
			for combined_name in joined:
				if name in joined[combined_name]:
					return combined_name
			return None
		
		idx_lookup = {}
		for name in samples:
			idx_lookup[name] = None
			if name in names_dist:
				idx_lookup[name] = names_dist.index(name)
			else:
				combined_name = _get_combined(name, joined)
				if combined_name is not None:
					idx_lookup[name] = names_dist.index(combined_name)
		
		ranges = {}
		for i in range(len(distributions)):
			ranges[i] = calc_range(self.model.years, distributions[i])
		
		# Update earlier_than based on probability distributions
		for i, s1 in enumerate(samples):
			d_i = idx_lookup[s1]
			if d_i is None:
				continue
			for j, s2 in enumerate(samples):
				if i == j:
					continue
				d_j = idx_lookup[s2]
				if d_j is None:
					continue
				if ranges[d_i][1] > ranges[d_j][0]:
					earlier_than[i][j] = True
		
		# Extend the earlier_than matrix to include computed relations
		earlier_than = extend_earlier_than(earlier_than)
		
		return earlier_than
	
	def update_earlier_than_by_clustering(self, earlier_than: np.ndarray, samples: List[str]) -> np.ndarray:
		"""
		Updates the "earlier than" matrix based on the clustering of the samples.
		
		A sample is considered earlier than another sample if it belongs to a cluster with an earlier mean dating.
		
		Parameters:
		earlier_than: matrix[n_samples x n_samples] = [True/False, ...]; sample in row is earlier than sample in column based on stratigraphy
		samples: The list of sample names.
		
		Returns:
		np.ndarray: The updated "earlier than" matrix.
		"""
		# Update earlier_than based on temporal clustering of the samples
		
		if self.model.cluster_opt_n is None:
			raise Exception("No clustering found")
		
		earlier_than = earlier_than.copy()
		clusters = self.model.clusters[self.model.cluster_opt_n]
		means = self.model.cluster_means[self.model.cluster_opt_n]
		
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
					if (self.model.samples[s1].group == self.model.samples[s2].group) and (self.model.samples[s1].phase is not None) and (
							self.model.samples[s2].phase is not None) and (self.model.samples[s1].phase > self.model.samples[s2].phase):
						errors.append(
							[s1, s2, phases_clu[s1], phases_clu[s2], self.model.samples[s1].phase, self.model.samples[s2].phase])
		if errors:
			print("Warning, collisions detected between stratigraphic phasing and clustering:")
			for s1, s2, clu1, clu2, ph1, ph2 in errors:
				print("%s (Strat. phase %s, Clu. phase %s), %s (Strat. phase %s, Clu. phase %s)" % (
					s1, ph1, clu1, s2, ph2, clu2))
			print()
		
		# Check if earlier_than has circular relationships
		if not check_circular_relationships(earlier_than, samples):
			# Visualize earlier_than as a DAG
			visualize_earlier_than(earlier_than, samples)
			raise Exception("Circular relationships detected")
		
		# Extend the earlier_than matrix to include computed relations
		earlier_than = extend_earlier_than(earlier_than)
		
		return earlier_than
	
	def find_dating_outliers(self) -> (List[str], List[str]):
		"""
		Finds outliers based on the stratigraphic relationships between samples.
		
		Returns:
		(outliers, candidates)
			- outliers: The list of outlier sample names selected from the candidates so that the amount of redeposited samples is maximized.
			- candidates: The list of candidate sample names.
		"""
		
		def _find_outliers_idxs(idxs: list or set, earlier_than: np.ndarray, samples: List[str],
								ranges: Dict[str, List[float]],
								check_only: bool = False):
			found = set()
			for i in idxs:
				for j in np.where(earlier_than[i])[0]:
					if ranges[samples[i]][0] < ranges[samples[j]][1]:
						if check_only:
							return True
						found.add(i)
						found.add(j)
			return list(found)
		
		def _pick_outliers(candidates: List[int], sample_idxs: List[int], earlier_than: np.ndarray, samples: List[str],
						   ranges: Dict[str, List[float]]):
			removed = set(candidates)
			base = set(sample_idxs).difference(removed)
			# Try adding back individual outliers first
			addable = set()
			for i in removed:
				if not _find_outliers_idxs(base.union({i}), earlier_than, samples, ranges, check_only=True):
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
						if not _find_outliers_idxs(base.union(added), earlier_than, samples, ranges, check_only=True):
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
		
		earlier_than, samples = self.create_earlier_than_matrix()
		if not earlier_than.sum():
			return [], []
		ranges = {}
		for name in self.model.samples:
			if self.model.samples[name].outlier:
				continue
			rng = self.model.samples[name].get_range()
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
			outliers = max(outliers, key=lambda row: sum([int(self.model.samples[samples[i]].redeposited) for i in row]))
		outliers = sorted([samples[i] for i in outliers])
		candidates = sorted([samples[i] for i in candidates])
		for name in self.model.outliers:
			if name not in candidates:
				candidates.append(name)
		print("Found %d outliers" % (len(outliers)))
		
		return outliers, candidates
	
	def process(self, by_clusters: bool = False, by_dates: bool = False, max_cpus: int = -1, max_queue_size: int = 10000, batch_size: int = 10000) -> bool:
		"""
		Updates the phasing of samples based on stratigraphic relations.

		Updates the groups and phases of samples based on their stratigraphic relations.

		:param by_clusters: If True, update the phasing by clustering sample dates. Defaults to False.
		:type by_clusters: bool, optional
		:param by_dates: If True, update the phasing by comparing sample dates. Defaults to False.
		:type by_dates: bool, optional
		:param max_cpus: The maximum number of CPUs to use for multiprocessing. Defaults to -1, which means all available CPUs will be used.
		:type max_cpus: int, optional
		:param max_queue_size: The maximum size of the queue for multiprocessing. -1 means the queue size is unlimited. Defaults to 10000.
		:type max_queue_size: int, optional
		:param batch_size: If set to >0, process data in batches. Higher values speed up processing but use more memory. Defaults to 10000.
		:type batch_size: int, optional
		
		:return: True if phasing has changed, False otherwise.
		:rtype: bool
		"""
		
		earlier_than, samples = self.create_earlier_than_matrix()
		if by_clusters and self.model.is_clustered:
			earlier_than = self.update_earlier_than_by_clustering(earlier_than, samples)
		
		if by_dates and self.model.is_modeled:
			earlier_than = self.update_earlier_than_by_dating(earlier_than, samples)
		
		ranges = [self.model.samples[name].get_range() for name in self.model.samples]		
		groups_phases = get_groups_and_phases(earlier_than, samples, ranges, max_cpus, max_queue_size, batch_size)
		
		# groups_phases = {sample: [group, phase], ...}
		for name in self.model.samples:
			if name in groups_phases:
				group, phase = groups_phases[name]
				self.model.samples[name].set_group(group)
				self.model.samples[name].set_phase(phase)
			else:
				self.model.samples[name].set_group(None)
				self.model.samples[name].set_phase(None)
		
		earlier_than_1, samples_1 = self.create_earlier_than_matrix()
		if not ((earlier_than.shape == earlier_than.shape) and np.allclose(earlier_than, earlier_than_1) and (
				samples == samples_1)):
			return True
		return False
