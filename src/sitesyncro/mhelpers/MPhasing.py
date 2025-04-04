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
		
		def _is_earlier_eap(eap1: [[int, int], [int, int]], eap2: [[int, int], [int, int]]):
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
		
		def _is_earlier_site_phase(ph1: [[int, int], [int, int]], ph2: [[int, int], [int, int]]):
			# Check if ph1 is earlier than ph2
			# ph = [[major from, minor from], [major to, minor to]]
			
			if ph1 is None:
				return False
			if ph2 is None:
				return False
			
			ph1 = eap_to_int(ph1)
			ph2 = eap_to_int(ph2)
			
			for val in [ph1, ph2]:
				if val is None:
					return False
				if val[0] is None:
					return False
			
			if (ph1[0][1] is None) or (ph2[1][1] is None):
				return (ph1[0][0] < ph2[1][0])
			return (ph1[0] < ph2[1])
		
		samples = sorted(list(self.model.samples.keys()))
		
		# Create a matrix of earlier-than relationships
		earlier_than = np.zeros((len(samples), len(samples)), dtype=bool)
		for i, s1 in enumerate(samples):
			for s2 in self.model.samples[s1].earlier_than:
				j = samples.index(s2)
#				earlier_than[i][j] = True
				earlier_than[j][i] = True  # DEBUG temp. fix to accomodate incorectly used Earlier-than column in input data
		
		# Check earlier-than relationships based on excavation area phases for consistency with stratigraphic relationships
		errors = []
		earlier_than_tst = extend_earlier_than(earlier_than)
		for i, s1 in enumerate(samples):
			if self.model.samples[s1].excavation_area_phase is None:
				continue
			for j, s2 in enumerate(samples):
				if s2 == s1:
					continue
				if self.model.samples[s1].area != self.model.samples[s2].area:
					continue
				is_e = _is_earlier_eap(self.model.samples[s1].excavation_area_phase, self.model.samples[s2].excavation_area_phase)
				if (is_e and earlier_than_tst[j][i]):
					errors.append((s1, s2))
		if errors:
			print()
			print("Excavation Area Phases not consistent with stratigraphic relationships:")
			for s1, s2 in errors:
				print("%s (EAP: %s) strat. earlier than %s (EAP: %s)" % (s2, self.model.samples[s2].excavation_area_phase, s1, self.model.samples[s1].excavation_area_phase))
			print()
			raise Exception("Modeling aborted: Inconsistent stratigraphic relations detected.")
		
		# Check earlier-than relationships based on site phases for consistency with stratigraphic relationships
		errors = []
		earlier_than_tst = extend_earlier_than(earlier_than)
		for i, s1 in enumerate(samples):
			if self.model.samples[s1].site_phase is None:
				continue
			for j, s2 in enumerate(samples):
				if s2 == s1:
					continue
				is_e = _is_earlier_site_phase(self.model.samples[s1].site_phase, self.model.samples[s2].site_phase)
				if (is_e and earlier_than_tst[j][i]):
					errors.append((s1, s2))
		if errors:
			print()
			print("Site Phases not consistent with stratigraphic relationships:")
			for s1, s2 in errors:
				print("%s (Phase: %s) strat. earlier than %s (Phase: %s)" % (s2, self.model.samples[s2].site_phase, s1, self.model.samples[s1].site_phase))
			print()
			raise Exception("Modeling aborted: Inconsistent stratigraphic relations detected.")
		
		# Update earlier-than relationships based on excavation area phases
		for i, s1 in enumerate(samples):
			if self.model.samples[s1].excavation_area_phase is None:
				continue
			for j, s2 in enumerate(samples):
				if s2 == s1:
					continue
				if self.model.samples[s1].area != self.model.samples[s2].area:
					continue
				earlier_than[i][j] = _is_earlier_eap(self.model.samples[s1].excavation_area_phase, self.model.samples[s2].excavation_area_phase)
		
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
		
		# Update earlier-than relationships based on site phases
		for i, s1 in enumerate(samples):
			if self.model.samples[s1].site_phase is None:
				continue
			for j, s2 in enumerate(samples):
				if s2 == s1:
					continue
				earlier_than[i][j] = _is_earlier_site_phase(self.model.samples[s1].site_phase, self.model.samples[s2].site_phase)
		
		# Check if earlier_than has circular relationships
		check_circular_relationships(earlier_than, samples)
		
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
		
		# Assign phase numbers to clusters
		phases_clu = dict([(label, n+1) for n, label in enumerate(labels)])
		
		collect = {}
		for label in phases_clu:
			for sample in clusters[label]:
				collect[sample] = phases_clu[label]
		phases_clu = collect
		# phases_clu = {sample: phase, ...}
		
		# Eliminate cluster pairs which violate chronological relations given by stratigraphy
		invalid = set()
		for i, s1 in enumerate(samples):
			if phases_clu[s1] is None:
				continue
			for j, s2 in enumerate(samples):
				if phases_clu[s2] is None:
					continue
				if phases_clu[s1] < phases_clu[s2]:
					if ((self.model.samples[s1].group == self.model.samples[s2].group) and (self.model.samples[s1].phase is not None) and (
							self.model.samples[s2].phase is not None) and (self.model.samples[s1].phase > self.model.samples[s2].phase)):
						invalid.add((phases_clu[s1], phases_clu[s2]))
						invalid.add((phases_clu[s2], phases_clu[s1]))
		
		# Update earlier-than relationships based on phases derived from clustering
		for i, s1 in enumerate(samples):
			if phases_clu[s1] is None:
				continue
			for j, s2 in enumerate(samples):
				if phases_clu[s2] is None:
					continue
				if (phases_clu[s1], phases_clu[s2]) in invalid:
					continue
				if phases_clu[s1] < phases_clu[s2]:
					earlier_than[i][j] = True
		
		# Check if earlier_than has circular relationships
		check_circular_relationships(earlier_than, samples)
		
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
					if (samples[i] not in ranges) or (samples[j] not in ranges): # already assigned as outliers
						continue
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
	
	def process(self, by_clusters: bool = False, by_dates: bool = False) -> bool:
		"""
		Updates the phasing of samples based on stratigraphic relations.

		Updates the groups and phases of samples based on their stratigraphic relations.

		:param by_clusters: If True, update the phasing by clustering sample dates. Defaults to False.
		:type by_clusters: bool, optional
		:param by_dates: If True, update the phasing by comparing sample dates. Defaults to False.
		:type by_dates: bool, optional		
		:return: True if phasing has changed, False otherwise.
		:rtype: bool
		"""
		
		phasing0 = set([(name, self.model.samples[name].group, self.model.samples[name].phasing_range[0], self.model.samples[name].phasing_range[1]) for name in self.model.samples])
		
		earlier_than, samples = self.create_earlier_than_matrix()
		if by_clusters and self.model.is_clustered:
			earlier_than = self.update_earlier_than_by_clustering(earlier_than, samples)
		
		if by_dates and self.model.is_modeled:
			earlier_than = self.update_earlier_than_by_dating(earlier_than, samples)
		
		ranges = [self.model.samples[name].get_range() for name in self.model.samples]		
		groups_phases, phasing_multi, phasing_ranges = get_groups_and_phases(earlier_than, samples, ranges)
		# groups_phases = {sample: [group, phase], ...}
		# phasing_multi = {sample: [phase, ...], ...}
		# phasing_ranges = {sample: [group, phase_min, phase_max], ...}
		
		print()  # DEBUG start
		print("Samples with phasing ranges:")
		for sample_id in phasing_ranges:
			gr, ph_min, ph_max = phasing_ranges[sample_id]
			if ph_min != ph_max:
				print("%s gr:%d ph:%d-%d" % (sample_id, gr, ph_min, ph_max))
		print()  # DEBUG end
		
		if phasing_multi:
			print("Warning! Some samples can be assigned to multiple phases:")
			for name in phasing_multi:
				print("%s: Phases %s" % (name, ", ".join([str(ph) for ph in phasing_multi[name]])))
		
		for name in self.model.samples:
			if name in groups_phases:
				group, phase = groups_phases[name]
				phasing_rng = phasing_ranges[name]
				self.model.samples[name].set_group(group)
				self.model.samples[name].set_phase(phase)
				self.model.samples[name].set_phasing_range(phasing_rng[1:])
			else:
				self.model.samples[name].set_group(None)
				self.model.samples[name].set_phase(None)
				self.model.samples[name].set_phasing_range(None)
		
		phasing1 = set([(name, self.model.samples[name].group, self.model.samples[name].phasing_range[0], self.model.samples[name].phasing_range[1]) for name in self.model.samples])
		
		if phasing0 != phasing1:
			return True
		return False
