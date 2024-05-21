from typing import List, Dict, Any

import numpy as np
from scipy.spatial.distance import squareform
from scipy.stats import norm
from tqdm import tqdm

from sitesyncro.utils.fnc_mp import (process_mp)
from sitesyncro.utils.fnc_simulate import (get_params, generate_random_distributions)
from sitesyncro.utils.fnc_stat import (calc_sum, samples_to_distributions)
from sitesyncro.utils.fnc_cluster import (calc_distance_matrix, calc_clusters_hca, calc_silhouette)

def worker_fnc(params: Any, dates_n: int, t_mean: float, t_std: float, uncertainties: List[float],
				uncertainty_base: float, curve: np.ndarray, uniform: bool) -> np.ndarray:
	
	distributions = generate_random_distributions(dates_n, t_mean, t_std, uncertainties, uncertainty_base, curve, uniform)
	D = calc_distance_matrix(distributions)
	return D


def collect_fnc(data: Any, D_pool: List[np.ndarray], pbar: tqdm) -> None:
	# data = D
	D_pool.append(data)
	pbar.n = len(D_pool)
	pbar.refresh()


class MCluster(object):
	"""
	A class implementing clustering functionality
	
	:param model: The parent Model object
	:type model: Model
	"""
	def __init__(self, model: 'Model'):
		
		self.model = model
	
	def get_clusterings(self) -> (Dict[int, Dict[int, List[str]]], Dict[int, Dict[int, float]], Dict[int, float], int, float, float):
		"""
		Cluster distributions using Hierarchical Cluster Analysis.
		
		Returns:
		(clusters, means, silhouette, t_mean, t_std)
			- clusters = {n: {label: [sample, ...], ...}, ...}; n = number of clusters
			- means = {n: {label: mean, ...}, ...}
			- sils = {n: Silhouette, ...}
			- distributions_n = number of distributions
			- t_mean = mean of the distributions
			- t_std = standard deviation of the distributions
		
		Raises:
		Exception: If the number of samples is less than 3, an exception is raised.
		"""
		
		distributions, samples, joined = samples_to_distributions(self.model.samples.values())
		# distributions = [[p, ...], ...]
		# samples = [sample name, ...] ordered by distributions
		# joined = {combined_name: [sample name, ...], ...}
		
		distributions_n = len(samples)
		
		if distributions_n < 3:
			raise Exception("Insufficient number samples for testing of clustering.")
		
		t_mean, t_std = get_params(distributions, self.model.curve, self.model.uniform)
		
		clusters = {}  # {cluster_n: {label: [idx, ...], ...}, ...}; idx = index in samples
		sils = {}  # {cluster_n: silhouette_score, ...}
		means = {}  # {cluster_n: {label: mean, ...}, ...}
		D = calc_distance_matrix(distributions)
		for cluster_n in range(2, distributions_n):
			clusters[cluster_n] = calc_clusters_hca(D, cluster_n)
			sils[cluster_n] = calc_silhouette(D, clusters[cluster_n])
			means[cluster_n] = {}
			for label in clusters[cluster_n]:
				means[cluster_n][label] = np.average(
					self.model.years,
					weights=calc_sum([distributions[idx] for idx in clusters[cluster_n][label]])
				)
			# convert clusters to {cluster_n: {label: [sample name, ...], ...}, ...}
			clusters[cluster_n] = dict(
				[(label, [samples[idx] for idx in clusters[cluster_n][label]]) for label in clusters[cluster_n]])
		
		# Split joined samples
		if joined:
			for cluster_n in clusters:
				for label in clusters[cluster_n]:
					for name in joined:
						if name in clusters[cluster_n][label]:
							clusters[cluster_n][label].remove(name)
							clusters[cluster_n][label] += joined[name]
		
		return clusters, means, sils, distributions_n, t_mean, t_std
	
	def find_opt_clusters_mcst(self, max_cpus: int = -1, max_queue_size: int = -1) -> (
			Dict[int, Dict[int, List[str]]], Dict[int, Dict[int, float]], Dict[int, float], Dict[int, float]):
		"""
		Test the clustering of distributions for randomness.

		This function takes a Model object as input, which contains the necessary data and parameters for clustering.
		The function returns a tuple containing four elements: a dictionary of clusters, a dictionary of means, a dictionary of silhouette scores, and a dictionary of p-values.

		Parameters:
		max_cpus (int, optional): The maximum number of CPUs to use for multiprocessing. Default is -1, which means all available CPUs will be used.
		max_queue_size (int, optional): The maximum size of the queue for multiprocessing. Default is -1, which means the queue size is unlimited.

		Returns:
		(clusters, means, sils, ps, opt_n)
			- clusters = {n: {label: [sample, ...], ...}, ...}; n = number of clusters
			- means = {n: {label: mean, ...}, ...}
			- sils = {n: Silhouette, ...}
			- ps = {n: p-value, ...}; p-value of the null hypothesis that the Silhouette for n clusters is the product of randomly distributed dates
			- opt_n = optimal number of clusters

		Raises:
		Exception: If the number of samples is less than 3, an exception is raised.
		"""
		
		def _find_opt(clusters, ps, sils, p_value):
			clu_ns = np.array(sorted(list(clusters.keys())), dtype=int)
			ps = np.array([ps[clu_n] for clu_n in clu_ns])
			sils = np.array([sils[clu_n] for clu_n in clu_ns])
			idxs = np.where(ps < p_value)[0]
			if not idxs.size:
				return None
			
			return int(clu_ns[idxs[np.argmax(sils[idxs])]])
		
		clusters, means, sils, distributions_n, t_mean, t_std = self.get_clusterings()
		
		clu_max = max(clusters.keys())
		
		ps = {}
		D_pool = []
		for cluster_n in clusters:
			sils_rnd = []
			sils_prev = None
			c = 0
			todo = self.model.npass
			with tqdm(total=todo * 2) as pbar:
				pbar.set_description("Clusters: %d/%d, Conv.: %0.3f" % (cluster_n, clu_max, c))
				i = -1
				while True:
					i += 1
					while i >= len(D_pool):
						n_dists = max(4, (todo - len(D_pool)) + 1)
						if n_dists > 50:						
							process_mp(worker_fnc, range(n_dists),
									   [distributions_n, t_mean, t_std, self.model.uncertainties, self.model.uncertainty_base, self.model.curve,
										self.model.uniform],
									   collect_fnc=collect_fnc, collect_args=[D_pool, pbar],
									   max_cpus=max_cpus, max_queue_size=max_queue_size)
						else:
							for j in range(n_dists):
								D_pool.append(worker_fnc(j, distributions_n, t_mean, t_std, self.model.uncertainties, self.model.uncertainty_base,
														 self.model.curve, self.model.uniform))
								pbar.n = len(D_pool)
								pbar.refresh()
					
					D = squareform(D_pool[i])
					sils_rnd.append(calc_silhouette(D, calc_clusters_hca(D, cluster_n)))
					pbar.n = len(sils_rnd)
					pbar.refresh()
					if len(sils_rnd) >= todo:
						sils_m = np.array(sils_rnd).mean()
						if sils_prev is not None:
							c = 1 - np.abs(sils_prev - sils_m) / sils_prev
							pbar.set_description("Clusters: %d/%d, Conv.: %0.3f" % (cluster_n, clu_max, c))
						sils_prev = sils_m
						if c >= self.model.convergence:
							break
						todo *= 2
						pbar.total = max(todo, self.model.npass * 2)
						pbar.refresh()
			sils_rnd = np.array(sils_rnd)
			s = sils_rnd.std()
			if s == 0:
				p = 0
			else:
				p = 1 - float(norm(sils_rnd.mean(), s).cdf(sils[cluster_n]))
			ps[cluster_n] = p
		
		opt_n = _find_opt(clusters, ps, sils, self.model.p_value)
		
		return clusters, means, sils, ps, opt_n
	
	def find_n_clusters(self) -> (Dict[int, Dict[int, List[str]]], Dict[int, Dict[int, float]], float):
		"""
		Cluster distributions into n clusters using Hierarchical Cluster Analysis.
		
		Returns:
		(clusters, means, silhouette)
			- clusters = {label: [sample, ...], ...}
			- means = {label: mean, ...}
			- silhouette = mean value for all clusters
		
		Raises:
		Exception: If the number of clusters is less than 2, an exception is raised.
		"""
		if self.model.cluster_n < 2:
			raise Exception("Number of clusters must be >1")
		
		clusters, means, sils, _, _, _ = self.get_clusterings()
		
		return clusters[self.model.cluster_n], means[self.model.cluster_n], sils[self.model.cluster_n]
	
	def find_opt_clusters_silhouette(self) -> (
			Dict[int, Dict[int, List[str]]], Dict[int, Dict[int, float]], int):
		"""
		Find the optimal number of clusters using the Silhouette method.
		
		Returns:
		(clusters, means, sils, opt_n)
			- clusters = {label: [sample, ...], ...}
			- means = {label: mean, ...}
			- sils = {n: Silhouette, ...}; n = number of clusters
			- opt_n = optimal number of clusters
		
		Raises:
		Exception: If the number of samples is less than 3, an exception is raised.
		
		"""
		
		clusters, means, sils, _, _, _ = self.get_clusterings()
		
		# Select cluster solution with highest Silhouette score
		opt_n = max(sils, key=sils.get)
		
		return clusters, means, sils, opt_n
	
	def process(self, max_cpus: int = -1, max_queue_size: int = -1) -> (
			Dict[int, Dict[int, List[str]]], Dict[int, float], Dict[int, float], int):
		"""
		Process the clustering of distributions.

		This function takes a Model object as input, which contains the necessary data and parameters for clustering.
		The function returns a tuple containing four elements: a dictionary of clusters, a dictionary of means, a dictionary of silhouette scores, and an integer representing the optimal number of clusters.

		Parameters:
		max_cpus (int, optional): The maximum number of CPUs to use for multiprocessing. Default is -1, which means all available CPUs will be used.
		max_queue_size (int, optional): The maximum size of the queue for multiprocessing. Default is -1, which means the queue size is unlimited.
		
		Returns:
		(clusters, means, sils, ps, opt_n)
			- clusters = {n: {label: [sample, ...], ...}, ...}; n = number of clusters
			- means = {n: {label: mean, ...}, ...}
			- sils = {n: Silhouette, ...}
			- ps = {n: p-value, ...}; p-value of the null hypothesis that the Silhouette for n clusters is the product of randomly distributed dates
			- opt_n = optimal number of clusters
		"""
		
		if self.model.cluster_n > -1:
			clusters, means, sil = self.find_n_clusters()
			clusters = {self.model.cluster_n: clusters}
			means = {self.model.cluster_n: means}
			sils = {self.model.cluster_n: sil}
			ps = {self.model.cluster_n: 1}
			opt_n = self.model.cluster_n
		
		elif self.model.cluster_selection == "silhouette":
			clusters, means, sils, opt_n = self.find_opt_clusters_silhouette()
			ps = dict([(n, 1) for n in clusters])
		
		elif self.model.cluster_selection == "mcst":
			clusters, means, sils, ps, opt_n = self.find_opt_clusters_mcst(max_cpus=max_cpus,
																	 max_queue_size=max_queue_size)
		
		else:
			raise Exception("Invalid cluster selection method: %s" % self.model.cluster_selection)
		
		if opt_n is None:
			opt_n = 0
			clusters[0] = {0: sorted(list(self.model.samples.keys()))}
			means[0] = {0: 0}
			sils[0] = 0
			ps[0] = 1
		
		return clusters, means, sils, ps, opt_n

