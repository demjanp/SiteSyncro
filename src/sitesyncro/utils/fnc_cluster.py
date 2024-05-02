from itertools import combinations
from typing import List, Dict, Any

import numpy as np
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform
from scipy.stats import norm
from sklearn.decomposition import PCA
from sklearn.metrics import silhouette_score
from tqdm import tqdm

from sitesyncro.utils.fnc_mp import (process_mp)
from sitesyncro.utils.fnc_simulate import (calculate_parameters, generate_random_distributions)
from sitesyncro.utils.fnc_stat import (calc_sum, samples_to_distributions)


def calc_distance_matrix(distributions: List[np.ndarray]) -> np.ndarray:
	"""
	Calculate a distance matrix of distributions of calibrated C-14 dates based on probabilities
	that they represent the same event.

	Parameters:
	distributions (List[np.ndarray]): A list of distributions. Each distribution is represented
	as a numpy array of probabilities.

	Returns:
	np.ndarray: A 1-d condensed distance matrix D. The value at index [i,j] in the squareform of D
	represents the inverse probability that distributions i and j represent the same event.

	Example:
	If we have three distributions, the function will calculate the pairwise distances between
	these distributions and return a 1-d condensed distance matrix. If we convert this matrix
	back to a 2-d form using scipy's squareform function, we'll get a matrix where the value at
	index [i,j] represents the inverse probability that distributions i and j represent the same event.
	"""
	dists_n = len(distributions)
	PS = np.zeros((dists_n * (dists_n - 1)) // 2, dtype=float)
	k = 0
	for d1, d2 in combinations(range(dists_n), 2):
		PS[k] = np.dot(distributions[d1], distributions[d2])
		k += 1
	D = 1 - PS
	D = D - D.min()
	
	return D


def calc_distances_pca(D: np.ndarray, n_components: int = None) -> np.ndarray:
	"""
	Calculate the Principal Component Analysis (PCA) scores for each distribution based on their distances.

	This function takes a distance matrix D and reduces its dimensionality using PCA. The number of components
	is chosen so that they explain 99% of variance. If the number of components is not provided, it will be
	determined automatically.

	Parameters:
	D (np.ndarray): A distance matrix. If D is a 1-d condensed distance matrix, it will be converted to a 2-d form.
	n_components (int, optional): The number of principal components to keep. If not provided, it will be determined
	automatically to explain 99% of variance.

	Returns:
	np.ndarray: The PCA scores for each distribution. Each row corresponds to a distribution and each column corresponds
	to a PCA component.

	Raises:
	Warning: If PCA cannot be calculated on the distance matrix, a warning message is printed and the original distance
	matrix is returned.
	"""
	if D.ndim == 1:
		D = squareform(D)
	if n_components is None:
		pca = PCA(n_components=None)
		try:
			pca.fit(D)
		except:
			print("\nWarning: could not calculate PCA of the distance matrix!\n")
			return D
		
		n_components = np.where(np.cumsum(pca.explained_variance_ratio_) > 0.99)[0]
		if not n_components.size:
			n_components = 1
		else:
			n_components = n_components.min() + 1
	pca = PCA(n_components=n_components)
	pca.fit(D)
	
	return pca.transform(D)


def calc_clusters_hca(D: np.ndarray, n: int, method: str = "average") -> Dict[int, List[int]]:
	"""
	Cluster distributions into n clusters based on the distance matrix D using Hierarchical Cluster Analysis (HCA).

	This function takes a distance matrix D and the desired number of clusters n, and applies HCA to cluster the
	distributions. The clustering method can be specified, with "average" being the default method.

	Parameters:
	D (np.ndarray): A distance matrix. If D is a 1-d condensed distance matrix, it will be converted to a 2-d form.
	n (int): The desired number of clusters.
	method (str, optional): The method to use for HCA. Default is "average".

	Returns:
		clusters = {label: [idx, ...], ...}; idx = index D
	
	Example:
	If we have a distance matrix for three distributions and we want to cluster them into two groups, the function
	might return {1: [0, 2], 2: [1]}. This means that distributions 0 and 2 are in one cluster and distribution 1 is
	in another cluster.
	"""
	if D.ndim == 1:
		D = squareform(D)
	d = calc_distances_pca(D)
	if d.shape[0] == d.shape[1]:
		d = squareform(d)
	clusters_l = fcluster(linkage(d, method=method, metric="euclidean"), n, criterion="maxclust")
	labels = np.unique(clusters_l)
	clusters = {}
	for label in labels:
		clusters[int(label)] = np.where(clusters_l == label)[0].tolist()
	
	return clusters


def calc_silhouette(D: np.ndarray, clusters: Dict[int, List[int]]) -> float:
	"""
	Calculate the mean Silhouette Coefficient of all samples.

	The Silhouette Coefficient is calculated using the mean intra-cluster distance (a) and the mean nearest-cluster
	distance (b) for each sample. The Silhouette Coefficient for a sample is (b - a) / max(a, b). To clarify, b is the
	distance between a sample and the nearest cluster that the sample is not a part of.

	Parameters:
	D (np.ndarray): A distance matrix. If D is a 1-d condensed distance matrix, it will be converted to a 2-d form.
	clusters: {label: [idx, ...], ...}; idx = index in the distance matrix D
	
	Returns:
	float: The mean Silhouette Coefficient over all samples. The best value is 1 and the worst value is -1. Values near
	0 indicate overlapping clusters. Negative values generally indicate that a sample has been assigned to the wrong
	cluster, as a different cluster is more similar.

	Raises:
	ValueError: If the number of clusters is less than 2, a ValueError is raised.
	"""
	if len(clusters) < 2:
		return -1
	if D.ndim == 1:
		D = squareform(D)
	clusters_l = np.zeros(D.shape[0], dtype=int)
	for li, label in enumerate(list(clusters.keys())):
		clusters_l[clusters[label]] = li + 1
	
	return silhouette_score(D, clusters_l, metric="precomputed")


def cluster_distributions(model: object) -> (Dict[int, Dict[int, List[str]]], Dict[int, Dict[int, float]], float):
	"""
	Cluster distributions into n clusters using Hierarchical Cluster Analysis.

	This function takes a Model object as input, which contains the necessary data and parameters for clustering.
	The function returns a tuple containing three elements: a dictionary of clusters, a dictionary of means, and a float
	representing the silhouette score.

	Parameters:
	model (Model): An object of class Model
	
	Returns:
	(clusters, means, silhouette)
		- clusters = {label: [sample, ...], ...}
		- means = {label: mean, ...}
		- silhouette = mean value for all clusters
	
	Raises:
	Exception: If the number of clusters is less than 2, an exception is raised.
	"""
	if model.cluster_n < 2:
		raise Exception("Number of clusters must be >1")
	
	samples = []
	distributions = []
	for name in model.samples:
		if model.samples[name].is_modeled:
			samples.append(name)
			distributions.append(model.samples[name].posterior)
	# distributions = [[p, ...], ...]; ordered by samples
	
	if not samples:
		raise Exception("No modeled samples found")
	
	D = calc_distance_matrix(distributions)
	clusters = calc_clusters_hca(D, model.cluster_n)
	means = {}
	for label in clusters:
		means[label] = np.average(
			model.years,
			weights=calc_sum([distributions[idx] for idx in clusters[label]])
		)
	sil = calc_silhouette(D, clusters)
	clusters = dict([(label, [samples[idx] for idx in clusters[label]]) for label in clusters])
	
	return clusters, means, sil


def worker_fnc(params: Any, dates_n: int, t_param1: float, t_param2: float, uncertainties: List[float],
               uncertainty_base: float, curve: np.ndarray, uniform: bool) -> np.ndarray:
	distributions = generate_random_distributions(dates_n, t_param1, t_param2, uncertainties, uncertainty_base, curve,
	                                              uniform)
	D = calc_distance_matrix(distributions)
	return D


def collect_fnc(data: Any, D_pool: List[np.ndarray], pbar: tqdm) -> None:
	# data = D
	D_pool.append(data)
	pbar.n = len(D_pool)
	pbar.refresh()


def test_distribution_clustering(model: object, max_cpus: int = -1, max_queue_size: int = -1) -> (
		Dict[int, Dict[int, List[str]]], Dict[int, Dict[int, float]], Dict[int, float], Dict[int, float]):
	"""
	Test the clustering of distributions for randomness.

	This function takes a Model object as input, which contains the necessary data and parameters for clustering.
	The function returns a tuple containing four elements: a dictionary of clusters, a dictionary of means, a dictionary of silhouette scores, and a dictionary of p-values.

	Parameters:
	model (Model): An object of class Model
	max_cpus (int, optional): The maximum number of CPUs to use for multiprocessing. Default is -1, which means all available CPUs will be used.
	max_queue_size (int, optional): The maximum size of the queue for multiprocessing. Default is -1, which means the queue size is unlimited.

	Returns:
	(clusters, means, sils, ps)
		- clusters = {n: {label: [sample, ...], ...}, ...}; n = number of clusters
		- means = {n: {label: mean, ...}, ...}
		- sils = {n: Silhouette, ...}
		- ps = {n: p-value, ...}; p-value of the null hypothesis that the Silhouette for n clusters is the product of randomly distributed dates

	Raises:
	Exception: If the number of samples is less than 3, an exception is raised.
	"""
	
	distributions, samples, joined = samples_to_distributions(model.samples.values())
	# distributions = [[p, ...], ...]
	# samples = [sample name, ...] ordered by distributions
	# joined = {combined_name: [sample name, ...], ...}
	
	dates_n = len(samples)
	
	if dates_n < 3:
		raise Exception("Insufficient number samples for testing of clustering.")
	
	print("Testing clustering of %d distributions" % dates_n)
	
	sum_obs = calc_sum(distributions)
	t_param1, t_param2 = calculate_parameters(model.years, sum_obs, model.uniform)
	
	# Get dating range of all samples
	rng_min, rng_max = np.inf, -np.inf
	for name in model.samples:
		rng = model.samples[name].get_range()
		if None in rng:
			continue
		rng_min = min(rng_min, min(rng))
		rng_max = max(rng_max, max(rng))
	rng = (rng_max - rng_min)
	if abs(rng) == np.inf:
		raise Exception("Invalid dating range")
	clu_max = min(max(2, int(round(rng / model.min_years_per_cluster))), dates_n - 1)
	
	clusters = {}  # {cluster_n: {label: [idx, ...], ...}, ...}; idx = index in samples
	sils = {}  # {cluster_n: silhouette_score, ...}
	means = {}  # {cluster_n: {label: mean, ...}, ...}
	D = calc_distance_matrix(distributions)
	for cluster_n in range(2, clu_max + 1):
		clusters[cluster_n] = calc_clusters_hca(D, cluster_n)
		sils[cluster_n] = calc_silhouette(D, clusters[cluster_n])
		means[cluster_n] = {}
		for label in clusters[cluster_n]:
			means[cluster_n][label] = np.average(
				model.years,
				weights=calc_sum([distributions[idx] for idx in clusters[cluster_n][label]])
			)
		# convert clusters to {cluster_n: {label: [sample name, ...], ...}, ...}
		clusters[cluster_n] = dict(
			[(label, [samples[idx] for idx in clusters[cluster_n][label]]) for label in clusters[cluster_n]])
	
	ps = {}
	D_pool = []
	for cluster_n in clusters:
		sils_rnd = []
		sils_prev = None
		c = 0
		todo = model.npass
		with tqdm(total=todo * 2) as pbar:
			pbar.set_description("Clusters: %d/%d, Conv.: %0.3f" % (cluster_n, clu_max, c))
			i = -1
			while True:
				i += 1
				while i >= len(D_pool):
					process_mp(worker_fnc, range(max(4, (todo - len(D_pool)) + 1)),
					           [dates_n, t_param1, t_param2, model.uncertainties, model.uncertainty_base, model.curve,
					            model.uniform],
					           collect_fnc=collect_fnc, collect_args=[D_pool, pbar],
					           max_cpus=max_cpus, max_queue_size=max_queue_size)
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
					if c >= model.convergence:
						break
					todo *= 2
					pbar.total = max(todo, model.npass * 2)
					pbar.refresh()
		sils_rnd = np.array(sils_rnd)
		s = sils_rnd.std()
		if s == 0:
			p = 0
		else:
			p = 1 - float(norm(sils_rnd.mean(), s).cdf(sils[cluster_n]))
		ps[cluster_n] = p
	
	# Split joined samples
	if joined:
		for cluster_n in clusters:
			for label in clusters[cluster_n]:
				for name in joined:
					if name in clusters[cluster_n][label]:
						clusters[cluster_n][label].remove(name)
						clusters[cluster_n][label] += joined[name]
	
	return clusters, means, sils, ps


def find_opt_clusters(clusters: Dict[int, Dict[int, List[int]]], ps: Dict[int, float], sils: Dict[int, float],
                      p_value: float = 0.05) -> int or None:
	"""
	Find the optimal number of clusters based on Silhouette scores and p-values of clustering solutions.

	This function takes dictionaries of clusters, p-values, and Silhouette scores, and a p-value threshold as input.
	It returns the optimal number of clusters that have a p-value less than the threshold and the highest Silhouette score.

	Parameters:
	clusters: {n: {label: [sample, ...], ...}, ...}; n = number of clusters
	ps: {n: p-value, ...}; p-value of the null hypothesis that the Silhouette for n clusters is the product of randomly distributed dates
	sils: {n: Silhouette score, ...}
	p_value (float, optional): The p-value threshold for determining the optimal number of clusters. Default is 0.05.

	Returns:
	int: The optimal number of clusters. If no clusters have a p-value less than the threshold, the function returns None.

	Example:
	If we have three clustering solutions with the following p-values and Silhouette scores:
	- 2 clusters: p-value = 0.04, Silhouette score = 0.5
	- 3 clusters: p-value = 0.03, Silhouette score = 0.6
	- 4 clusters: p-value = 0.06, Silhouette score = 0.7
	And the p-value threshold is 0.05, the function will return 3 as the optimal number of clusters because the p-value for 3 clusters is less than the threshold and its Silhouette score is the highest among the solutions with a p-value less than the threshold.
	"""
	
	clu_ns = np.array(sorted(list(clusters.keys())), dtype=int)
	ps = np.array([ps[clu_n] for clu_n in clu_ns])
	sils = np.array([sils[clu_n] for clu_n in clu_ns])
	idxs = np.where(ps < p_value)[0]
	if not idxs.size:
		return None
	
	return int(clu_ns[idxs[np.argmax(sils[idxs])]])


def proc_clustering(model: object, max_cpus: int = -1, max_queue_size: int = -1) -> (
		Dict[int, Dict[int, List[str]]], Dict[int, float], Dict[int, float], int):
	"""
	Process the clustering of distributions.

	This function takes a Model object as input, which contains the necessary data and parameters for clustering.
	The function returns a tuple containing four elements: a dictionary of clusters, a dictionary of means, a dictionary of silhouette scores, and an integer representing the optimal number of clusters.

	Parameters:
	model (object): An object of class Model.
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
	
	if model.cluster_n > -1:
		clusters, means, sil = cluster_distributions(model)
		clusters = {model.cluster_n: clusters}
		means = {model.cluster_n: means}
		sils = {model.cluster_n: sil}
		ps = {model.cluster_n: 1}
		opt_n = model.cluster_n
	else:
		clusters, means, sils, ps = test_distribution_clustering(model, max_cpus=max_cpus,
		                                                         max_queue_size=max_queue_size)
		opt_n = find_opt_clusters(clusters, ps, sils, model.p_value)
	if opt_n is None:
		opt_n = 0
		clusters[0] = {0: sorted(list(model.samples.keys()))}
		means[0] = {0: 0}
		sils[0] = 0
		ps[0] = 1
	
	return clusters, means, sils, ps, opt_n
