from sitesyncro.lib.fnc_sum import (sum_distributions, calc_mean_std)
from sitesyncro.lib.fnc_simulate import (calculate_parameters, generate_random_distributions)
from sitesyncro.lib.fnc_mp import (process_mp)

import numpy as np
from scipy.stats import norm
from itertools import combinations
from sklearn.decomposition import PCA
from sklearn.metrics import silhouette_score
from scipy.spatial.distance import squareform
from scipy.cluster.hierarchy import linkage, fcluster

def p_same_event(dist1, dist2):
	# Calculate probability that distributions dist1 and dist2 represent the same event
	#
	# Returns probability value
	
	return (4*(dist1 * dist2).sum()) / (dist1.sum() + dist2.sum())**2

def calc_distance_matrix(distributions):
	# Calculate a distance matrix of distributions of calibrated C-14 dates based on probabilities that they represent the same event
	# distributions = [[p, ...], ...]
	#
	# Returns D[i,j] = d; i,j = index in distributions; d = inverse probability that distributions i and j represent the same event
	
	dists_n = len(distributions)
	PS = np.zeros((dists_n, dists_n), dtype = float)
	for d1, d2 in combinations(range(dists_n), 2):
		p = p_same_event(distributions[d1], distributions[d2])
		PS[d1,d2] = p
		PS[d2,d1] = p
	D = 1 - PS
	D = D - D.min()
	for i in range(D.shape[0]):
		D[i, i] = 0
	
	return D

def calc_distances_pca(D, n_components = None):
	# Principal Component Analysis scores for each distribution based on their distances as defined in calc_distance_matrix
	# Number of components is chosen so that they explain 99% of variance
	# D[i,j] = d; i,j = index in distributions; d = inverse probability that distributions i and j represent the same event
	#
	# Returns S[i,k] = PCA score; i = index in distributions; k = index of PCA component
	
	if n_components is None:
		pca = PCA(n_components = None)
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
	pca = PCA(n_components = n_components)
	pca.fit(D)
	
	return pca.transform(D)

def calc_clusters_hca(D, n, method = "average"):
	# Cluster distributions into n clusters based on the distance matrix D using Hierarchical Cluster Analysis
	#
	# Returns clusters = {label: [idx, ...], ...}; idx = index D
	
	d = calc_distances_pca(D)
	if d.shape[0] == d.shape[1]:
		d = squareform(d)
	clusters_l = fcluster(linkage(d, method = method, metric = "euclidean"), n, criterion = "maxclust")
	labels = np.unique(clusters_l)
	clusters = {}
	for label in labels:
		clusters[int(label)] = np.where(clusters_l == label)[0].tolist()
	
	return clusters

def calc_silhouette(D, clusters):
	# Calculate mean Silhouette of clustered distributions as defined by 
	#	Rousseeuw (1987, https://doi.org/10.1016/0377-0427(87)90125-7)
	# clusters = {label: [idx, ...], ...}; idx = index in the distance matrix D
	#
	# Returns silhouette_score
	
	if len(clusters) < 2:
		return -1
	
	clusters_l = np.zeros(D.shape[0], dtype = int)
	for li, label in enumerate(list(clusters.keys())):
		clusters_l[clusters[label]] = li + 1
	
	return silhouette_score(D, clusters_l, metric = "precomputed")

def cluster_distributions(distributions, years, clusters_n):
	# Cluster distributions into n clusters using Hierarchical Cluster Analysis
	# distributions = {sample: [p, ...], ...}
	#
	# Returns clusters, means, silhouette
	# 	clusters = {label: [sample, ...], ...}
	#   means = {label: mean, ...}
	
	samples = sorted(list(distributions.keys()))
	distributions = [distributions[sample] for sample in samples]
	
	D = calc_distance_matrix(distributions)
	clusters = calc_clusters_hca(D, clusters_n)
	means = {}
	for label in clusters:
		means[label] = np.average(
			years,
			weights=sum_distributions([distributions[idx] for idx in clusters[label]])
		)
	sil = calc_silhouette(D, clusters)
	clusters = dict([(label, [samples[idx] for idx in clusters[label]]) for label in clusters])
	
	return clusters, means, sil

def worker_fnc(params, clusters_n, dates_n, t_param1, t_param2, uncert_base, curve, uniform):
	
	distributions = generate_random_distributions(dates_n, t_param1, t_param2, uncert_base, curve, uniform)
	D = calc_distance_matrix(distributions)
	sil = calc_silhouette(D, calc_clusters_hca(D, clusters_n))
	
	return sil

def collect_fnc(data, results):
	
	# data = sil
	results.append(data)

def progress_fnc(done, todo, clusters_n, clu_max, all_done, all_todo, c):
	
	print("\rClusters: %d/%d, Iteration: %d/%d, Conv: %0.4f         " % (clusters_n, clu_max, all_done + done, all_todo, c), end = "")

def test_distribution_clustering(distributions, curve, 
		uncert_base = 15, uniform = False,
		npass = 100, convergence = 0.99, max_cpus = -1, max_queue_size = 100):
	# Test the clustering of distributions for randomness
	#
	# distributions = {sample: [p, ...], ...}
	# curve: [[CalBP, ConvBP, CalSigma], ...]
	# uncert_base: base uncertainty to simulate the radiocarbon dates
	# uniform: flag indicating whether to use a uniform distribution for the calendar ages
	#
	# returns clusters, means, sils, ps
	# 	clusters = {n: {label: [sample, ...], ...}, ...}; n = number of clusters
	#   means = {n: {label: mean, ...}, ...}
	#	sils = {n: Silhouette, ...}
	#	ps = {n: p-value, ...}; p-value of the null hypothesis that the 
	#		Silhouette for n clusters is the product of randomly distributed dates
	
	samples = sorted(list(distributions.keys()))
	distributions = [distributions[sample] for sample in samples]
	dates_n = len(distributions)
	
	years = curve[:, 0]
	sum_obs = sum_distributions(distributions)
	t_param1, t_param2 = calculate_parameters(years, sum_obs, uniform)
	
	clusters = {}  # {clusters_n: {label: [idx, ...], ...}, ...}; idx = index in samples
	sils = {} # {clusters_n: silhouette_score, ...}
	means = {} # {clusters_n: {label: mean, ...}, ...}
	D = calc_distance_matrix(distributions)
	for clusters_n in range(2, dates_n - 1):
		clusters[clusters_n] = calc_clusters_hca(D, clusters_n)
		sils[clusters_n] = calc_silhouette(D, clusters[clusters_n])
		means[clusters_n] = {}
		for label in clusters[clusters_n]:
			means[clusters_n][label] = np.average(
				years,
				weights=sum_distributions([distributions[idx] for idx in clusters[clusters_n][label]])
			)
		clusters[clusters_n] = dict([(label, [samples[idx] for idx in clusters[clusters_n][label]]) for label in clusters[clusters_n]])
	
	clu_max = dates_n - 1
	
	# Calculate the mean or median and standard deviation or range of the summed distribution
	ps = {}
	for clusters_n in clusters:
		sils_rnd = []
		sils_prev = None
		c = 0
		params_list = list(range(npass))
		todo = npass
		while True:
			process_mp(worker_fnc, params_list, [clusters_n, dates_n, t_param1, t_param2, uncert_base, curve, uniform],
				collect_fnc = collect_fnc, collect_args = [sils_rnd],
				progress_fnc = progress_fnc, progress_args = [clusters_n, clu_max, len(sils_rnd), npass*2, c],
				max_cpus = max_cpus, max_queue_size = max_queue_size)
			if len(sils_rnd) >= todo:
				sils_m = np.array(sils_rnd).mean()
				if sils_prev is not None:
					c = 1 - np.abs(sils_prev - sils_m) / sils_prev
				sils_prev = sils_m
				if c >= convergence:
					print("\nConverged at:", c)
					break
				todo *= 2
		sils_rnd = np.array(sils_rnd)
		s = sils_rnd.std()
		if s == 0:
			p = 0
		else:
			p = 1 - float(norm(sils_rnd.mean(), s).cdf(sils[clusters_n]))
		ps[clusters_n] = p				
	
	return clusters, means, sils, ps

def find_opt_clusters(clusters, ps, sils, p_value = 0.05):
	# Find optimal number of clusters based on Silhouette scores and p-values of clustering solutions
	#
	# clusters = {n: {label: [sample, ...], ...}, ...}; n = number of clusters
	# ps = {n: p-value, ...}; p-value of the null hypothesis that the Silhouette for n clusters is the product of randomly distributed dates
	# sils = {n: Silhouette score, ...}
	#
	# Returns number of clusters
	
	clu_ns = np.array(sorted([n for n in clusters]), dtype = int)
	ps = np.array([ps[clu_n] for clu_n in clu_ns])
	sils = np.array([sils[clu_n] for clu_n in clu_ns])
	idxs = np.where(ps < p_value)[0]
	if not idxs.size:
		return None
	
	return int(clu_ns[idxs[np.argmax(sils[idxs])]])

