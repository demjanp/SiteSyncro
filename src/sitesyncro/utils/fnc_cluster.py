from sitesyncro.utils.fnc_mp import (process_mp)
from sitesyncro.utils.fnc_stat import (calc_sum, calc_mean_std, samples_to_distributions)
from sitesyncro.utils.fnc_simulate import (calculate_parameters, generate_random_distributions)

import numpy as np
from tqdm import tqdm
from scipy.stats import norm
from itertools import combinations
from collections import defaultdict
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

def cluster_distributions(model):
	# Cluster distributions into n clusters using Hierarchical Cluster Analysis
	#
	# Returns clusters, means, silhouette
	# 	clusters = {label: [sample, ...], ...}
	#   means = {label: mean, ...}
	#	silhouette = mean value for all clusters
	
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

def worker_fnc(params, cluster_n, dates_n, t_param1, t_param2, uncertainties, uncertainty_base, curve, uniform):
	
	distributions = generate_random_distributions(dates_n, t_param1, t_param2, uncertainties, uncertainty_base, curve, uniform)
	D = calc_distance_matrix(distributions)
	sil = calc_silhouette(D, calc_clusters_hca(D, cluster_n))
	
	return sil

def collect_fnc(data, results, pbar):
	
	# data = sil
	pbar.update(1)
	results.append(data)

def test_distribution_clustering(model, max_cpus = -1, max_queue_size = 10000):
	# Test the clustering of distributions for randomness
	#
	# returns clusters, means, sils, ps
	# 	clusters = {n: {label: [sample, ...], ...}, ...}; n = number of clusters
	#   means = {n: {label: mean, ...}, ...}
	#	sils = {n: Silhouette, ...}
	#	ps = {n: p-value, ...}; p-value of the null hypothesis that the 
	#		Silhouette for n clusters is the product of randomly distributed dates
	
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
	
	clusters = {}  # {cluster_n: {label: [idx, ...], ...}, ...}; idx = index in samples
	sils = {} # {cluster_n: silhouette_score, ...}
	means = {} # {cluster_n: {label: mean, ...}, ...}
	D = calc_distance_matrix(distributions)
	for cluster_n in range(2, dates_n):
		clusters[cluster_n] = calc_clusters_hca(D, cluster_n)
		sils[cluster_n] = calc_silhouette(D, clusters[cluster_n])
		means[cluster_n] = {}
		for label in clusters[cluster_n]:
			means[cluster_n][label] = np.average(
				model.years,
				weights=calc_sum([distributions[idx] for idx in clusters[cluster_n][label]])
			)
		# convert clusters to {cluster_n: {label: [sample name, ...], ...}, ...}
		clusters[cluster_n] = dict([(label, [samples[idx] for idx in clusters[cluster_n][label]]) for label in clusters[cluster_n]])
	
	clu_max = dates_n - 1	
	ps = {}
	for cluster_n in clusters:
		sils_rnd = []
		sils_prev = None
		c = 0
		todo = model.npass
		params_list = list(range(model.npass))
		with tqdm(total=todo) as pbar:
			pbar.total = model.npass*2
			pbar.set_description("Clusters: %d/%d, Conv.: %0.3f" % (cluster_n, clu_max, c))
			while True:
				process_mp(worker_fnc, params_list, [cluster_n, dates_n, t_param1, t_param2, model.uncertainties, model.uncertainty_base, model.curve, model.uniform],
		           collect_fnc = collect_fnc, collect_args = [sils_rnd, pbar],
		           max_cpus = max_cpus, max_queue_size = max_queue_size)
				if len(sils_rnd) >= todo:
					sils_m = np.array(sils_rnd).mean()
					if sils_prev is not None:
						c = 1 - np.abs(sils_prev - sils_m) / sils_prev
						pbar.set_description("Clusters: %d/%d, Conv.: %0.3f" % (cluster_n, clu_max, c))
					sils_prev = sils_m
					if c >= model.convergence:
						break
					todo *= 2
					pbar.total = max(todo, model.npass*2)
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

def find_opt_clusters(clusters, ps, sils, p_value = 0.05):
	# Find optimal number of clusters based on Silhouette scores and p-values of clustering solutions
	#
	# clusters = {n: {label: [sample, ...], ...}, ...}; n = number of clusters
	# ps = {n: p-value, ...}; p-value of the null hypothesis that the Silhouette for n clusters is the product of randomly distributed dates
	# sils = {n: Silhouette score, ...}
	#
	# Returns number of clusters
	
	clu_ns = np.array(sorted(list(clusters.keys())), dtype = int)
	ps = np.array([ps[clu_n] for clu_n in clu_ns])
	sils = np.array([sils[clu_n] for clu_n in clu_ns])
	idxs = np.where(ps < p_value)[0]
	if not idxs.size:
		return None
	
	return int(clu_ns[idxs[np.argmax(sils[idxs])]])

def proc_clustering(model, max_cpus = -1, max_queue_size = 10000):
	if model.cluster_n > -1:
		clusters, means, sil = cluster_distributions(model)
		clusters = {model.cluster_n: clusters}
		means = {model.cluster_n: means}
		sils = {model.cluster_n: sil}
		ps = {model.cluster_n: 1}
		opt_n = model.cluster_n
	else:
		clusters, means, sils, ps = test_distribution_clustering(model, max_cpus, max_queue_size)
		opt_n = find_opt_clusters(clusters, ps, sils, model.p_value)
	if opt_n is None:
		opt_n = 0
		clusters[0] = {0: sorted(list(model.samples.keys()))}
		means[0] = {0: 0}
		sils[0] = 0
		ps[0] = 1
		
	return clusters, means, sils, ps, opt_n

