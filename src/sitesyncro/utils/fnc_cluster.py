from itertools import combinations
from typing import List, Dict

import numpy as np
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform
from scipy.stats import wasserstein_distance
from sklearn.decomposition import PCA
from sklearn.metrics import silhouette_score

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

def calc_distance_matrix_new(distributions: List[np.ndarray]) -> np.ndarray:
	"""
	Calculate a distance matrix of distributions of calibrated C-14 dates based on their Wasserstein distance
	(Earth Mover's Distance).
	
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
	D = np.zeros((dists_n * (dists_n - 1)) // 2, dtype=float)
	
	time_axis = np.arange(len(distributions[0]))
	
	k = 0
	for d1, d2 in combinations(range(dists_n), 2):
		D[k] = wasserstein_distance(time_axis, time_axis, distributions[d1], distributions[d2])
		k += 1
	
	# Normalize distances between 0 and 1
	D = (D - D.min()) / (D.max() - D.min()) if D.max() != D.min() else np.zeros_like(D)
	
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
