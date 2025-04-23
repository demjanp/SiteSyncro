import os
import json
import numpy as np
import matplotlib.pyplot as plt

DIRECTORY = "c:\\documents_synced\\SiteSyncro\\sitesyncro\\silhouette"
FORIG = os.path.join(DIRECTORY, "test_clu_orig\\test_results_orig.json")
FNEW = os.path.join(DIRECTORY, "test_clu_new\\test_results_new.json")

def calc_clustering_index(clusters_n_expected, clusters_n_observed):
	if clusters_n_observed <= 0:
		return 0
	if clusters_n_observed == clusters_n_expected:
		return 1.0
	return max(0, 1 - abs(clusters_n_observed - clusters_n_expected) / clusters_n_expected)

if __name__ == '__main__':
	with open(FORIG, 'r') as file:
		data_orig = json.load(file)
	with open(FNEW, 'r') as file:
		data_new = json.load(file)

	clustering_idx_orig = {}
	clustering_idx_new = {}
	phasing_idx_orig = {}
	phasing_idx_new = {}

	for clusters_n in data_orig:
		clustering_idx_orig[clusters_n] = {}
		clustering_idx_new[clusters_n] = {}
		phasing_idx_orig[clusters_n] = {}
		phasing_idx_new[clusters_n] = {}	
		for sample_size in data_orig[clusters_n]:
			_, samples_n, _ = data_orig[clusters_n][sample_size][0]
			clustering_idx_orig[clusters_n][samples_n] = []
			clustering_idx_new[clusters_n][samples_n] = []
			phasing_idx_orig[clusters_n][samples_n] = []
			phasing_idx_new[clusters_n][samples_n] = []
			for ph_good, samples_n, cluster_opt_n in data_orig[clusters_n][sample_size]:
				clustering_idx_orig[clusters_n][samples_n].append(
					calc_clustering_index(int(clusters_n), cluster_opt_n)
				)
				phasing_idx_orig[clusters_n][samples_n].append(ph_good / samples_n)
			for ph_good, samples_n, cluster_opt_n in data_new[clusters_n][sample_size]:
				clustering_idx_new[clusters_n][samples_n].append(
					calc_clustering_index(int(clusters_n), cluster_opt_n)
				)
				phasing_idx_new[clusters_n][samples_n].append(ph_good / samples_n)

	# Prepare the figure canvas dynamically based on the number of unique combinations
	plot_keys = []
	for clusters_n in sorted(clustering_idx_orig.keys(), key=int):  # clusters_n is still numeric as a string
		for samples_n in sorted(clustering_idx_orig[clusters_n].keys(), key=int):  # no key=int here
			plot_keys.append((clusters_n, samples_n))
	
	nrows = len(plot_keys)
	fig, axs = plt.subplots(nrows, 2, figsize=(12, 4 * nrows))
	if nrows == 1:
		axs = np.expand_dims(axs, axis=0)  # Ensure axs is always 2D

	for idx, (clusters_n, samples_n) in enumerate(plot_keys):
		title = f"{clusters_n} clusters, {samples_n} samples"
		clu_data = [
			clustering_idx_orig[clusters_n][samples_n],
			clustering_idx_new[clusters_n][samples_n]
		]
		pha_data = [
			phasing_idx_orig[clusters_n][samples_n],
			phasing_idx_new[clusters_n][samples_n]
		]

		axs[idx, 0].boxplot(clu_data, tick_labels=["Original", "New"])
		axs[idx, 0].set_title(f"Clustering - {title}")
		axs[idx, 0].set_ylim(-0.05, 1.05)

		axs[idx, 1].boxplot(pha_data, tick_labels=["Original", "New"])
		axs[idx, 1].set_title(f"Phasing - {title}")
		axs[idx, 1].set_ylim(-0.05, 1.05)

	plt.tight_layout()
	pdf_path = os.path.join(DIRECTORY, "comparison_boxplots.pdf")
	plt.savefig(pdf_path)
	print(f"Saved box plot comparison to: {pdf_path}")
