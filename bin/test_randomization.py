from sitesyncro import Model
from sitesyncro.utils.fnc_radiocarbon import (get_curve)

import os
import json
import numpy as np

UNCERTAINTY_BASE = 30
GAP = 300
CENTER_DATE = 3500
DATES_N = {'small': 10, 'large': 50}
ITERATIONS = 10

FRESULTS = "test_results.json"
DIRECTORY = "test_rnd"

def sim_event(cal_age_bp, curve):
	# Simulate a C-14 dated event
	
	# Find the closest index in the calibration curve
	idx = np.argmin(np.abs(curve[:, 0] - cal_age_bp))
	# Get the corresponding radiocarbon age
	return curve[idx, 1]

def sim_clustered(cal_bp_mean, gap, clusters_n, dates_n, uncertainty_base, curve):
	# Simulate dates_n dates from clusters_n of events with gaps of gap years around cal_bp_mean
	# returns [[phase, c14age, uncert], ...]
	
	if clusters_n == 0:
		# Generate uniformly distributed events around cal_bp_mean +-gap
		means = np.random.uniform(cal_bp_mean - gap, cal_bp_mean + gap, dates_n)
		means.sort()
		dates = []
		for cal_age in means:
			c14age = sim_event(cal_age, curve)
			uncert = uncertainty_base * np.exp(c14age / (2 * 8033))
			dates.append([1, c14age, np.random.normal(uncert, 2)])
		return dates
	
	if clusters_n == 1:
		# Generate normally distributed events around cal_bp_mean +-gap
		means = np.random.normal(cal_bp_mean, gap / 2, dates_n)
		means.sort()
		dates = []
		for cal_age in means:
			c14age = sim_event(cal_age, curve)
			uncert = uncertainty_base * np.exp(c14age / (2 * 8033))
			dates.append([1, c14age, np.random.normal(uncert, 2)])
		return dates
	
	# Randomly split dates_n between the clusters
	ns = np.random.multinomial(dates_n, [1 / clusters_n] * clusters_n)
	
	if clusters_n == 2:
		means = [cal_bp_mean - gap / 2, cal_bp_mean + gap / 2]
		means.sort()
	else:
		rng = 1.2 * gap * (clusters_n - 1)
		L = cal_bp_mean - rng / 2
		U = cal_bp_mean + rng / 2
		means = [L] + sorted(np.random.uniform(L, U, clusters_n - 2).tolist()) + [U]
		for i in range(1, clusters_n):
			if means[i] - means[i - 1] < 0.8 * gap:
				means[i] = means[i - 1] + 0.8 * gap
		means = np.array(means)
	means = means[::-1]
	dates = []
	for i, mean in enumerate(means):
		for cal_age_bp in np.random.normal(mean, gap / 2, ns[i]):
			c14age = sim_event(cal_age_bp, curve)
			uncert = uncertainty_base * np.exp(c14age / (2 * 8033))
			dates.append([i + 1, c14age, np.random.normal(uncert, 2)])
	return dates

if __name__ == '__main__':
	
	curve = get_curve()
	
	results = {}
	for clusters_n in [0, 1, 4]:
		results[clusters_n] = {}
		for sample_size in DATES_N:
			results[clusters_n][sample_size] = []
			for i in range(ITERATIONS):
				print("\nClusters: %d, Sample size: %s, Iter. %d/%d\n" % (clusters_n, sample_size, i + 1, ITERATIONS))
				dates = sim_clustered(CENTER_DATE, GAP, clusters_n, DATES_N[sample_size], UNCERTAINTY_BASE, curve)
				model = Model(directory = DIRECTORY, uniform = (clusters_n == 0), overwrite = True)
				n = 1
				for phase, age, uncert in dates:
					name = "%d_%d" % (phase, n)
					n += 1
					model.add_sample(name, age, uncert)
#				model.process_randomization()
#				model.plot_randomized(fname = os.path.join(DIRECTORY, "rnd_%d_%s_%03d.pdf" % (clusters_n, sample_size, i)))
				model.process_clustering()
				model.plot_clusters(fname = os.path.join(DIRECTORY, "clu_%d_%s_%03d.pdf" % (clusters_n, sample_size, i)))
				model.process_phasing(by_clusters = True)
#				distr_check = None
				distr_check = 1   # DEBUG
				clustering_check = None
				phasing_check = None
				
				if clusters_n == 0:
					# uniform distribution of events
#					distr_check = float(model.random_p >= model.p_value)
					# no clustering of events
					clustering_check = float(model.cluster_opt_n < 2)
					# no phasing
					phasing_check = 1 / len(set([model.samples[name].phase for name in model.samples]))
				
				elif clusters_n == 1:
					# normal distribution of events
#					distr_check = float(model.random_p >= model.p_value)
					# no clustering of events
					clustering_check = float(model.cluster_opt_n < 2)
					# no phasing
					phasing_check = 1 / len(set([model.samples[name].phase for name in model.samples]))
				
				else:
					# non-normal distribution of events
#					distr_check = float(model.random_p < model.p_value)
					# 4 clusters
					if model.cluster_opt_n <= 0:
						clustering_check = 0
					elif model.cluster_opt_n == clusters_n:
						clustering_check = 1.0
					else:
						clustering_check = max(0, 1 - abs(model.cluster_opt_n - clusters_n) / clusters_n)
					# phasing present
					ph_good = 0
					for name in model.samples:
						if model.samples[name].phase == int(name.split("_")[0]):
							ph_good += 1
					phasing_check = ph_good / len(model.samples)
				
				results[clusters_n][sample_size].append([distr_check, clustering_check, phasing_check, model.random_p, model.cluster_opt_n])
				with open(os.path.join(DIRECTORY, FRESULTS), 'w') as file:
					json.dump(results, file)
				
				print("Distr: %0.2f, Clust: %0.2f, Phasing: %0.2f" % (distr_check, clustering_check, phasing_check))
