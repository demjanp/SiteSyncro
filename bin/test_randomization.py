from sitesyncro import Model
from sitesyncro.utils.fnc_radiocarbon import (get_curve)

import json
import numpy as np

uncertainty_base = 30
gap = 300
center_date = 3500
dates_n = {'small': 10, 'large': 50}
iterations = 10

fresults = "test_results.json"

def sim_event(cal_age_bp, n_dates, curve):
	
	c14_age, error = curve[np.argmin(np.abs(curve[:,0] - cal_age_bp))][1:3]
	return np.random.normal(c14_age, error, n_dates)

def sim_clustered(cal_bp_mean, gap, clusters_n, dates_n, uncertainty_base, curve):
	# Simulate dates_n dates from clusters_n of events with gaps of gap years around cal_bp_mean
	# returns [[phase, c14age, uncert], ...]
	
	if clusters_n == 0:
		# Generate uniformly distributed events around cal_bp_mean +-gap
		means = np.random.uniform(cal_bp_mean - gap, cal_bp_mean + gap, dates_n)
		means.sort()
		dates = []
		for cal_age in means:
			c14age = sim_event(cal_age, 1, curve)[0]
			dates.append([1, c14age, np.random.normal(uncertainty_base, 5)])
		return dates
	
	if clusters_n == 1:
		dates = []
		for c14age in sim_event(cal_bp_mean, dates_n, curve):
			dates.append([1, c14age, np.random.normal(uncertainty_base, 5)])
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
		for c14age in sim_event(mean, ns[i], curve):
			dates.append([i + 1, c14age, np.random.normal(uncertainty_base, 5)])
	return dates

if __name__ == '__main__':
	
	curve = get_curve()
	
	results = {}
	for clusters_n in [0, 1, 4]:
		results[clusters_n] = {}
		for sample_size in dates_n:
			results[clusters_n][sample_size] = []
			directory = "sim_%d_%s" % (clusters_n, sample_size)
			for iter in range(iterations):
				print("\nClusters: %d, Sample size: %s, Iter. %d/%d\n" % (clusters_n, sample_size, iter + 1, iterations))
				dates = sim_clustered(center_date, gap, clusters_n, dates_n[sample_size], uncertainty_base, curve)
				model = Model(directory = directory, uniform = (clusters_n == 0), overwrite = True)
				n = 1
				for phase, age, uncert in dates:
					name = "%d_%d" % (phase, n)
					n += 1
					model.add_sample(name, age, uncert, context = "A.%d" % (n), area = "A", area_excavation_phase = 1)
				model.process_phasing()
				model.process_randomization()
				model.process_clustering()
				distr_check = None
				clustering_check = None
				if clusters_n == 0:
					# uniform distribution of events
					distr_check = float(model.random_p >= model.p_value)
					# no clustering of events
					clustering_check = float(model.cluster_opt_n < 2)
				elif clusters_n == 1:
					# normal distribution of events
					distr_check = float(model.random_p >= model.p_value)
					# no clustering of events
					clustering_check = float(model.cluster_opt_n < 2)
				else:
					# non-normal distribution of events
					distr_check = float(model.random_p < model.p_value)
					# 4 clusters
					if model.cluster_opt_n <= 0:
						clustering_check = 0
					elif model.cluster_opt_n == clusters_n:
						clustering_check = 1.0
					else:
						clustering_check = max(0, 1 - abs(model.cluster_opt_n - clusters_n) / clusters_n)
				
				results[clusters_n][sample_size].append([distr_check, clustering_check, model.random_p, model.cluster_opt_n])
				with open(fresults, 'w') as file:
					json.dump(results, file)