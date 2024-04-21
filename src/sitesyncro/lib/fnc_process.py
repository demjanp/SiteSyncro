from sitesyncro.lib.fnc_radiocarbon import (load_calibration_curve)
from sitesyncro.lib.fnc_oxcal import (download_oxcal, r_dates_to_dates, run_model, update_model_by_clustering)
from sitesyncro.lib.fnc_data import (OxCalData, RandomizeData, ClusterData)
from sitesyncro.lib.fnc_sum import (unify_distributions)
from sitesyncro.lib.fnc_simulate import (test_distributions)
from sitesyncro.lib.fnc_cluster import (test_distribution_clustering, cluster_distributions, find_opt_clusters)
from sitesyncro.lib.fnc_phase import (create_earlier_than_matrix, get_groups_and_phases)

import os

def get_curve(curve_name):

	fcurve = os.path.join("OxCal\\bin", curve_name)
	
	if not os.path.isfile(fcurve):
		raise ValueError("Calibration curve not found")
	
	return load_calibration_curve(fcurve)

def proc_model(samples, context_samples, context_area, areas, groups, phases, earlier_than, long_lived, r_dates, model, curve_name, fmodel):
	# samples = [sample, ...]
	# context_samples = {context: [sample, ...], ...}
	# context_area = {context: area, ...}
	# areas = [area, ...]
	# groups = {group: [sample, ...], ...}
	# phases = {group: {sample: phase, ...}, ...}
	# context_phase = {context: phase, ...}
	# earlier_than: matrix[n_samples x n_samples] = [True/False, ...]; sample in row is earlier than sample in column based on stratigraphy
	# long_lived = {sample: True/False, ...}
	# r_dates = {sample: (age, uncertainty), ...}
	# model = 'sequence' / 'contiguous' / 'overlapping' / 'none'
	
	result_path = os.path.dirname(fmodel)
	
	dates = r_dates_to_dates(long_lived, r_dates, curve_name, result_path)
	# dates = {sample: (age, uncertainty, 'R'), sample: (CE from, CE to, 'U'), ...}
	
	oc_data = run_model(dates, phases, model, curve_name, result_path)
	
	# Save results
	oc_data.set_priors(samples, context_samples, context_area, areas, groups, phases, earlier_than, long_lived, r_dates, dates, model, curve_name)
	oc_data.save_json(fmodel)
	
	# oc_data = OxCalData()
	return oc_data
	
def proc_randomize(oc_data, curve_name, frandomized, 
	p_value = 0.05, uncert_base = 15, uniform = False, 
	npass = 100, convergence = 0.99, max_cpus = -1, max_queue_size = 100):
	
	curve = get_curve(curve_name)
	
	# Unify distributions
	distributions = unify_distributions(oc_data.get_posteriors(), curve)
	
	# Randomization testing of the null hypothesis that the observed C-14 dates represent a normal / uniform distribution
	years, sum_obs, distributions_rnd, sums_rnd_lower, sums_rnd_upper, p = test_distributions(distributions, curve,
		   p_value = p_value, uncert_base = uncert_base, uniform = uniform, npass = npass, convergence = convergence,
		   max_cpus = max_cpus, max_queue_size = max_queue_size)
	
	# Store results
	rnd_data = RandomizeData(years, sum_obs, uniform, p, p_value, sums_rnd_lower, sums_rnd_upper)
	rnd_data.save_json(frandomized)
	
	return rnd_data

def proc_clustering(oc_data, curve_name, clusters_n, fclusters,
	p_value = 0.05, uncert_base = 15, uniform = False, 
	npass = 100, convergence = 0.99, max_cpus = -1, max_queue_size = 100):
	
	# clusters_n = number of clusters to form (-1 = automatic) = number of clusters to form (-1 = automatic)
	
	curve = get_curve(curve_name)
	
	# Unify distributions
	distributions = unify_distributions(oc_data.get_posteriors(), curve)
	
	if clusters_n > -1:
		clusters, means, sil = cluster_distributions(distributions, curve[:, 0], clusters_n)
		clusters = {clusters_n: clusters}
		means = {clusters_n: means}
		sils = {clusters_n: sil}
		ps = {clusters_n: 1}
		
	else:
		clusters, means, sils, ps = test_distribution_clustering(distributions, curve,
			uncert_base=uncert_base, uniform=uniform, npass=npass, convergence=convergence,
			max_cpus=max_cpus, max_queue_size=max_queue_size)
	
	# Find optimal number of clusters
	opt_n = find_opt_clusters(clusters, ps, sils, p_value)
		
	# Store results
	clu_data = ClusterData(clusters, means, sils, ps, p_value, opt_n)
	clu_data.save_json(fclusters)
	
	return clu_data

def proc_update_model(oc_data, clu_data, fclu_model):
	# oc_data = OxCalData()
	# clu_data = ClusterData()
	
	result_path = os.path.dirname(fclu_model)
	
	# Update model phasing based on temporal clustering
	oc_clu_data = update_model_by_clustering(oc_data, clu_data, result_path)
	
	# Store results
	oc_clu_data.save_json(fclu_model)
	
	return oc_clu_data

def process(
	samples, context_samples, context_area, areas, context_phase, earlier_than_rel, 
	long_lived, r_dates, model, 
	curve_name = "intcal20.14c", result_path = "result", existing = False, clusters_n = -1,
	p_value = 0.05, uncert_base = 15, uniform = False, 
	npass = 100, convergence = 0.99, max_cpus = -1, max_queue_size = 100):
	
	# samples = [sample, ...]
	# context_samples = {context: [sample, ...], ...}
	# context_area = {context: area, ...}
	# areas = [area, ...]
	# context_phase = {context: phase, ...}
	# earlier_than_rel = {sample: [sample, ...], ...}
	# long_lived = {sample: True/False, ...}
	# r_dates = {sample: (age, uncertainty), ...}
	# model = 'sequence' / 'contiguous' / 'overlapping' / 'none'
	
	if not download_oxcal():
		raise Exception("OxCal not found")
	
	print("Modeling C-14 data based on stratigraphy\n")
	earlier_than = create_earlier_than_matrix(context_phase, earlier_than_rel, samples, context_samples, context_area)
	groups, phases = get_groups_and_phases(earlier_than, samples)
	fmodel = os.path.join(result_path, "model.json")
	if existing and os.path.isfile(fmodel):
		print("Loading from %s\n" % (fmodel))
		oc_data = OxCalData()
		oc_data.load_json(fmodel)
	else:
		oc_data = proc_model(samples, context_samples, context_area, areas, groups, phases, earlier_than, long_lived, r_dates, model, curve_name, fmodel)
		print("Results saved to %s\n" % (fmodel))
	
	print("Testing the distribution of dates for randomness\n")
	frandomized = os.path.join(result_path, "randomized.json.gz")
	if existing and os.path.isfile(frandomized):
		print("Loading from %s\n" % (frandomized))
		rnd_data = RandomizeData()
		rnd_data.load_json(frandomized)
	else:
		rnd_data = proc_randomize(oc_data, curve_name, frandomized, 
			p_value = p_value, uncert_base = uncert_base, uniform = uniform, 
			npass = npass, convergence = convergence, max_cpus = max_cpus, max_queue_size = max_queue_size)
		print("Results saved to %s" % (frandomized))
	
	p = rnd_data.get_p()
	if p < p_value:
		null_hypothesis_txt = "Dates are not %s distributed." % ("uniformly" if uniform else "normally")
	else:
		null_hypothesis_txt = "Dates are %s distributed." % ("uniformly" if uniform else "normally")
	print("\np:", p)
	print(null_hypothesis_txt)
	print()
	
	print("Clustering temporal distributions\n")
	fclusters = os.path.join(result_path, "clusters.json")
	clu_data = proc_clustering(oc_data, curve_name, clusters_n, fclusters,
		p_value = p_value, uncert_base = uncert_base, uniform = uniform, 
		npass = npass, convergence = convergence, max_cpus = max_cpus, max_queue_size = max_queue_size)
	print("\n\nOptimal number of clusters:", clu_data.get_opt_n())
	print("Results saved to %s\n" % (fclusters))
	
	print("Updating model phasing based on temporal clustering\n")
	fclu_model = os.path.join(result_path, "clu_model.json")
	oc_clu_data = proc_update_model(oc_data, clu_data, fclu_model)
	print("Results saved to %s\n" % (fclu_model))
	
	# oc_data = OxCalData()
	# rnd_data = RandomizeData()
	# clu_data = ClusterData()
	# oc_clu_data = OxCalData()
	
	return oc_data, rnd_data, clu_data, oc_clu_data
