import numpy as np
from collections import defaultdict

def calc_sum(distributions):
	# distributions = [[p, ...], ...]
	summed = np.zeros(len(distributions[0]), dtype = np.float64)
	for dist in distributions:
		summed += np.array(dist, dtype = np.float64)
	s = summed.sum()
	if s != 0:
		summed /= s
	return summed

def calc_mean(years, distribution):
	if not distribution.sum():
		return None
	return float(np.average(years, weights=distribution))

def calc_mean_std(years, distribution):
	mean = calc_mean(years, distribution)
	if mean is None:
		return None, None
	std = np.sqrt(np.average((years - mean) ** 2, weights=distribution))
	return mean, std

def calc_range(years, distribution, p = 0.9545, p_threshold = 0.0001, max_multiplier = 32):
	# returns [from, to] in calendar years BP
	
	def _find_rng(values, weights, v_min, v_max, mul):
		
		vs = np.linspace(values.min(), values.max(), values.shape[0] * mul)
		weights = np.interp(vs, values, weights)
		values = vs
		weights /= weights.sum()
		idx_start = max(0, int(np.where(values <= v_min)[0].max()) - 1)
		idx_end = min(values.shape[0] - 1, int(np.where(values >= v_max)[0].min()) + 1)
		values = values[idx_start:idx_end]
		weights = weights[idx_start:idx_end]
		levels = np.unique(weights)
		levels.sort()
		lvl_prev = 0
		for lvl in levels[::-1]:
			mask = (weights >= lvl)
			p_sum = weights[mask].sum()
			if p_sum >= p:
				mask = (weights >= lvl_prev)
				p_sum = weights[mask].sum()
				idxs = np.where(mask)[0]
				return values[idxs.min()], values[idxs.max()], p - p_sum
			lvl_prev = lvl
		idxs = np.where(weights > 0)[0]
		return values[idxs.min()], values[idxs.max()], 1 - p	
	
	years = np.array(years)
	distribution = np.array(distribution)
	if not distribution.sum():
		return [None, None]
	v_min, v_max = years.min(), years.max()
	p_diff_opt = np.inf
	v_opt = [v_min, v_max]
	mul = 2
	while True:
		v_min, v_max, p_diff = _find_rng(years, distribution, v_min, v_max, mul)
		if p_diff >= p_diff_opt:
			mul *= 2
		if p_diff < p_diff_opt:
			p_diff_opt = p_diff
			v_opt = [v_min, v_max]
		if p_diff <= p_threshold:
			break
		if mul > max_multiplier:
			break
	
	return [float(v_opt[1]), float(v_opt[0])]

def calc_percentiles(distributions, perc_lower, perc_upper):
	# distributions = [[p, ...], ...]
	combined = np.array(distributions, dtype=np.float64)
	percentiles = np.percentile(combined, [perc_lower, perc_upper], axis=0)
	return percentiles

def samples_to_distributions(samples):
	# Get a list of probability distributions from samples
	# Combine samples from the same context by summation
	# Use only posterior (modeled) distributions for long-lived samples and outliers
	# samples = [Sample, ...]
	#
	# returns distributions, samples_collected, joined
	#	distributions = [[p, ...], ...]
	#	samples_collected = [sample name, ...] ordered by distributions
	#	joined = {combined_name: [sample name, ...], ...}
	distributions_dict = {}  # {name: np.array([p, ...]), ...}
	contexts = defaultdict(list)
	for sample in samples:
		context = sample.context
		if context is None:
			context = sample.name
		if sample.is_modeled:
			distributions_dict[sample.name] = sample.posterior
			contexts[context].append(sample.name)
		elif sample.outlier:
			continue
		elif sample.is_calibrated and not (sample.long_lived or sample.outlier):
			distributions_dict[sample.name] = sample.likelihood
			contexts[context].append(sample.name)
	joined = {}  # {combined_sample: [sample_name, ...], ...}
	samples_collected = []
	distributions = []
	for context, sample_names in contexts.items():
		if len(sample_names) > 1:
			combined_sample = '+'.join(sample_names)
			joined[combined_sample] = sample_names
			distributions.append(calc_sum([distributions_dict[name] for name in sample_names]))
			samples_collected.append(combined_sample)
		else:
			distributions.append(distributions_dict[sample_names[0]])
			samples_collected.append(sample_names[0])
	
	return distributions, samples_collected, joined
