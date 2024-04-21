import numpy as np
from scipy.interpolate import interp1d

def unify_distributions(distributions, curve):
	# distributions = {name: [p, ...], ...}
	# Years are converted from CE to BP
	
	# Modify distributions so all lists have the same length
	all_years = curve[:,0]
	result = {}
	for name in distributions:
		dist, years = distributions[name]['prob'], distributions[name]['years'] # years are CE
		years = 1950 - np.array(years)
		
		# Check if the years are in the correct order
		if np.sign(years[-1] - years[0]) != np.sign(all_years[-1] - all_years[0]):
			years = years[::-1]
			dist = dist[::-1]
		
		years = np.concatenate(([all_years[0]], years, [all_years[-1]]))
		dist = np.concatenate(([0], dist, [0]))
		s = dist.sum()
		if s > 0:
			dist /= s
		new_prob = interp1d(years, dist, kind="linear")(all_years)
		s = new_prob.sum()
		if s > 0:
			new_prob /= s
		result[name] = new_prob
	
	return result

def sum_distributions(distributions):
	# distributions = [[p, ...], ...]
	
	# Create a list of all probabilities
	summed = np.zeros(len(distributions[0]), dtype = float)
	
	# Sum probabilities
	for dist in distributions:
		summed += np.array(dist, dtype = float)
	
	# Normalize probabilities
	s = summed.sum()
	if s != 0:
		summed /= s
	
	return summed

def calc_mean_std(distribution):
	# distribution = [years, prob]
	
	years, prob = distribution
	
	mean = np.average(years, weights=prob)
	std = np.sqrt(np.average((years - mean) ** 2, weights=prob))
	
	return mean, std

def calc_range(distribution, p = 0.9545, p_threshold = 0.0001, max_multiplier = 32):
	# distribution = [years, prob]
	
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
	
	years, prob = distribution
	years = np.array(years)
	prob = np.array(prob)
	v_min, v_max = years.min(), years.max()
	p_diff_opt = np.inf
	v_opt = [v_min, v_max]
	mul = 2
	while True:
		v_min, v_max, p_diff = _find_rng(years, prob, v_min, v_max, mul)
		if p_diff >= p_diff_opt:
			mul *= 2
		if p_diff < p_diff_opt:
			p_diff_opt = p_diff
			v_opt = [v_min, v_max]
		if p_diff <= p_threshold:
			break
		if mul > max_multiplier:
			break
	
	return v_opt

def calc_percentiles(distributions, perc_lower, perc_upper):
	# distributions = [[p, ...], ...]
	
	combined = np.array(distributions)
	
	percentiles = np.percentile(combined, [perc_lower, perc_upper], axis=0)
	
	return percentiles
