from collections import defaultdict
from typing import List, Dict

import numpy as np


def calc_sum(distributions: List[np.ndarray]) -> np.ndarray:
	"""
	Calculate the sum of the given distributions.

	Parameters:
	distributions (List[np.ndarray]): List of distributions.

	Returns:
	np.ndarray: The summed distribution.
	"""
	
	summed = np.zeros(len(distributions[0]), dtype=np.float64)
	for dist in distributions:
		summed += np.array(dist, dtype=np.float64)
	s = summed.sum()
	if s != 0:
		summed /= s
	return summed


def calc_mean(years: np.ndarray, distribution: np.ndarray) -> float or None:
	"""
	Calculate the mean of the given distribution.

	Parameters:
	years (np.ndarray): Array of years.
	distribution (np.ndarray): Array of probabilities for each year.

	Returns:
	float: The mean of the distribution, or None if the distribution sum is zero.
	"""
	if not distribution.sum():
		return None
	return float(np.average(years, weights=distribution))


def calc_mean_std(years: np.ndarray, distribution: np.ndarray) -> (float, float) or (None, None):
	"""
	Calculate the mean and standard deviation of the given distribution.

	Parameters:
	years (np.ndarray): Array of years.
	distribution (np.ndarray): Array of probabilities for each year.

	Returns:
	(mean, std): The mean and standard deviation of the distribution, or (None, None) if the mean is None.
	"""
	mean = calc_mean(years, distribution)
	if mean is None:
		return None, None
	std = np.sqrt(np.average((years - mean) ** 2, weights=distribution))
	return mean, std


def calc_range(years: np.ndarray, distribution: np.ndarray, p: float = 0.9545, p_threshold: float = 0.0001,
               max_multiplier: int = 32) -> [float, float] or [None, None]:
	"""
	Calculate the range of the given distribution.

	Parameters:
	years (np.ndarray): Array of years.
	distribution (np.ndarray): Array of probabilities for each year.
	p (float): The desired percentile. Default is 0.9545.
	p_threshold (float): The threshold for the percentile difference. Default is 0.0001.
	max_multiplier (int): The maximum multiplier for the number of values. Default is 32.

	Returns:
	[from, to]: The range of the distribution, or [None, None] if the distribution sum is zero.
	"""
	
	def _find_rng(values: np.ndarray, weights: np.ndarray, v_min: float, v_max: float, mul: int):
		
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


def calc_percentiles(distributions: List[np.ndarray] or np.ndarray, perc_lower: float, perc_upper: float) -> np.ndarray:
	"""
	Calculate the lower and upper percentiles of the given distributions.

	Parameters:
	distributions (List[np.ndarray]): List of distributions.
	perc_lower (float): The lower percentile.
	perc_upper (float): The upper percentile.

	Returns:
	np.ndarray: A 2D array of the lower and upper percentiles.
	"""
	combined = np.array(distributions, dtype=np.float64)
	percentiles = np.percentile(combined, [perc_lower, perc_upper], axis=0)
	return percentiles


def samples_to_distributions(samples: List[object]) -> (List[np.ndarray], List[str], Dict[str, List[str]]):
	"""
	Convert a list of samples to a list of probability distributions.
	
	Combine samples from the same context by summation. Use only posterior (modeled) distributions for long-lived samples and outliers.

	Parameters:
	samples (List[Sample]): List of Sample objects.

	Returns:
	(distributions, samples_collected, joined)
		- distributions: List of probability distributions.
		- samples_collected: List of names of the samples or combined samples.
		- joined: {combined_name: [sample name, ...], ...}; Dictionary mapping combined samples to their original samples.
	"""
	
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
