from typing import Any, List

import numpy as np
from scipy.stats import norm
from tqdm import tqdm

from sitesyncro.utils.fnc_mp import (process_mp)
from sitesyncro.utils.fnc_radiocarbon import (calibrate)
from sitesyncro.utils.fnc_stat import (calc_sum, calc_mean_std, calc_range, calc_range_approx, calc_percentiles, samples_to_distributions)


def get_params(distributions, curve, uniform):
	"""
	Calculate parameters of the summed distributions
	
	Returns:
	(mean_sum, std_sum)
	mean_sum: Mean (normal model) or center of 1-sigma range (uniform model)
	std_sum: Standard deviation (normal model) or 1/2 range (uniform model)
	"""
	# Sum the distributions
	sum_dist = calc_sum(distributions)
	if uniform:
		rng = calc_range(curve[:,0], sum_dist, p = 0.6827)
		mean_sum = np.mean(rng)
		std_sum = abs(np.diff(rng)[0]) / 2
	else:
		# Calculate the mean and standard deviation of the summed distribution
		mean_sum, std_sum = calc_mean_std(curve[:, 0], sum_dist)
	return mean_sum, std_sum


def gen_random_dists_uniform(dates_n: int, t_mean: float, t_std: float, uncertainties: List[float], uncertainty_base: float, curve: np.ndarray) -> List[np.ndarray]:
	
	distributions = []
	sum_dists = np.zeros(curve.shape[0], dtype = float)
	while len(distributions) < dates_n:
		cal_age = np.random.uniform(t_mean - t_std*2, t_mean + t_std*2)
		age = curve[np.argmin(np.abs(curve[:, 0] - cal_age)), 1]
		if uncertainties:
			uncert = np.random.choice(uncertainties)
		else:
			uncert = uncertainty_base * np.exp(age / (2 * 8033))
		dist = calibrate(age, uncert, curve)
		s_d = sum_dists + dist
		s = s_d.sum()
		if s > 0:
			s_d /= s
		rng = calc_range_approx(curve[:,0], s_d, p = 0.6827)
		if (rng[0] <= t_mean + t_std) and (rng[1] >= t_mean - t_std):
			distributions.append(dist)
			sum_dists += dist
	return distributions


def gen_random_dists_normal(dates_n: int, t_mean: float, t_std: float, uncertainties: List[float],
								  uncertainty_base: float, curve: np.ndarray) -> List[np.ndarray]:
	"""
	Generate sets of randomized distributions based on observed distributions.

	Parameters:
	dates_n (int): Number of dates to generate.
	t_mean (float): Mean of the calendar ages.
	t_std (float): Standard deviation of the calendar ages.
	uncertainties (list): Pool of uncertainties to simulate C-14 dates.
	uncertainty_base (float): Base uncertainty to calculate uncertainty based on C-14 age if pool is not available.
	curve (np.ndarray): A 2D array containing the calibration curve data. Each row represents a calendar year BP, C-14 year, and uncertainty.
	
	Returns:
	List[np.ndarray]: List of generated distributions.
	"""
	
	def _get_params(dates, curve):
		
		sum_dist = calc_sum([calibrate(age, uncert, curve) for age, uncert in dates])
		mean_sum, std_sum = calc_mean_std(curve[:, 0], sum_dist)
		return mean_sum, std_sum
	
	def _gen_dates(dates_n, t_mean, t_std, uncertainties, uncertainty_base, curve):
		
		cal_ages = np.random.normal(t_mean, t_std, dates_n)
		dates = []
		for cal_age in cal_ages:
			age = curve[np.argmin(np.abs(curve[:, 0] - cal_age)), 1]
			if uncertainties:
				uncert = np.random.choice(uncertainties)
			else:
				uncert = uncertainty_base * np.exp(age / (2 * 8033))
			dates.append([age, uncert])
		return np.array(dates)
	
	dates = _gen_dates(dates_n, t_mean, t_std, uncertainties, uncertainty_base, curve)
	
	# Adjust the dates until the mean and standard deviation of the summed distribution are close to the desired values
	c_threshold = t_std * 0.001
	m, s = 0, 0
	curve_min = curve[:, 1].min()
	curve_max = curve[:, 1].max()
	while True:
		# Adjust mean
		reset = False
		scaling_factor = 1
		while not reset:
			m, _ = _get_params(dates, curve)
			if abs(m - t_mean) <= c_threshold:
				break
			dates[:,0] += scaling_factor*(t_mean - m)
			scaling_factor *= 0.99
			if (scaling_factor < 0.05) or (dates[:,0].min() < curve_min) or (dates[:,0].max() > curve_max):
				reset = True
		# Adjust stdev
		scaling_factor = 1
		while not reset:
			_, s = _get_params(dates, curve)
			if abs(s - t_std) <= c_threshold:
				break
			# Scale dates around center
			sf = 1 + scaling_factor*((t_std / s) - 1)
			d_m = dates[:,0].mean()
			dates[:,0] = d_m + (dates[:,0] - d_m)*sf
			scaling_factor *= 0.99
			if (scaling_factor < 0.05) or (dates[:,0].min() < curve_min) or (dates[:,0].max() > curve_max):
				reset = True
		if reset:
			# Re-initialize
			dates = _gen_dates(dates_n, t_mean, t_std, uncertainties, uncertainty_base, curve)
			continue
		
		m, s = _get_params(dates, curve)
		if max(abs(m - t_mean), abs(s - t_std)) <= c_threshold:
			break
	distributions = [calibrate(age, uncert, curve) for age, uncert in dates]
	
	return distributions


def generate_random_distributions(dates_n: int, t_mean: float, t_std: float, uncertainties: List[float],
								  uncertainty_base: float, curve: np.ndarray, uniform: bool) -> List[np.ndarray]:
	
	if uniform:
		return gen_random_dists_uniform(dates_n, t_mean, t_std, uncertainties, uncertainty_base, curve)
	return gen_random_dists_normal(dates_n, t_mean, t_std, uncertainties, uncertainty_base, curve)


def test_distributions_worker(params: Any, dates_n: int, t_mean: float, t_std: float, uncertainties: List[float],
			   uncertainty_base: float, curve: np.ndarray, uniform: bool) -> np.ndarray:
	
	return calc_sum(generate_random_distributions(dates_n, t_mean, t_std, uncertainties, uncertainty_base, curve, uniform))


def test_distributions_collect(data: np.ndarray, results: List[np.ndarray], pbar: tqdm) -> None:
	# data = dist_sum
	# dist_sum = [p, ...]
	pbar.update(1)
	results.append(data)


def test_distributions(model: object, max_cpus: int = -1, max_queue_size: int = -1) -> (
		np.ndarray, np.ndarray, np.ndarray, float):
	"""
	Test if the samples in the model are distributed uniformly / normally in time

	Parameters:
	model (Model): The Model object containing the samples.
	max_cpus (int): Maximum number of CPUs to use for multiprocessing. Default is -1 (use all available CPUs).
	max_queue_size (int): Maximum size of the queue for multiprocessing. Default is -1 (no limit).

	Returns:
	(sum_obs, sums_rnd_lower, sums_rnd_upper, p)
		- sum_obs: Sum of the observed probability distributions of the sample dates.
		- sums_rnd_lower: Lower bound of the randomized sum of the probability distributions.
		- sums_rnd_upper: Upper bound of the randomized sum of the probability distributions.
		- p: P-value for the test.
	"""
	
	distributions, _, _ = samples_to_distributions(model.samples.values())
	
	dates_n = len(distributions)
	
	if dates_n < 3:
		raise Exception("Insufficient number samples for testing.")
	
	print("Testing %d distributions" % dates_n)
	
	t_mean, t_std = get_params(distributions, model.curve, model.uniform)
	
	sums = []
	sums_prev = None
	c = 0
	todo = model.npass
	with tqdm(total=todo * 2) as pbar:
		pbar.set_description("Convergence: %0.3f" % (c))
		while True:
			n_dists = max(4, (todo - len(sums)) + 1)
			if n_dists > 50:
				process_mp(test_distributions_worker, range(n_dists),
						[dates_n, t_mean, t_std, model.uncertainties, model.uncertainty_base, model.curve, model.uniform],
						collect_fnc=test_distributions_collect, collect_args=[sums, pbar],
						max_cpus=max_cpus, max_queue_size=max_queue_size)
			else:
				for i in range(n_dists):
					sums.append(test_distributions_worker(i, dates_n, t_mean, t_std, model.uncertainties, model.uncertainty_base, model.curve, model.uniform))
					pbar.update(1)
			
			if len(sums) >= todo:
				sums_m = np.array(sums).mean(axis=0)
				if sums_prev is not None:
					c = ((sums_prev * sums_m).sum() ** 2) / ((sums_prev ** 2).sum() * (sums_m ** 2).sum())
					pbar.set_description("Convergence: %0.3f" % (c))
				sums_prev = sums_m.copy()
				if c >= model.convergence:
					break
				todo *= 2
				pbar.total = max(todo, model.npass * 2)
				pbar.refresh()
	
	sums_rnd = np.array(sums)
	
	sum_obs = calc_sum(distributions)
	mask = (sums_rnd.std(axis=0) > 0)
	p = (1 - norm(sums_rnd[:, mask].mean(axis=0), sums_rnd[:, mask].std(axis=0)).cdf(sum_obs[mask])).min()
	
	perc_lower = (model.p_value * 100) / 2
	perc_upper = 100 - perc_lower
	sums_rnd_lower, sums_rnd_upper = calc_percentiles(sums_rnd, perc_lower, perc_upper)
	
	return sum_obs, sums_rnd_lower, sums_rnd_upper, p

