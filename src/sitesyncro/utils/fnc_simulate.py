from typing import Any, List

import numpy as np
from scipy.stats import norm
from tqdm import tqdm

from sitesyncro.utils.fnc_mp import (process_mp)
from sitesyncro.utils.fnc_radiocarbon import (calibrate)
from sitesyncro.utils.fnc_stat import (calc_sum, calc_mean_std, calc_percentiles, samples_to_distributions)


def generate_random_distributions(dates_n: int, t_mean: float, t_std: float, uncertainties: List[float],
                                  uncertainty_base: float, curve: np.ndarray, uniform: bool,
                                  max_iterations: int = 10000) -> List[np.ndarray]:
	"""
	Generate sets of randomized distributions based on observed distributions.

	Parameters:
	dates_n (int): Number of dates to generate.
	t_mean (float): Mean of the calendar ages.
	t_std (float): Standard deviation of the calendar ages.
	uncertainties (list): Pool of uncertainties to simulate C-14 dates.
	uncertainty_base (float): Base uncertainty to calculate uncertainty based on C-14 age if pool is not available.
	curve (np.ndarray): A 2D array containing the calibration curve data. Each row represents a calendar year BP, C-14 year, and uncertainty.
	uniform (bool): Flag indicating whether to use a uniform distribution for the calendar ages.
	max_iterations (int): Maximum number of iterations for the adjustment process. Default is 10000.

	Returns:
	List[np.ndarray]: List of generated distributions.
	"""
	
	def _sim_age(t_mean, t_std, uniform, curve):
		# Generate a random calendar age
		if uniform:
			t = np.random.uniform(t_mean - t_std, t_mean + t_std)
		else:
			t = np.random.normal(t_mean, t_std)
		# Find the closest index in the calibration curve
		idx = np.argmin(np.abs(curve[:, 0] - t))
		# Get the corresponding radiocarbon age
		age = curve[idx, 1]
		return age
	
	def _gen_distributions(dates_n, uncertainties, uncertainty_base, t_mean, t_std, uniform, curve):
		dates = []
		distributions = []
		for i in range(dates_n):
			# Generate a random calendar age
			age = _sim_age(t_mean, t_std, uniform, curve)
			# Calculate the standard deviation for the date
			if uncertainties:
				stdev = np.random.choice(uncertainties)
			else:
				stdev = uncertainty_base * np.exp(age / (2 * 8033))
			# Append the date and its standard deviation to the dates list
			dates.append([age, stdev])
			# Calibrate the date and append the distribution to the distributions list
			distribution = calibrate(age, stdev, curve)
			distributions.append(distribution)
		return np.array(dates), distributions
	
	def _analyze(distributions, curve):
		# Sum the distributions
		dist_summed = calc_sum(distributions)
		# Calculate the mean and standard deviation of the summed distribution
		mean_sum, std_sum = calc_mean_std(curve[:, 0], dist_summed)
		return mean_sum, std_sum
	
	dates, distributions = _gen_distributions(dates_n, uncertainties, uncertainty_base, t_mean, t_std, uniform, curve)
	
	curve_min = curve[:, 1].min()
	curve_max = curve[:, 1].max()
	
	# Adjust the dates until the mean and standard deviation of the summed distribution are close to the desired values
	iteration = 0
	scaling_factor = 0.5
	c_threshold = t_std * 0.001
	m, s = 0, 0
	while True:
		# Adjust mean
		reset = False
		while not reset:
			m, _ = _analyze(distributions, curve)
			if abs(m - t_mean) <= c_threshold:
				break
			d = scaling_factor*(t_mean - m)
			for i in range(len(dates)):
				dates[i][0] += d
				if (dates[i][0] < curve_min) or (dates[i][0] > curve_max):
					reset = True
					break
				distributions[i] = calibrate(dates[i][0], dates[i][1], curve)
		# Adjust stdev
		while not reset:
			m, s = _analyze(distributions, curve)
			if abs(s - t_std) <= c_threshold:
				break
			# Scale dates around center
			sf = 1 + scaling_factor*((t_std / s) - 1)
			d_m = dates[:,0].mean()
			for i in range(len(dates)):
				dates[i][0] = d_m + (dates[i][0] - d_m)*sf
				if (dates[i][0] < curve_min) or (dates[i][0] > curve_max):
					reset = True
					break
				distributions[i] = calibrate(dates[i][0], dates[i][1], curve)
		
		iteration += 1
		if (iteration > max_iterations) or reset:
			# Re-initialize
			dates, distributions = _gen_distributions(dates_n, uncertainties, uncertainty_base, t_mean, t_std, uniform, curve)
			iteration = 0
			continue
		
		if max(abs(m - t_mean), abs(s - t_std)) <= c_threshold:
			break
	
	return distributions


def worker_fnc(params: Any, dates_n: int, t_mean: float, t_std: float, uncertainties: List[float],
               uncertainty_base: float, curve: np.ndarray, uniform: bool) -> np.ndarray:
	return calc_sum(
		generate_random_distributions(dates_n, t_mean, t_std, uncertainties, uncertainty_base, curve, uniform))


def collect_fnc(data: np.ndarray, results: List[np.ndarray], pbar: tqdm) -> None:
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
	
	sum_obs = calc_sum(distributions)
	years = model.curve[:, 0]
	t_mean, t_std = calc_mean_std(years, sum_obs)
	
	sums = []
	sums_prev = None
	c = 0
	todo = model.npass
	with tqdm(total=todo * 2) as pbar:
		pbar.set_description("Convergence: %0.3f" % (c))
		while True:
			process_mp(worker_fnc, range(max(4, (todo - len(sums)) + 1)),
			           [dates_n, t_mean, t_std, model.uncertainties, model.uncertainty_base, model.curve,
			            model.uniform],
			           collect_fnc=collect_fnc, collect_args=[sums, pbar],
			           max_cpus=max_cpus, max_queue_size=max_queue_size)
			'''
			for i in range(max(4, (todo - len(sums)) + 1)):
				sums.append(worker_fnc(i, dates_n, t_mean, t_std, model.uncertainties, model.uncertainty_base, model.curve, model.uniform))
				pbar.update(1)
			'''  # DEBUG
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
	
	mask = (sums_rnd.std(axis=0) > 0)
	p = (1 - norm(sums_rnd[:, mask].mean(axis=0), sums_rnd[:, mask].std(axis=0)).cdf(sum_obs[mask])).min()
	
	perc_lower = (model.p_value * 100) / 2
	perc_upper = 100 - perc_lower
	sums_rnd_lower, sums_rnd_upper = calc_percentiles(sums_rnd, perc_lower, perc_upper)
	
	return sum_obs, sums_rnd_lower, sums_rnd_upper, p
