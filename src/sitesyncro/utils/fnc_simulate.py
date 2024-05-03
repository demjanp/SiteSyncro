from typing import Any, List

import numpy as np
from scipy.stats import norm
from tqdm import tqdm

from sitesyncro.utils.fnc_mp import (process_mp)
from sitesyncro.utils.fnc_radiocarbon import (calibrate)
from sitesyncro.utils.fnc_stat import (calc_sum, calc_mean_std, calc_range, calc_range_approx, calc_percentiles, samples_to_distributions)


def get_params(distributions, curve, uniform, approx = False):
	"""
	Calculate parameters of the summed distributions
	
	Returns:
	(mean_sum, std_sum)
	mean_sum: Mean (normal model) or center of range (uniform model)
	std_sum: Standard deviation (normal model) or 1/2 range (uniform model)
	"""
	# Sum the distributions
	sum_dist = calc_sum(distributions)
	if uniform:
		range_fnc = calc_range_approx if approx else calc_range
		rng = range_fnc(curve[:,0], sum_dist)
		mean_sum = np.mean(rng)
		std_sum = abs(np.diff(rng)[0]) / 2
	else:
		# Calculate the mean and standard deviation of the summed distribution
		mean_sum, std_sum = calc_mean_std(curve[:, 0], sum_dist)
	return mean_sum, std_sum


def get_range_pool_worker(params, curve, t1, t2):
	
	age, uncert = params
	dist = calibrate(age, uncert, curve)
	r1, r2 = calc_range(curve[:,0], dist)
	if (r1 <= t1) and (r2 >= t2):
		return [r1, r2, age, uncert]
	return None


def get_range_collect(data, range_pool, pbar):
	
	pbar.update(1)
	if data is not None:
		range_pool.append(data)


def get_range_pool(t_mean: float, t_std: float, uncertainties: List[float], uncertainty_base: float, curve: np.ndarray, max_cpus: int = -1, max_queue_size: int = -1):
	
	print("Generating pool of dating ranges")
	t1 = t_mean + t_std
	t2 = t_mean - t_std
	ages = np.unique(curve[(curve[:,0] <= t1) & (curve[:,0] >= t2)][:,1])
	uncerts = []
	if uncertainties:
		uncerts = np.unique(np.round(np.array(uncertainties)/2)*2).tolist()
	params_list = []
	for age in ages:
		if uncerts:
			for uncert in uncerts:
				params_list.append([age, uncert])
		else:
			params_list.append([age, uncertainty_base * np.exp(age / (2 * 8033))])
	
	range_pool = []
	with tqdm(total=len(params_list)) as pbar:
		process_mp(get_range_pool_worker, params_list, [curve, t1, t2],
					collect_fnc=get_range_collect, collect_args=[range_pool, pbar],
					max_cpus=max_cpus, max_queue_size=max_queue_size)
	
	return np.array(range_pool)


def gen_random_dists_uniform(dates_n: int, t_mean: float, t_std: float, uncertainties: List[float], curve: np.ndarray, range_pool: np.ndarray) -> List[np.ndarray]:
	
	distributions = []
	while len(distributions) < dates_n:
		t = np.random.uniform(t_mean - t_std, t_mean + t_std)
		slice = range_pool[(range_pool[:,0] >= t) & (range_pool[:,1] <= t)]
		if not slice.size:
			continue
		if uncertainties:
			uncert = np.random.choice(uncertainties)
			slice = slice[np.abs(slice[:,3] - uncert) == np.abs(slice[:,3] - uncert).min()]
			if not slice.size:
				continue
		_, _, age, uncert = slice[np.random.randint(len(slice))]
		distributions.append(calibrate(age, uncert, curve))
	
	return distributions


def gen_random_dists_normal(dates_n: int, t_mean: float, t_std: float, uncertainties: List[float],
								  uncertainty_base: float, curve: np.ndarray, max_iterations: int = 10000) -> List[np.ndarray]:
	"""
	Generate sets of randomized distributions based on observed distributions.

	Parameters:
	dates_n (int): Number of dates to generate.
	t_mean (float): Mean of the calendar ages.
	t_std (float): Standard deviation of the calendar ages.
	uncertainties (list): Pool of uncertainties to simulate C-14 dates.
	uncertainty_base (float): Base uncertainty to calculate uncertainty based on C-14 age if pool is not available.
	curve (np.ndarray): A 2D array containing the calibration curve data. Each row represents a calendar year BP, C-14 year, and uncertainty.
	max_iterations (int): Maximum number of iterations for the adjustment process. Default is 10000.

	Returns:
	List[np.ndarray]: List of generated distributions.
	"""
	
	def _sim_age(t_mean, t_std, curve):
		# Generate a random calendar age
		t = np.random.normal(t_mean, t_std)
		# Find the closest index in the calibration curve
		idx = np.argmin(np.abs(curve[:, 0] - t))
		# Get the corresponding radiocarbon age
		age = curve[idx, 1]
		return age
	
	def _gen_distributions(dates_n, uncertainties, uncertainty_base, t_mean, t_std, curve):
		dates = []
		distributions = []
		for i in range(dates_n):
			# Generate a random calendar age
			age = _sim_age(t_mean, t_std, curve)
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
	
	dates, distributions = _gen_distributions(dates_n, uncertainties, uncertainty_base, t_mean, t_std, curve)
	
	curve_min = curve[:, 1].min()
	curve_max = curve[:, 1].max()
	
	# Adjust the dates until the mean and standard deviation of the summed distribution are close to the desired values
	iteration = 0
	c_threshold = t_std * 0.001
	m, s = 0, 0
	while True:
		# Adjust mean
		reset = False
		scaling_factor = 1
		while not reset:
			m, _ = get_params(distributions, curve, False)
			if abs(m - t_mean) <= c_threshold:
				break
			sfactor = scaling_factor
			d = sfactor*(t_mean - m)
			for i in range(len(dates)):
				dates[i][0] += d
				if (dates[i][0] < curve_min) or (dates[i][0] > curve_max):
					reset = True
					break
				distributions[i] = calibrate(dates[i][0], dates[i][1], curve)
			scaling_factor *= 0.99
			if scaling_factor < 0.01:
				reset = True
		# Adjust stdev
		scaling_factor = 1
		while not reset:
			_, s = get_params(distributions, curve, False)
			if abs(s - t_std) <= c_threshold:
				break
			# Scale dates around center
			sfactor = scaling_factor
			sf = 1 + sfactor*((t_std / s) - 1)
			d_m = dates[:,0].mean()
			for i in range(len(dates)):
				dates[i][0] = d_m + (dates[i][0] - d_m)*sf
				if (dates[i][0] < curve_min) or (dates[i][0] > curve_max):
					reset = True
					break
				distributions[i] = calibrate(dates[i][0], dates[i][1], curve)
			scaling_factor *= 0.99
			if scaling_factor < 0.01:
				reset = True
		iteration += 1
		if (iteration > max_iterations) or reset:
			# Re-initialize
			dates, distributions = _gen_distributions(dates_n, uncertainties, uncertainty_base, t_mean, t_std, curve)
			iteration = 0
			continue
		
		m, s = get_params(distributions, curve, False)
		if max(abs(m - t_mean), abs(s - t_std)) <= c_threshold:
			break
	
	return distributions


def generate_random_distributions(dates_n: int, t_mean: float, t_std: float, uncertainties: List[float],
								  uncertainty_base: float, curve: np.ndarray, uniform: bool, range_pool: np.ndarray = None,
								  max_iterations: int = 10000, max_cpus: int = -1, max_queue_size: int = -1) -> List[np.ndarray]:
	
	if uniform:
		if range_pool is None:
			range_pool = get_range_pool(t_mean, t_std, uncertainties, uncertainty_base, curve, max_cpus=max_cpus, max_queue_size=max_queue_size)
			if not range_pool.size:
				raise Exception("Could not generate random dates")
		return gen_random_dists_uniform(dates_n, t_mean, t_std, uncertainties, curve, range_pool)
	return gen_random_dists_normal(dates_n, t_mean, t_std, uncertainties, uncertainty_base, curve, max_iterations)


def test_distributions_worker(params: Any, dates_n: int, t_mean: float, t_std: float, uncertainties: List[float],
			   uncertainty_base: float, curve: np.ndarray, uniform: bool, range_pool: np.ndarray) -> np.ndarray:
	
	return calc_sum(generate_random_distributions(dates_n, t_mean, t_std, uncertainties, uncertainty_base, curve, uniform, range_pool))


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
	
	range_pool = None
	if model.uniform:
		range_pool = model._get_range_pool(t_mean, t_std)
		if range_pool is None:
			range_pool = get_range_pool(t_mean, t_std, model.uncertainties, model.uncertainty_base, model.curve, max_cpus=max_cpus, max_queue_size=max_queue_size)
		if not range_pool.size:
			raise Exception("Could not generate random dates")
		model._set_range_pool(range_pool, t_mean, t_std)
	
	sums = []
	sums_prev = None
	c = 0
	todo = model.npass
	with tqdm(total=todo * 2) as pbar:
		pbar.set_description("Convergence: %0.3f" % (c))
		while True:
			n_dists = max(4, (todo - len(sums)) + 1)
			if n_dists > 200:
				process_mp(test_distributions_worker, range(n_dists),
						[dates_n, t_mean, t_std, model.uncertainties, model.uncertainty_base, model.curve,
						model.uniform, range_pool],
						collect_fnc=test_distributions_collect, collect_args=[sums, pbar],
						max_cpus=max_cpus, max_queue_size=max_queue_size)
			else:
				for i in range(n_dists):
					sums.append(test_distributions_worker(i, dates_n, t_mean, t_std, model.uncertainties, model.uncertainty_base, model.curve, model.uniform, range_pool))
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

