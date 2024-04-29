from sitesyncro.utils.fnc_radiocarbon import (calibrate)
from sitesyncro.utils.fnc_stat import (calc_sum, calc_mean_std, calc_percentiles, samples_to_distributions)
from sitesyncro.utils.fnc_mp import (process_mp)

import numpy as np
from tqdm import tqdm
from scipy.stats import norm

def generate_random_distributions(dates_n, t_param1, t_param2, uncertainties, uncertainty_base, curve, uniform, max_iterations = 10000):
	# Generate sets of randomized distributions based on observed distributions
	#
	# dates_n: number of dates to generate
	# t_param1: mean or median of the calendar ages
	# t_param2: standard deviation or range of the calendar ages
	# uncertainties: [uncertainty, ...] pool of uncertainties to simulate C-14 dates
	# uncertainty_base: base uncertainty to calculate uncertainty based on C-14 age if pool is not available
	# curve: [[CalBP, ConvBP, CalSigma], ...]
	# uniform: flag indicating whether to use a uniform distribution for the calendar ages
	
	def _sim_age(t_param1, t_param2, uniform, curve):
		# Generate a random calendar age
		if uniform:
			t = np.random.uniform(t_param1 - t_param2, t_param1 + t_param2)
		else:
			t = np.random.normal(t_param1, t_param2)
		# Find the closest index in the calibration curve
		idx = np.argmin(np.abs(curve[:, 0] - t))
		# Get the corresponding radiocarbon age
		age = curve[idx, 1]
		return age
	
	def _gen_distributions(dates_n, uncertainties, uncertainty_base, t_param1, t_param2, uniform, curve):
		dates = []
		distributions = []
		for i in range(dates_n):
			# Generate a random calendar age
			age = _sim_age(t_param1, t_param2, uniform, curve)
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
		return dates, distributions
	
	def _analyze(distributions, uniform, curve):
		# Sum the distributions
		dist_summed = calc_sum(distributions)
		median_sum, range_sum, mean_sum, std_sum = None, None, None, None
		if uniform:
			# Calculate weighted median and range
			median_sum = np.average(curve[:, 0], weights=dist_summed)
			range_sum = np.sqrt(np.average((curve[:, 0] - median_sum) ** 2, weights=dist_summed))
		else:
			# Calculate the mean and standard deviation of the summed distribution
			mean_sum, std_sum = calc_mean_std(curve[:, 0], dist_summed)
		return median_sum, range_sum, mean_sum, std_sum
	
	
	dates, distributions = _gen_distributions(dates_n, uncertainties, uncertainty_base, t_param1, t_param2, uniform, curve)
	median_sum, range_sum, mean_sum, std_sum = _analyze(distributions, uniform, curve)
	
	curve_min = curve[:,1].min()
	curve_max = curve[:,1].max()
	
	# Adjust the dates until the mean and standard deviation of the summed distribution are close to the desired values
	iteration = 0
	scaling_factor = 1.0
	if uniform:
		c_threshold = t_param2 * 0.01
	else:
		c_threshold = t_param2 * 0.01
	while True:
		if uniform:
			c = max(abs(median_sum - t_param1), abs(range_sum - t_param2))
		else:
			c = max(abs(mean_sum - t_param1), abs(std_sum - t_param2))
		
		if c <= c_threshold:
			break
		
		iteration += 1
		if iteration > max_iterations:
			# Re-initialize
			dates, distributions = _gen_distributions(dates_n, uncertainties, uncertainty_base, t_param1, t_param2, uniform, curve)
			median_sum, range_sum, mean_sum, std_sum = _analyze(distributions, uniform, curve)
			iteration = 0
			scaling_factor = 1.0
			continue
		
		# Randomly select a date to adjust
		i = np.random.choice(dates_n)
		age_new = dates[i][0]
		if uniform:
			age_new += 1 * (t_param1 - median_sum)
			age_new *= 1 + scaling_factor * ((t_param2 / range_sum) - 1)
		else:
			age_new += (t_param1 - mean_sum)
			age_new *= 1 + scaling_factor * ((t_param2 / std_sum) - 1)
		if (age_new < curve_min) or (age_new > curve_max):
			age_new = _sim_age(t_param1, t_param2, uniform, curve)
		dates[i][0] = age_new
		distributions[i] = calibrate(dates[i][0], dates[i][1], curve)
		# Recalculate the sum of the distributions and their mean or median
		median_sum, range_sum, mean_sum, std_sum = _analyze(distributions, uniform, curve)
		# Decrease the scaling factor
		scaling_factor *= 0.999
	
	return distributions

def calculate_parameters(years, distribution, uniform):
	if uniform:
		# Calculate weighted median and range of the summed distribution
		t_median = np.average(years, weights=distribution)
		t_range = np.sqrt(np.average((years - t_median) ** 2, weights=distribution))
		t_param1, t_param2 = t_median, t_range
	else:
		# Calculate the mean and standard deviation of the summed distribution
		t_mean, t_std = calc_mean_std(years, distribution)
		t_param1, t_param2 = t_mean, t_std

	return t_param1, t_param2

def worker_fnc(params, dates_n, t_param1, t_param2, uncertainties, uncertainty_base, curve, uniform):
	
	return calc_sum(generate_random_distributions(dates_n, t_param1, t_param2, uncertainties, uncertainty_base, curve, uniform))

def collect_fnc(data, results, pbar):
	
	# data = dist_sum
	# dist_sum = [p, ...]
	pbar.update(1)
	results.append(data)

def test_distributions(model, max_cpus = -1, max_queue_size = -1):
	
	distributions, _, _ = samples_to_distributions(model.samples.values())
	
	dates_n = len(distributions)
	
	if dates_n < 3:
		raise Exception("Insufficient number samples for testing.")
	
	print("Testing %d distributions" % dates_n)
	
	sum_obs = calc_sum(distributions)
	years = model.curve[:, 0]
	t_param1, t_param2 = calculate_parameters(years, sum_obs, model.uniform)
	
	sums = []
	sums_prev = None
	c = 0
	todo = model.npass
	with tqdm(total=todo*2) as pbar:
		pbar.set_description("Convergence: %0.3f" % (c))
		while True:
			process_mp(worker_fnc, range(max(4, (todo - len(sums)) + 1)), [dates_n, t_param1, t_param2, model.uncertainties, model.uncertainty_base, model.curve, model.uniform],
		           collect_fnc = collect_fnc, collect_args = [sums, pbar],
		           max_cpus = max_cpus, max_queue_size = max_queue_size)
			if len(sums) >= todo:
				sums_m = np.array(sums).mean(axis = 0)
				if sums_prev is not None:
					c = ((sums_prev * sums_m).sum()**2) / ((sums_prev**2).sum() * (sums_m**2).sum())
					pbar.set_description("Convergence: %0.3f" % (c))
				sums_prev = sums_m.copy()
				if c >= model.convergence:
					break
				todo *= 2
				pbar.total = max(todo, model.npass*2)
				pbar.refresh()
	
	sums_rnd = np.array(sums)
	
	mask = (sums_rnd.std(axis = 0) > 0)
	p = (1 - norm(sums_rnd[:,mask].mean(axis = 0), sums_rnd[:,mask].std(axis = 0)).cdf(sum_obs[mask])).min()
	
	perc_lower = (model.p_value * 100) / 2
	perc_upper = 100 - perc_lower
	sums_rnd_lower, sums_rnd_upper = calc_percentiles(sums_rnd, perc_lower, perc_upper)
	
	return sum_obs, sums_rnd_lower, sums_rnd_upper, p
