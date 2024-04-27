from sitesyncro.utils.fnc_radiocarbon import (calibrate)
from sitesyncro.utils.fnc_stat import (calc_sum, calc_mean_std, calc_percentiles, samples_to_distributions)
from sitesyncro.utils.fnc_mp import (process_mp)

import numpy as np
from tqdm import tqdm
from scipy.stats import norm

def generate_random_distributions(dates_n, t_param1, t_param2, uncertainties, uncertainty_base, curve, uniform):
	# Generate sets of randomized distributions based on observed distributions
	#
	# dates_n: number of dates to generate
	# t_param1: mean or median of the calendar ages
	# t_param2: standard deviation or range of the calendar ages
	# uncertainties: [uncertainty, ...] pool of uncertainties to simulate C-14 dates
	# uncertainty_base: base uncertainty to calculate uncertainty based on C-14 age if pool is not available
	# curve: [[CalBP, ConvBP, CalSigma], ...]
	# uniform: flag indicating whether to use a uniform distribution for the calendar ages
	
	def _sim_age():
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
	
	dates = []
	distributions = []
	for i in range(dates_n):
		# Generate a random calendar age
		age = _sim_age()
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
	
	# Sum the distributions of the generated dates
	dist_summed = calc_sum(distributions)
	
	if uniform:
		# Calculate weighted median and range
		median_sum = np.average(curve[:, 0], weights=dist_summed)
		range_sum = np.sqrt(np.average((curve[:, 0] - median_sum) ** 2, weights=dist_summed))
	else:
		# Calculate the mean and standard deviation of the summed distribution
		mean_sum, std_sum = calc_mean_std(curve[:, 0], dist_summed)
	
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
			age_new = _sim_age()
		dates[i][0] = age_new
		distributions[i] = calibrate(dates[i][0], dates[i][1], curve)
		
		# Recalculate the sum of the distributions and their mean or median
		dist_summed = calc_sum(distributions)
		if uniform:
			median_sum = np.average(curve[:, 0], weights=dist_summed)
			range_sum = np.sqrt(np.average((curve[:, 0] - median_sum) ** 2, weights=dist_summed))
		else:
			mean_sum, std_sum = calc_mean_std(curve[:, 0], dist_summed)
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
	
	return generate_random_distributions(dates_n, t_param1, t_param2, uncertainties, uncertainty_base, curve, uniform)

def collect_fnc(data, results, pbar):
	
	# data = distributions
	# distributions = [[p, ...], ...]
	
	pbar.update(1)
	results.append(data)

def progress_fnc(done, todo, all_done, all_todo, c, pbar):
	
	print("\rIteration: %d/%d, Conv: %0.4f         " % (all_done + done, all_todo, c), end = '')

def test_distributions(model, max_cpus = -1, max_queue_size = 100):
	
	distributions, _, _ = samples_to_distributions(model.samples.values())
	
	dates_n = len(distributions)
	
	if dates_n < 3:
		raise Exception("Insufficient number samples for testing.")
	
	print("Testing %d distributions" % dates_n)
	
	sum_obs = calc_sum(distributions)
	years = model.curve[:, 0]
	t_param1, t_param2 = calculate_parameters(years, sum_obs, model.uniform)
	
	distributions_rnd = []
	sums = []
	sums_prev = None
	c = 0
	todo = model.npass
	params_list = list(range(model.npass))
	with tqdm(total=todo) as pbar:
		pbar.total = model.npass*2
		pbar.set_description("Convergence: %0.3f" % (c))
		while True:
			process_mp(worker_fnc, params_list, [dates_n, t_param1, t_param2, model.uncertainties, model.uncertainty_base, model.curve, model.uniform],
		           collect_fnc = collect_fnc, collect_args = [distributions_rnd, pbar],
		           max_cpus = max_cpus, max_queue_size = max_queue_size)
			if len(distributions_rnd) >= todo:
				sums += [calc_sum(dists) for dists in distributions_rnd[-(len(distributions_rnd) - len(sums)):]]
				sums_m = np.array(sums).mean(axis = 0)
				if sums_prev is not None:
					c = ((sums_prev * sums_m).sum()**2) / ((sums_prev**2).sum() * (sums_m**2).sum())
					pbar.set_description("Convergence: %0.3f" % (c))
				sums_prev = sums_m.copy()
				if c >= model.convergence:
					break
				todo *= 2
				pbar.total = max(todo, model.npass*2)
	
	sums_rnd = np.array(sums)
	
	mask = (sums_rnd.std(axis = 0) > 0)
	p = (1 - norm(sums_rnd[:,mask].mean(axis = 0), sums_rnd[:,mask].std(axis = 0)).cdf(sum_obs[mask])).min()
	
	perc_lower = (model.p_value * 100) / 2
	perc_upper = 100 - perc_lower
	sums_rnd_lower, sums_rnd_upper = calc_percentiles(sums_rnd, perc_lower, perc_upper)
	
	return sum_obs, sums_rnd_lower, sums_rnd_upper, p
