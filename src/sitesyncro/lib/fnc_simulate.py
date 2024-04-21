from sitesyncro.lib.fnc_radiocarbon import (calibrate)
from sitesyncro.lib.fnc_sum import (sum_distributions, calc_range, calc_mean_std, calc_percentiles)
from sitesyncro.lib.fnc_mp import (process_mp)

import numpy as np
from scipy.stats import norm

def generate_random_distributions(dates_n, t_param1, t_param2, uncert_base, curve, uniform):
	# Generate sets of randomized distributions based on observed distributions
	#
	# dates_n: number of dates to generate
	# t_param1: mean or median of the calendar ages
	# t_param2: standard deviation or range of the calendar ages
	# uncert_base: base uncertainty to simulate radiocarbon dates
	# curve: [[CalBP, ConvBP, CalSigma], ...]
	# uniform: flag indicating whether to use a uniform distribution for the calendar ages
	
	dates = []
	distributions = []
	for i in range(dates_n):
		# Generate a random calendar age
		if uniform:
			t = np.random.uniform(t_param1 - t_param2, t_param1 + t_param2)
		else:
			t = np.random.normal(t_param1, t_param2)
		# Find the closest index in the calibration curve
		idx = np.argmin(np.abs(curve[:, 0] - t))
		# Get the corresponding radiocarbon age
		age = curve[idx, 1]
		# Calculate the standard deviation for the date
		stdev = uncert_base * np.exp(age / (2 * 8033))
		# Append the date and its standard deviation to the dates list
		dates.append([age, stdev])
		# Calibrate the date and append the distribution to the distributions list
		distribution = calibrate(age, stdev, curve, normalize=True)
		distributions.append(distribution)
	
	# Sum the distributions of the generated dates
	dist_summed = sum_distributions(distributions)
	
	if uniform:
		# Calculate weighted median and range
		median_sum = np.average(curve[:, 0], weights=dist_summed)
		range_sum = np.sqrt(np.average((curve[:, 0] - median_sum) ** 2, weights=dist_summed))
	else:
		# Calculate the mean and standard deviation of the summed distribution
		mean_sum, std_sum = calc_mean_std([curve[:, 0], dist_summed])
	
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
		if uniform:
			dates[i][0] += 1 * (t_param1 - median_sum)
			dates[i][0] *= 1 + scaling_factor * ((t_param2 / range_sum) - 1)
		else:
			dates[i][0] += (t_param1 - mean_sum)
			dates[i][0] *= 1 + scaling_factor * ((t_param2 / std_sum) - 1)
		distributions[i] = calibrate(dates[i][0], dates[i][1], curve, normalize=True)
		
		# Recalculate the sum of the distributions and their mean or median
		dist_summed = sum_distributions(distributions)
		if uniform:
			median_sum = np.average(curve[:, 0], weights=dist_summed)
			range_sum = np.sqrt(np.average((curve[:, 0] - median_sum) ** 2, weights=dist_summed))
		else:
			mean_sum, std_sum = calc_mean_std([curve[:, 0], dist_summed])
		# Decrease the scaling factor
		scaling_factor *= 0.999
	
	return distributions

def worker_fnc(params, dates_n, t_param1, t_param2, uncert_base, curve, uniform):
	
	return generate_random_distributions(dates_n, t_param1, t_param2, uncert_base, curve, uniform)

def collect_fnc(data, results):
	
	# data = distributions
	# distributions = [[p, ...], ...]
	
	results.append(data)

def progress_fnc(done, todo, all_done, all_todo, c):
	
	print("\rIteration: %d/%d, Conv: %0.4f         " % (all_done + done, all_todo, c), end = '')

def calculate_parameters(years, distribution, uniform):
	if uniform:
		# Calculate weighted median and range of the summed distribution
		t_median = np.average(years, weights=distribution)
		t_range = np.sqrt(np.average((years - t_median) ** 2, weights=distribution))
		t_param1, t_param2 = t_median, t_range
	else:
		# Calculate the mean and standard deviation of the summed distribution
		t_mean, t_std = calc_mean_std([years, distribution])
		t_param1, t_param2 = t_mean, t_std

	return t_param1, t_param2

def test_distributions(distributions, curve, p_value = 0.05, uncert_base = 15, uniform = False,
                   npass = 100, convergence = 0.99, max_cpus = -1, max_queue_size = 100):
	# Randomization test of the null hypothesis that the observed distributions
	# represent a normal / uniform distribution of events
	#
	# distributions = {name: [p, ...], ...}
	# curve: [[CalBP, ConvBP, CalSigma], ...]
	# uncert_base: base uncertainty to simulate radiocarbon dates
	# uniform: flag indicating whether to use a uniform distribution for the calendar ages
	#
	# returns: [[age, stdev], ...], [distributions, ...]
	
	# Number of dates to generate
	dates_n = len(distributions)

	# Sum the distributions of all samples
	sum_obs = sum_distributions(list(distributions.values()))
	
	years = curve[:, 0]
	
	# Calculate the mean or median and standard deviation or range of the summed distribution
	t_param1, t_param2 = calculate_parameters(years, sum_obs, uniform)
	
	# Generate sets of randomized distributions and their sums
	distributions_rnd = []
	sums = []
	sums_prev = None
	c = 0
	params_list = list(range(npass))
	todo = npass
	while True:
		process_mp(worker_fnc, params_list, [dates_n, t_param1, t_param2, uncert_base, curve, uniform],
		           collect_fnc = collect_fnc, collect_args = [distributions_rnd],
		           progress_fnc = progress_fnc, progress_args = [len(distributions_rnd), npass*2, c],
		           max_cpus = max_cpus, max_queue_size = max_queue_size)
		if len(distributions_rnd) >= todo:
			sums += [sum_distributions(dists) for dists in distributions_rnd[-(len(distributions_rnd) - len(sums)):]]
			sums_m = np.array(sums).mean(axis = 0)
			if sums_prev is not None:
				c = ((sums_prev * sums_m).sum()**2) / ((sums_prev**2).sum() * (sums_m**2).sum())
			sums_prev = sums_m.copy()
			if c >= convergence:
				print("\nConverged at:", c)
				break
			todo *= 2
	
	sums_rnd = np.array(sums)
	
	# Calculate p-value
	mask = (sums_rnd.std(axis = 0) > 0)
	p = (1 - norm(sums_rnd[:,mask].mean(axis = 0), sums_rnd[:,mask].std(axis = 0)).cdf(sum_obs[mask])).min()
	
	# Find the range of years with non-zero probabilities
	threshold = 0.001 * min(sum_obs.max(), sums_rnd.max())
	j1 = np.argwhere(sum_obs > threshold).min()
	j2 = np.argwhere(sum_obs > threshold).max()
	for i in range(sums_rnd.shape[1]):
		if sums_rnd[:,i].max() > threshold:
			j1 = min(j1, i)
			j2 = max(j2, i)
	for distrs in distributions_rnd:
		for i in range(len(distrs)):
			if distrs[i].max() > threshold:
				j1 = min(j1, i)
				j2 = max(j2, i)
	
	# Trim the years, observed distribution, and randomized distributions
	years = years[j1:j2+1]
	sum_obs = sum_obs[j1:j2+1]
	sums_rnd = sums_rnd[:,j1:j2+1]
	for i in range(len(distributions_rnd)):
		distributions_rnd[i] = [distr[j1:j2+1] for distr in distributions_rnd[i]]
	
	perc_lower = (p_value * 100) / 2
	perc_upper = 100 - perc_lower
	sums_rnd_lower, sums_rnd_upper = calc_percentiles(sums_rnd, perc_lower, perc_upper)
	
	return years, sum_obs, distributions_rnd, sums_rnd_lower, sums_rnd_upper, p

