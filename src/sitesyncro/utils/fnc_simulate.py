from typing import List

import numpy as np

from sitesyncro.utils.fnc_radiocarbon import (calibrate)
from sitesyncro.utils.fnc_stat import (calc_sum, calc_mean_std, calc_range, calc_range_approx)


def get_params(distributions, curve, uniform):
	"""
	Calculate parameters of the summed distributions
	
	Returns:
	(mean_sum, std_sum)
	mean_sum: Mean (normal model) or the center of the 2-sigma range (uniform model)
	std_sum: Standard deviation (normal model) or 1/2 of the 2-sigma range (uniform model)
	"""
	# Sum the distributions
	sum_dist = calc_sum(distributions)
	if uniform:
		rng = calc_range(curve[:,0], sum_dist, p = 0.9545)
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
		cal_age = np.random.uniform(t_mean - t_std*2.5, t_mean + t_std*2.5)
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
		rng = calc_range_approx(curve[:,0], s_d, p = 0.9545)
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

