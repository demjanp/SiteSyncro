
from sitesyncro.utils.fnc_radiocarbon import (get_curve, calibrate)
from sitesyncro.utils.fnc_stat import (calc_sum, calc_mean, calc_mean_std, calc_range, calc_range_approx)
from sitesyncro.utils.fnc_simulate import (get_params, gen_random_dists_uniform, gen_random_dists_normal)
from sitesyncro.utils.fnc_load import (load_data)

import os
import time
import numpy as np
from matplotlib import pyplot
from scipy.interpolate import interp1d
from typing import List, Dict

T_MEAN = 3500
T_STD = 200
DATES_N = 10
UNCERTAINTY = 20

def get_params_alt(distributions, curve, uniform):
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
		rng = calc_range(curve[:,0], sum_dist, p = 0.9545)
		mean_sum = np.mean(rng)
		std_sum = abs(np.diff(rng)[0]) / 2
	else:
		# Calculate the mean and standard deviation of the summed distribution
		mean_sum, std_sum = calc_mean_std(curve[:, 0], sum_dist)
	return mean_sum, std_sum

def gen_random_dists_uniform_alt(dates_n: int, t_mean: float, t_std: float, uncertainties: List[float], uncertainty_base: float, curve: np.ndarray) -> List[np.ndarray]:
	
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
		rng = calc_range_approx(curve[:,0], s_d, p = 0.9545)
#		rng = calc_range(curve[:,0], s_d, p = 0.6827)
		print("\r%d/%d      " % (len(distributions), dates_n), end = '')
		if (rng[0] <= t_mean + t_std) and (rng[1] >= t_mean - t_std):
#			rng = calc_range(curve[:,0], s_d, p = 0.6827)
#			if (rng[0] <= t_mean + t_std) and (rng[1] >= t_mean - t_std):
			if 1:
				distributions.append(dist)
				sum_dists += dist
	return distributions


def gen_random_dists_normal_alt(dates_n: int, t_mean: float, t_std: float, uncertainties: List[float],
								  uncertainty_base: float, curve: np.ndarray) -> List[np.ndarray]:
	
	distributions = []
	sum_dists = np.zeros(curve.shape[0], dtype = float)
	
	c_threshold = t_std * 0.001
	
	while len(distributions) < dates_n:
		cal_age = np.random.normal(t_mean, t_std)
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
		m, s = calc_mean_std(curve[:, 0], s_d)
		c = max(abs(m - t_mean), abs(s - t_std))
		print("\r%d/%d %.2f     " % (len(distributions), dates_n, c), end = '')
		if c <= c_threshold:
			distributions.append(dist)
			sum_dists += dist
	return distributions


def sim_event(cal_age_bp, curve):
	# Simulate a C-14 dated event
	
	# Find the closest index in the calibration curve
	idx = np.argmin(np.abs(curve[:, 0] - cal_age_bp))
	# Get the corresponding radiocarbon age
	return curve[idx, 1]

def sim_uniform(cal_bp_mean, gap, dates_n, uncertainty_base, curve):
	# Generate uniformly distributed events around cal_bp_mean +-gap
	means = np.random.uniform(cal_bp_mean - gap, cal_bp_mean + gap, dates_n)
	means.sort()
	dates = []
	for cal_age in means:
		c14age = sim_event(cal_age, curve)
		uncert = uncertainty_base * np.exp(c14age / (2 * 8033))
		dates.append([c14age, np.random.normal(uncert, 2)])
	return dates

def sim_normal(cal_bp_mean, gap, dates_n, uncertainty_base, curve):
	# Generate uniformly distributed events around cal_bp_mean +-gap
	means = np.random.normal(cal_bp_mean, gap, dates_n)
	means.sort()
	dates = []
	for cal_age in means:
		c14age = sim_event(cal_age, curve)
		uncert = uncertainty_base * np.exp(c14age / (2 * 8033))
		dates.append([c14age, np.random.normal(uncert, 2)])
	return dates

if __name__ == '__main__':
	
	curve = get_curve()
	
	years = curve[:,0]
	
	'''
	samples, contexts, context_area, long_lived, redeposited, outlier, r_dates, context_phase, earlier_than = load_data('data_20240329.csv')
	dists_obs = []
	uncertainties = []
	for name in samples:
		age, uncert = r_dates[name]
		dists_obs.append(calibrate(age, uncert, curve))
		uncertainties.append(uncert)
	'''
	
	np.random.seed(1978)
	
	dates = sim_uniform(T_MEAN, T_STD, DATES_N, UNCERTAINTY, curve)
#	dates = sim_normal(T_MEAN, T_STD, DATES_N, UNCERTAINTY, curve)
	dists_obs = []
	uncertainties = []
	for age, uncert in dates:
		dists_obs.append(calibrate(age, uncert, curve))
		uncertainties.append(uncert)
	sum_obs = calc_sum(dists_obs)
#	mean_sum, std_sum = get_params_alt(dists_obs, curve, True)
	mean_sum, std_sum = get_params(dists_obs, curve, True)
#	mean_sum, std_sum = get_params(dists_obs, curve, False)
	
	dates_n = len(dists_obs)
	
	sums = []
	ts = []
	n_iters = 5
	for n in range(n_iters):
		print("\n\n%d/%d" % (n+1, n_iters))
		t0 = time.time()
#		dists_rnd = gen_random_dists_uniform_alt(dates_n, mean_sum, std_sum, uncertainties, UNCERTAINTY, curve)
		dists_rnd = gen_random_dists_uniform(dates_n, mean_sum, std_sum, uncertainties, UNCERTAINTY, curve)
#		dists_rnd = gen_random_dists_normal_alt(dates_n, mean_sum, std_sum, uncertainties, UNCERTAINTY, curve)
#		dists_rnd = gen_random_dists_normal(dates_n, mean_sum, std_sum, uncertainties, UNCERTAINTY, curve)
		ts.append(time.time() - t0)
		sum_rnd = calc_sum(dists_rnd)
		rng_rnd = calc_range(years, sum_rnd)
		sums.append([sum_rnd, rng_rnd])
	
	rng_obs = calc_range(years, sum_obs)
	
	pyplot.plot(years, sum_obs, color = "k")
	for sum_rnd, rng_rnd in sums:
		pyplot.plot(years, sum_rnd, color = "green", alpha = 0.3)
	pyplot.xlim(mean_sum + 4*std_sum, mean_sum - 4*std_sum)
	pyplot.show()
	
	print()
	print()
	print("time:", np.mean(ts))
	print()
	print(np.round(rng_obs, 2))
	print("generated")
	for sum_rnd, rng_rnd in sums:
		print("\t", np.round(rng_rnd, 2))
	print()
	
