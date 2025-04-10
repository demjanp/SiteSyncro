import numpy as np
from itertools import product
import itertools
from functools import reduce
import operator
import time
import copy

from sitesyncro.utils.fnc_mp import (process_mp)

def check_solution(phasing, phasing_ranges):
	
	for i, ph in enumerate(phasing):
		if ph < phasing_ranges[i][0] or ph > phasing_ranges[i][1]:
			return False
	return True

def calc_rng(phasing, ranges_gr):
	"""
	Calculate the range of the grouped ranges for the given phasing values.
	Args:
		phasing: List of integers. The phasing values for each sample.
		ranges_gr: List of lists of floats. Each list contains the maximum and minimum values of the grouped ranges for each sample.

	Returns:
		Float. The range of the grouped ranges.
	"""
	
	ranges_found = dict([(ph, [-np.inf, np.inf]) for ph in set(phasing) if ph > -1])
	for i, ph in enumerate(phasing):
		if ph == -1:
			continue
		ranges_found[ph][0] = max(ranges_found[ph][0], ranges_gr[i][0])
		ranges_found[ph][1] = min(ranges_found[ph][1], ranges_gr[i][1])
	rng = sum((ranges_found[ph][0] - ranges_found[ph][1]) for ph in ranges_found) / len(ranges_found)
	return rng

def calc_dice(phasing, ranges_gr):
	
	def dice_coefficient(r1, r2, s1, s2):
		L_r = r2 - r1
		L_s = s2 - s1
		intersection_length = max(0, min(r2, s2) - max(r1, s1))
		dice_coef = (2 * intersection_length) / (L_r + L_s)
		return dice_coef
	
	ranges_found = dict([(ph, [-np.inf, np.inf]) for ph in set(phasing) if ph > -1])
	for i, ph in enumerate(phasing):
		if ph == -1:
			continue
		ranges_found[ph][0] = max(ranges_found[ph][0], ranges_gr[i][0])
		ranges_found[ph][1] = min(ranges_found[ph][1], ranges_gr[i][1])
	d_sum = 0
	for i, ph in enumerate(phasing):
		if ph == -1:
			continue
		r2, r1 = ranges_gr[i]
		s2, s1 = ranges_found[ph]
		d_sum += dice_coefficient(r1, r2, s1, s2)
	return d_sum / len(phasing)

def reduce_phasing_permutation(phasing_ranges, ranges_gr):
	"""
	Optimize phasing ranges to minimize the range of the grouped ranges.
	Args:
		phasing_ranges: List of lists of integers. Each list contains the minimum and maximum phasing values for each sample.
		ranges_gr: List of lists of floats. Each list contains the maximum and minimum values of the grouped ranges for each sample.

	Returns:
		List of integers. The optimized phasing values for each sample.
	"""
	
	phasing_opt = None
	d_opt = -np.inf
	values = [list(range(phase_min, phase_max+1)) for phase_min, phase_max in phasing_ranges]
	cmax = reduce(operator.mul, (len(value) for value in values), 1)
	cnt = 1
	for phasing in product(*values):
		if cnt%1000 == 0:
			print("\rOptimizing phasing: %d/%d, d:%.3f" % (cnt, cmax, d_opt), end="")
		cnt += 1
		d = calc_dice(phasing, ranges_gr)
		if d > d_opt:
			d_opt = d
			phasing_opt = phasing
	print()
	phasing = list(phasing_opt)
	return phasing

def reduce_phasing_permutation_alt(phasing_ranges, ranges_gr):
	phase_max = max(rng[1] for rng in phasing_ranges)
	idxs_all = list(range(len(phasing_ranges)))
	idxs_moving = [i for i in idxs_all if phasing_ranges[i][0] != phasing_ranges[i][1]]
	phasing_opt = [-1]*len(phasing_ranges)
	for idx in idxs_all:
		if idx not in idxs_moving:
			phasing_opt[idx] = phasing_ranges[idx][0]
	while -1 in phasing_opt:
		ds0 = []
		res0 = []
		for i0 in idxs_moving:
			phasing0 = np.array(phasing_opt, dtype = int)
			for ph0 in range(phasing_ranges[i0][0], phasing_ranges[i0][1]+1):
				phasing0[i0] = ph0
				while (phasing0 == -1).any():
					ds1 = []
					res1 = []
					for i1 in np.where(phasing0 == -1)[0]:
						phasing1 = phasing0.copy()
						for ph1 in range(phasing_ranges[i1][0], phasing_ranges[i1][1]+1):
							phasing1[i1] = ph1
							ds1.append(calc_dice(phasing1, ranges_gr))
							res1.append([i1, ph1])
					i1, ph1 = res1[np.argmax(ds1)]
					phasing0[i1] = ph1
				ds0.append(calc_dice(phasing0, ranges_gr))
				res0.append(phasing0.tolist())
		if ds0:
			phasing_opt = res0[np.argmax(ds0)]
	d_opt = calc_dice(phasing_opt, ranges_gr)
	changed = True
	while changed:
		changed = False
		for idx in idxs_moving:
			phs = list(range(phasing_ranges[idx][0], phasing_ranges[idx][1]+1))
			phasing = copy.copy(phasing_opt)
			ds = []
			for ph in phs:
				phasing[idx] = ph
				ds.append(calc_dice(phasing, ranges_gr))
			d_max = max(ds)
			if d_max > d_opt:
				phasing_opt[idx] = phs[np.argmax(ds)]
				d_opt = d_max
				changed = True
	
	phasing_multi = {}
	for idx in idxs_moving:
		ph_collect = [phasing_opt[idx]]
		phasing = copy.copy(phasing_opt)
		for ph in range(phasing_ranges[idx][0], phasing_ranges[idx][1]+1):
			if ph in ph_collect:
				continue
			phasing[idx] = ph
			d = calc_dice(phasing, ranges_gr)
			if d>=d_opt:
				ph_collect.append(ph)
		if len(ph_collect) > 1:
			phasing_multi[idx] = ph_collect
	
	if phasing_multi:
		print("PHS max:", max(len(phs) for phs in phasing_multi.values()))
	
	return phasing_opt


if __name__ == '__main__':
	
#	np.random.seed(1978)
	
	n_samples = 20
	ph_max = 8
	range_min = 2500
	range_max = 3500
	
	cnt_ok = 0
	cnt_fail = 0
	for r in range(1000):
	
		phasing_ranges = []
		for i in range(n_samples):
			phase_min = np.random.randint(0, ph_max+1)
			phase_max = np.random.randint(phase_min, ph_max+1)
			phasing_ranges.append([phase_min, phase_max])
		
		ranges_gr = []
		for i in range(n_samples):
			rng_min = np.random.uniform(range_min, range_max)
			rng_max = np.random.uniform(rng_min, range_max)
			ranges_gr.append([rng_max, rng_min])
		
		phasing1 = reduce_phasing_permutation_alt(phasing_ranges, ranges_gr)
		'''
		phasing2 = reduce_phasing_permutation(phasing_ranges, ranges_gr)
		d1, sol_ok = calc_dice(phasing1, ranges_gr), check_solution(phasing1, phasing_ranges)
		d2 = calc_dice(phasing2, ranges_gr)
		print()
		print("iter:", r+1, "ok:", cnt_ok, "fail:", cnt_fail)
		if (not sol_ok) or (d1 != d2):
			cnt_fail += 1
			print()
			print(phasing1)
			print(phasing2)
			print("Sol ok:", sol_ok)
			print("d1:", d1)
			print("d2:", d2)
			print()
		else:
			cnt_ok += 1
		'''
	
	'''
	t0 = time.time()
	phasing1 = reduce_phasing_permutation_alt(phasing_ranges, ranges_gr)
	t1 = time.time() - t0
	
	t0 = time.time()
	phasing2 = reduce_phasing_permutation(phasing_ranges, ranges_gr)
	t2 = time.time() - t0
	
	print(phasing1, calc_rng(phasing1, ranges_gr), check_solution(phasing1, phasing_ranges), t1)
	print(phasing2, calc_rng(phasing2, ranges_gr), check_solution(phasing2, phasing_ranges), t2)
	'''
	