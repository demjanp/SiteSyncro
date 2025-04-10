import numpy as np

from sitesyncro import Model
from sitesyncro.utils.fnc_phase import (prob_earlier_than)
from sitesyncro.utils.fnc_stat import (calc_range)

def prob_earlier_than_alt(dist1: np.ndarray, dist2: np.ndarray, n_iters: int = 1000) -> float:
   
   years = np.arange(dist1.shape[0])
   years1 = np.random.choice(years, n_iters, p=dist1)
   years2 = np.random.choice(years, n_iters, p=dist2)
   p = (years1 > years2).sum() / n_iters
   return p

P_VALUE = 0.01

if __name__ == '__main__':
	
	model = Model(directory = "stage_1")
	samples = list(model.samples.keys())
	
	s1 = [s for s in model.samples if s.startswith('108.522.51.12')][0]
	s2 = [s for s in model.samples if s.startswith('109.523.586.5')][0]
	
	rng1 = calc_range(model.years, model.samples[s1].likelihood)
	rng2 = calc_range(model.years, model.samples[s2].likelihood)
	
	earlier_than, samples = model.mphasing.create_earlier_than_matrix()
	i1, i2 = samples.index(s1), samples.index(s2)
	earlier_than2 = model.mphasing.update_earlier_than_by_dating(earlier_than, samples)
	
	print()
	print(rng1)
	print(rng2)
	print(earlier_than[i1,i2])
	print(earlier_than2[i1,i2])
	print()
	raise SystemExit()
	
	
	ranges = {}
	for s1 in samples:
		ranges[s1] = calc_range(model.years, model.samples[s1].likelihood, p = 1-P_VALUE)
		ranges[s1] = [round(ranges[s1][0]), round(ranges[s1][1])]
	
	cnt1 = 0
#	cnt2 = 0
	cnt_all = 0
	
	for s1 in samples:
		dist1 = model.samples[s1].likelihood
#		rng1 = model.samples[s1].likelihood_range
		rng1 = ranges[s1]
		for s2 in samples:
			if s1 == s2:
				continue
			dist2 = model.samples[s2].likelihood
#			rng2 = model.samples[s2].likelihood_range
			rng2 = ranges[s2]
			p1 = prob_earlier_than(dist1, dist2)
#			p2 = prob_earlier_than_alt(dist1, dist2)
			is_earlier_1 = (p1 >= (1 - P_VALUE))
#			is_earlier_2 = (p2 >= (1 - P_VALUE))
			is_earlier_3 = (rng1[1] > rng2[0])
			cnt_all += 1
			if (is_earlier_1 != is_earlier_3):
				cnt1 += 1
#			if (is_earlier_2 != is_earlier_3):
#				cnt2 += 1
			
			if (is_earlier_1 != is_earlier_3):
				print()
				print("prob1:", is_earlier_1, round(p1,3),">",round(1 - P_VALUE, 3))
#				print("prob2:", is_earlier_2, round(p2,3),">",round(1 - P_VALUE, 3))
				print("rng:", is_earlier_3, rng1,">",rng2, rng1[1] - rng2[0])
			
	print("prob_earlier_than:", cnt1/cnt_all)
#	print("prob_earlier_than_alt:", cnt2/cnt_all)
	