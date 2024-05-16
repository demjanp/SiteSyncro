from typing import Any, List

import numpy as np
from scipy.stats import norm
from tqdm import tqdm

from sitesyncro.utils.fnc_simulate import (get_params, generate_random_distributions)
from sitesyncro.utils.fnc_mp import (process_mp)
from sitesyncro.utils.fnc_stat import (calc_sum, calc_percentiles, samples_to_distributions)

def test_distributions_worker(params: Any, dates_n: int, t_mean: float, t_std: float, uncertainties: List[float],
			   uncertainty_base: float, curve: np.ndarray, uniform: bool) -> np.ndarray:
	
	return calc_sum(generate_random_distributions(dates_n, t_mean, t_std, uncertainties, uncertainty_base, curve, uniform))


def test_distributions_collect(data: np.ndarray, results: List[np.ndarray], pbar: tqdm) -> None:
	# data = dist_sum
	# dist_sum = [p, ...]
	pbar.update(1)
	results.append(data)


class MRandomization(object):
	"""
	A class implementing randomization testing functionality
	
	:param model: The parent Model object
	:type model: Model
	"""
	def __init__(self, model: 'Model'):
		
		self.model = model
	
	def test_distributions(self, max_cpus: int = -1, max_queue_size: int = -1) -> (
			np.ndarray, np.ndarray, np.ndarray, float):
		"""
		Test if the samples in the model are distributed uniformly / normally in time
		
		Parameters:
		max_cpus (int): Maximum number of CPUs to use for multiprocessing. Default is -1 (use all available CPUs).
		max_queue_size (int): Maximum size of the queue for multiprocessing. Default is -1 (no limit).
		
		Returns:
		(sum_obs, sums_rnd_lower, sums_rnd_upper, p)
			- sum_obs: Sum of the observed probability distributions of the sample dates.
			- sums_rnd_lower: Lower bound of the randomized sum of the probability distributions.
			- sums_rnd_upper: Upper bound of the randomized sum of the probability distributions.
			- p: P-value for the test.
		"""
		
		distributions, _, _ = samples_to_distributions(self.model.samples.values())
		
		dates_n = len(distributions)
		
		if dates_n < 3:
			raise Exception("Insufficient number samples for testing.")
		
		print("Testing %d distributions" % dates_n)
		
		t_mean, t_std = get_params(distributions, self.model.curve, self.model.uniform)
		
		sums = []
		sums_prev = None
		c = 0
		todo = self.model.npass
		with tqdm(total=todo * 2) as pbar:
			pbar.set_description("Convergence: %0.3f" % (c))
			while True:
				n_dists = max(4, (todo - len(sums)) + 1)
				if n_dists > 50:
					process_mp(test_distributions_worker, range(n_dists),
							[dates_n, t_mean, t_std, self.model.uncertainties, self.model.uncertainty_base, self.model.curve, self.model.uniform],
							collect_fnc=test_distributions_collect, collect_args=[sums, pbar],
							max_cpus=max_cpus, max_queue_size=max_queue_size)
				else:
					for i in range(n_dists):
						sums.append(test_distributions_worker(i, dates_n, t_mean, t_std, self.model.uncertainties, self.model.uncertainty_base, self.model.curve, self.model.uniform))
						pbar.update(1)
				
				if len(sums) >= todo:
					sums_m = np.array(sums).mean(axis=0)
					if sums_prev is not None:
						c = ((sums_prev * sums_m).sum() ** 2) / ((sums_prev ** 2).sum() * (sums_m ** 2).sum())
						pbar.set_description("Convergence: %0.3f" % (c))
					sums_prev = sums_m.copy()
					if c >= self.model.convergence:
						break
					todo *= 2
					pbar.total = max(todo, self.model.npass * 2)
					pbar.refresh()
		
		sums_rnd = np.array(sums)
		
		sum_obs = calc_sum(distributions)
		mask = (sums_rnd.std(axis=0) > 0)
		p = (1 - norm(sums_rnd[:, mask].mean(axis=0), sums_rnd[:, mask].std(axis=0)).cdf(sum_obs[mask])).min()
		
		perc_lower = (self.model.p_value * 100) / 2
		perc_upper = 100 - perc_lower
		sums_rnd_lower, sums_rnd_upper = calc_percentiles(sums_rnd, perc_lower, perc_upper)
		
		return sum_obs, sums_rnd_lower, sums_rnd_upper, p
