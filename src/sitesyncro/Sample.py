from sitesyncro.utils.fnc_stat import (calc_mean, calc_range)
from sitesyncro.utils.fnc_oxcal import (oxcal_date)
from sitesyncro.utils.fnc_radiocarbon import (calibrate)
from sitesyncro.utils.fnc_data import (dict_keys_to_int, dict_np_to_list)

import copy
import numpy as np

class Sample(object):
	
	def __init__(self, *args, **kwargs):
		# positional arguments:
		#	name: unique name (ID) of the sample
		# 	age: C-14 age (years BP) for date_type 'R'; mean calendar age (years BP) for date_type 'U'
		# 	uncertainty: uncertainty (years BP) for date_type 'R'; 1/2 range (years BP) for date_type 'U'
		# keyword arguments:
		# 	date_type: 'R' for radiocarbon date; 'U' for calendar date as a uniform distribution
		# 	long_lived: True if sample is older than the examined deposition event due to e.g. old wood effect or redeposition from older strata
		# 	context: name of the context where sample was found
		# 	area: excavation area
		# 	area_excavation_phase: chronological phase of the context within the excavation area (integer, lower = earlier (older) phase)
		# 	earlier_than: list of names of samples which are stratigraphically later (younger) than this sample
		# 	curve: [[calendar year BP, C-14 year, uncertainty], ...]; Radiocarbon calibration curve
		#
		# Alternatively, if only one argument of type dictionary is provided, the object is initialized using this data
		
		def _from_arguments(name, age, uncertainty,
				date_type = 'R',
				long_lived = False,
				context = None,
				area = None,
				area_excavation_phase = None,
				earlier_than = [],
				curve = None,
			):
			
			if date_type not in ['R','U']:
				raise Exception("Invalid date type specified: %s (must be 'R' or 'U')" % (date_type))
			
			self._data = dict(
				name = name,
				age = age,
				uncertainty = uncertainty,
				date_type = date_type,
				long_lived = long_lived,
				context = context,
				area = area,
				area_excavation_phase = area_excavation_phase,
				earlier_than = earlier_than,
				
				group = None,
				phase = None,
				years = None,
				likelihood = None,
				posterior = None,
				likelihood_range = None,
				likelihood_mean = None,
				posterior_range = None,
				posterior_mean = None,
				posterior_agreement = None,
			)
			
			if curve is not None:
				self.calibrate(curve)
		
		if len(args) == 1 and isinstance(args[0], dict):
			self.from_dict(args[0])
		else:
			_from_arguments(*args, **kwargs)
	
	def _calc_range(self, level = 0.9545):
		if self.years is None:
			return
		if self.likelihood is not None:
			self._data['likelihood_range'] = calc_range(self.years, self.likelihood, p = level)
		if self.posterior is not None:
			self._data['posterior_range'] = calc_range(self.years, self.posterior, p = level)
	
	def _calc_mean(self):
		if self.years is None:
			return
		if self.likelihood is not None:
			self._data['likelihood_mean'] = calc_mean(self.years, self.likelihood)
		if self.posterior is not None:
			self._data['posterior_mean'] = calc_mean(self.years, self.posterior)
	
	
	# Assigned properties
	
	@property
	def name(self):
		return self._data['name']
	
	@property
	def age(self):
		return self._data['age']
	
	@property
	def uncertainty(self):
		return self._data['uncertainty']
	
	@property
	def date_type(self):
		return self._data['date_type']
	
	@property
	def long_lived(self):
		return self._data['long_lived']
	
	@property
	def context(self):
		return self._data['context']
	
	@property
	def area(self):
		return self._data['area']
	
	@property
	def area_excavation_phase(self):
		return self._data['area_excavation_phase']
	
	@property
	def earlier_than(self):
		if self._data['earlier_than'] is None:
			return None
		return copy.copy(self._data['earlier_than'])
	
	
	# Calculated properties
	
	@property
	def group(self):
		# Group that the sample belongs to based on stratigraphic interconnection with other samples
		return self._data['group']
	
	@property
	def phase(self):
		# Stratigraphic phase of the sample within the group (lower = earlier (older) phase)
		return self._data['phase']
	
	@property
	def years(self):
		# Calendar years BP corresponding to the probability distributions
		# years = np.array([year, ...])
		if self._data['years'] is None:
			return None
		return self._data['years'].copy()
	
	@property
	def likelihood(self):
		# Probability distribution of the dating of the sample before Bayesian modeling
		# likelihood = np.array([p, ...]), where p is the probability of the calendar year
		if self._data['likelihood'] is None:
			return None
		return self._data['likelihood'].copy()
	
	@property
	def posterior(self):
		# Probability distribution of the dating of the sample after Bayesian modeling
		# posterior = np.array([p, ...]), where p is the probability of the calendar year
		if self._data['posterior'] is None:
			return None
		return self._data['posterior'].copy()
	
	@property
	def likelihood_range(self):
		# 2-sigma (95.45%) range of the dating of the sample in calendar years BP before Bayesian modeling
		# likelihood_range = [from, to]
		if self._data['likelihood_range'] is None:
			self._calc_range()
		if self._data['likelihood_range'] is None:
			return [None, None]
		return copy.copy(self._data['likelihood_range'])
	
	@property
	def likelihood_mean(self):
		# Mean calendar age of the sample in calendar years BP before Bayesian modeling
		if self._data['likelihood_mean'] is None:
			self._calc_mean()
		return self._data['likelihood_mean']
	
	@property
	def posterior_range(self):
		# 2-sigma (95.45%) range of the dating of the sample in calendar years BP after Bayesian modeling
		# posterior_range = [from, to]
		if self._data['posterior_range'] is None:
			self._calc_range()
		if self._data['posterior_range'] is None:
			return [None, None]
		return copy.copy(self._data['posterior_range'])
	
	@property
	def posterior_mean(self):
		# Mean calendar age of the sample in calendar years BP after Bayesian modeling
		if self._data['posterior_mean'] is None:
			self._calc_mean()
		return self._data['posterior_mean']
	
	@property
	def posterior_agreement(self):
		return self._data['posterior_agreement']
	
	@property
	def is_calibrated(self):
		return (self._data['likelihood'] is not None)
	
	@property
	def is_modeled(self):
		return (self._data['posterior'] is not None)
	
	
	# Methods
	
	def calibrate(self, curve):
		
		self._data['years'] = curve[:,0].copy()
		self._data['likelihood'] = calibrate(self.age, self.uncertainty, curve, self.date_type)
	
	def set_group(self, group):
		
		self._data['group'] = group
	
	def set_phase(self, phase):
		
		self._data['phase'] = phase
	
	def set_likelihood(self, distribution, mean = None, rng = None):
		# distribution = np.array([p, ...])
		# rng = [min, max]; 2-sigma (95.45%) range in calendar years BP
		self._data['likelihood'] = distribution
		self._data['likelihood_mean'] = mean
		if rng is None:
			self._data['likelihood_range'] = None
		else:
			self._data['likelihood_range'] = copy.copy(rng)
	
	def set_posterior(self, distribution, mean = None, rng = None, agreement = 0):
		# distribution = np.array([p, ...])
		# rng = [min, max]; 2-sigma (95.45%) range in calendar years BP
		# agreement = agreement index generated by OxCal modeling
		self._data['posterior'] = distribution
		self._data['posterior_mean'] = mean
		self._data['posterior_agreement'] = agreement
		if rng is None:
			self._data['posterior_range'] = None
		else:
			self._data['posterior_range'] = copy.copy(rng)
	
	def to_oxcal(self):
		if not self.is_calibrated:
			raise Exception("Sample must be calibrated to generate OxCal definition")
		if self.long_lived and (self.date_type == 'R'):
			# If the sample is long-lived, change it to a uniform distribution
			# Use the oldest date in the 2-sigma range as lower boundary and 0 BP as upper boundary
			lower = self.likelihood_range[0]
			return oxcal_date(self.name, lower / 2, lower / 2, 'U')
		return oxcal_date(self.name, self.age, self.uncertainty, self.date_type)
	
	def to_dict(self):
		data = copy.deepcopy(self._data)
		return dict_np_to_list(data)
	
	def from_dict(self, data):
		self._data = dict(
			name = None,
			age = None,
			uncertainty = None,
			date_type = None,
			long_lived = None,
			context = None,
			area = None,
			area_excavation_phase = None,
			earlier_than = [],
			
			group = None,
			phase = None,
			years = None,
			likelihood = None,
			posterior = None,
			likelihood_range = None,
			likelihood_mean = None,
			posterior_range = None,
			posterior_mean = None,
			posterior_agreement = None,
		)
		self._data.update(dict_keys_to_int(data))
		for key in ['years', 'likelihood', 'posterior']:
			if self._data[key] is not None:
				self._data[key] = np.array(self._data[key], dtype = np.float64)
	
	def copy(self):
		return Sample(self.to_dict())
		
