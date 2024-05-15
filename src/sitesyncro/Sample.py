import copy
from typing import List

import numpy as np

from sitesyncro.utils.fnc_data import (dict_keys_to_int, dict_np_to_list)
from sitesyncro.utils.fnc_oxcal import (oxcal_date)
from sitesyncro.utils.fnc_radiocarbon import (calibrate)
from sitesyncro.utils.fnc_stat import (calc_mean, calc_range)


class Sample(object):
	"""
	A class representing a single radiocarbon sample.
	
	The constructor accepts a single dictionary argument or multiple arguments to initialize the object with the provided data.
	
	Parameters:
	:param name: Sample ID (required, unique identifier)
	:type name: str
	
	:param age: C-14 age (years BP) for date_type 'R'; mean calendar age (years BP) for date_type 'U' (required)
	:type age: float
	
	:param uncertainty: Uncertainty (years BP) for date_type 'R'; 1/2 range (years BP) for date_type 'U' (required)
	:type uncertainty: float
	
	:param date_type: 'R' for radiocarbon date; 'U' for calendar date as a uniform distribution
	:type date_type: str, optional
	
	:param long_lived: True if sample could be older than the examined deposition event due to e.g. old wood effect
	:type long_lived: bool, optional
	
	:param redeposited: True if sample could be redeposited from a different context
	:type redeposited: bool, optional
	
	:param outlier: True if sample is an outlier and should not be used for modeling
	:type outlier: bool, optional
	
	:param context: Name of the context where sample was found
	:type context: str, optional
	
	:param area: Excavation area
	:type area: str, optional
	
	:param excavation_area_phase: Chronological phase of the context within the excavation area (format: "1" or "1a" or "1-2" or "1a-b" or "1a-2b", higher = earlier (older) phase)
	:type excavation_area_phase: str, optional
	
	:param earlier_than: List of names of samples which are stratigraphically earlier (older) than this sample
	:type earlier_than: List[str], optional
	
	:param curve: Radiocarbon calibration curve in format `np.array([[calendar year BP, C-14 year, uncertainty], ...])`
	:type curve: np.ndarray, optional
	"""
	
	def __init__(self, *args, **kwargs):
		
		# Alternatively, if only one argument of type dictionary is provided, the object is initialized using this data
		
		def _from_arguments(name: str, age: float, uncertainty: float,
							date_type: str = 'R',
							long_lived: bool = False,
							redeposited: bool = False,
							outlier: bool = False,
							context: str = None,
							area: str = None,
							excavation_area_phase: str = None,
							earlier_than: List[str] = [],
							curve: np.ndarray = None,
							) -> None:
			
			if date_type not in ['R', 'U']:
				raise Exception("Invalid date type specified: %s (must be 'R' or 'U')" % (date_type))
			
			self._data = dict(
				name=name,
				age=age,
				uncertainty=uncertainty,
				date_type=date_type,
				long_lived=long_lived,
				redeposited=redeposited,
				outlier=outlier,
				context=context,
				area=area,
				excavation_area_phase=excavation_area_phase,
				earlier_than=earlier_than,
				
				group=None,
				phase=None,
				years=None,
				likelihood=None,
				posterior=None,
				likelihood_range=None,
				likelihood_mean=None,
				posterior_range=None,
				posterior_mean=None,
				posterior_agreement=None,
			)
			
			if curve is not None:
				self.calibrate(curve)
		
		self._data = {}
		if len(args) == 1 and isinstance(args[0], dict):
			self.from_dict(args[0])
		else:
			_from_arguments(*args, **kwargs)
	
	def _calc_range(self, level: float = 0.9545) -> None:
		if self.years is None:
			return
		if self.likelihood is not None:
			self._data['likelihood_range'] = calc_range(self.years, self.likelihood, p=level)
		if self.posterior is not None:
			self._data['posterior_range'] = calc_range(self.years, self.posterior, p=level)
	
	def _calc_mean(self) -> None:
		if self.years is None:
			return
		if self.likelihood is not None:
			self._data['likelihood_mean'] = calc_mean(self.years, self.likelihood)
		if self.posterior is not None:
			self._data['posterior_mean'] = calc_mean(self.years, self.posterior)
	
	# Assigned properties
	
	@property
	def name(self) -> str or None:
		"""
		Gets the unique identifier (ID) of the sample.

		:return: The unique identifier of the sample or None if it's not set.
		:rtype: str or None
		"""
		
		return self._data['name']
	
	@property
	def age(self) -> float or None:
		"""
		The C-14 age or mean calendar age of the sample in years BP.

		:return: The C-14 age or mean calendar age of the sample or None if it's not set.
		:rtype: float or None
		"""
		return self._data['age']
	
	@property
	def uncertainty(self) -> float or None:
		"""
		The uncertainty of the sample's C-14 age in years BP.

		:return: The uncertainty of the sample's age or None if it's not set.
		:rtype: float or None
		"""
		return self._data['uncertainty']
	
	@property
	def date_type(self) -> str or None:
		"""
		The type of the date for the sample.

		:return: The type of the date ('R' for radiocarbon date; 'U' for calendar date as a uniform distribution) or None if it's not set.
		:rtype: str or None
		"""
		
		return self._data['date_type']
	
	@property
	def long_lived(self) -> bool or None:
		"""
		Indicates if the sample is long-lived.

		True if the sample could be older than the examined deposition event due to factors like the old wood effect.

		:return: True if the sample is long-lived. Returns None if it's not set.
		:rtype: bool or None
		"""
		
		return self._data['long_lived']
	
	@property
	def redeposited(self) -> bool or None:
		"""
		Indicates if the sample could be redeposited from a different context.

		:return: True if the sample could be redeposited, False otherwise. Returns None if it's not set.
		:rtype: bool or None
		"""
		
		return self._data['redeposited']
	
	@property
	def outlier(self) -> bool or None:
		"""
		Indicates if the sample is an outlier and should not be used for modeling.

		:return: True if the sample is an outlier, False otherwise. Returns None if it's not set.
		:rtype: bool or None
		"""
		
		return self._data['outlier']
	
	@property
	def context(self) -> str or None:
		"""
		The context where the sample was found.

		:return: The name of the context where the sample was found or None if it's not set.
		:rtype: str or None
		"""
		
		return self._data['context']
	
	@property
	def area(self) -> str or None:
		"""
		The excavation area where the sample was found.

		:return: The name of the excavation area where the sample was found or None if it's not set.
		:rtype: str or None
		"""
		
		return self._data['area']
	
	@property
	def excavation_area_phase(self) -> str or None:
		"""
		The chronological phase of the context within the excavation area.
		
		:return: The chronological phase of the context within the excavation area (format: "1" or "1a" or "1-2" or "1a-b" or "1a-2b", higher = earlier (older) phase) or None if it's not set.
		:rtype: str or None
		"""
		
		return self._data['excavation_area_phase']
	
	@property
	def earlier_than(self) -> List[str] or None:
		"""
		List of names of samples which are stratigraphically earlier (older) than this sample.

		:return: List of names of samples which are stratigraphically earlier (older) than this sample or None if it's not set.
		:rtype: List[str] or None
		"""
		
		if self._data['earlier_than'] is None:
			return None
		return copy.copy(self._data['earlier_than'])
	
	# Calculated properties
	
	@property
	def group(self) -> str or None:
		"""
		The group that the sample belongs to based on stratigraphic interconnection with other samples.

		:return: The group that the sample belongs to or None if it's not set.
		:rtype: str or None
		"""
		
		return self._data['group']
	
	@property
	def phase(self) -> float or None:
		"""
		The stratigraphic phase of the sample within the group.

		:return: The stratigraphic phase of the sample within the group (lower = earlier (older) phase) or None if it's not set.
		:rtype: float or None
		"""
		
		return self._data['phase']
	
	@property
	def years(self) -> np.ndarray or None:
		"""
		Calendar years BP corresponding to the probability distributions.

		:return: An array of calendar years BP corresponding to the probability distributions or None if it's not set.
		:rtype: np.ndarray or None
		"""
		
		if self._data['years'] is None:
			return None
		return self._data['years'].copy()
	
	@property
	def likelihood(self) -> np.ndarray or None:
		"""
		The probability distribution of the dating of the sample before Bayesian modeling.

		:return: An array representing the probability distribution of the dating of the sample before Bayesian modeling, where each element is the probability of the corresponding calendar year. Returns None if it's not set.
		:rtype: np.ndarray or None
		"""
		
		if self._data['likelihood'] is None:
			return None
		return self._data['likelihood'].copy()
	
	@property
	def posterior(self) -> np.ndarray or None:
		"""
		The probability distribution of the dating of the sample after Bayesian modeling.

		:return: An array representing the probability distribution of the dating of the sample after Bayesian modeling, where each element is the probability of the corresponding calendar year. Returns None if it's not set.
		:rtype: np.ndarray or None
		"""
		
		if self._data['posterior'] is None:
			return None
		return self._data['posterior'].copy()
	
	@property
	def likelihood_range(self) -> [float, float] or [None, None]:
		"""
		The 2-sigma (95.45%) range of the dating of the sample in calendar years BP before Bayesian modeling.

		:return: A list containing the lower and upper bounds of the 2-sigma range, or [None, None] if it's not set.
		:rtype: list
		"""
		
		if self._data['likelihood_range'] is None:
			self._calc_range()
		if self._data['likelihood_range'] is None:
			return [None, None]
		return copy.copy(self._data['likelihood_range'])
	
	@property
	def likelihood_mean(self) -> float or None:
		"""
		The mean calendar age of the sample in calendar years BP before Bayesian modeling.

		:return: The mean calendar age of the sample before Bayesian modeling or None if it's not set.
		:rtype: float or None
		"""
		
		if self._data['likelihood_mean'] is None:
			self._calc_mean()
		return self._data['likelihood_mean']
	
	@property
	def posterior_range(self) -> [float, float] or [None, None]:
		"""
		The 2-sigma (95.45%) range of the dating of the sample in calendar years BP after Bayesian modeling.

		:return: A list containing the lower and upper bounds of the 2-sigma range, or [None, None] if it's not set.
		:rtype: list
		"""
		
		if self._data['posterior_range'] is None:
			self._calc_range()
		if self._data['posterior_range'] is None:
			return [None, None]
		return copy.copy(self._data['posterior_range'])
	
	@property
	def posterior_mean(self) -> float or None:
		"""
		The mean calendar age of the sample in calendar years BP after Bayesian modeling.

		:return: The mean calendar age of the sample after Bayesian modeling or None if it's not set.
		:rtype: float or None
		"""
		
		if self._data['posterior_mean'] is None:
			self._calc_mean()
		return self._data['posterior_mean']
	
	@property
	def posterior_agreement(self) -> float or None:
		"""
		The agreement index generated by OxCal modeling.

		:return: The agreement index generated by OxCal modeling or None if it's not set.
		:rtype: float or None
		"""
		
		return self._data['posterior_agreement']
	
	@property
	def is_calibrated(self) -> bool:
		"""
		Checks if the sample has been calibrated.

		:return: True if the sample has been calibrated, False otherwise.
		:rtype: bool
		"""
		
		return (self._data['likelihood'] is not None)
	
	@property
	def is_modeled(self) -> bool:
		"""
		Checks if the posterior has been calculated using Bayesian modeling.

		:return: True if the posterior has been calculated, False otherwise.
		:rtype: bool
		"""
		return (self._data['posterior'] is not None)
	
	# Methods
	
	def calibrate(self, curve: np.ndarray) -> None:
		"""
		Calibrates the sample using the provided calibration curve.
		
		:param curve: The calibration curve to be used for calibration. Format: `np.array([[calendar year BP, C-14 year, uncertainty], ...])`.
		:type curve: np.ndarray
		:raises Exception: If the calibration curve is not provided or if the necessary sample data (age, uncertainty, date_type) is not set.
		"""
		
		if curve is None:
			raise Exception("Calibration curve not provided")
		if None in [self.age, self.uncertainty, self.date_type]:
			return
		self._data['years'] = curve[:, 0].copy()
		self._data['likelihood'] = calibrate(self.age, self.uncertainty, curve, self.date_type)
	
	def set_outlier(self, state: bool) -> None:
		"""
		Sets the outlier status of the sample.

		:param state: The outlier status to be set. True if the sample is an outlier and should not be used for modeling, False otherwise.
		:type state: bool
		:raises Exception: If the provided state is not a boolean value.
		"""
		
		if not isinstance(state, bool):
			raise Exception("Outlier state must be a boolean value")
		self._data['outlier'] = state
	
	def set_group(self, group: str) -> None:
		"""
		Sets the group that the sample belongs to based on stratigraphic interconnection with other samples.

		:param group: The group that the sample belongs to.
		:type group: str
		"""
		
		self._data['group'] = group
	
	def set_phase(self, phase: float) -> None:
		"""
		Sets the stratigraphic phase of the sample within the group.

		:param phase: The stratigraphic phase of the sample within the group (lower = earlier (older) phase).
		:type phase: float
		"""
		
		self._data['phase'] = phase
	
	def set_likelihood(self, distribution: np.ndarray, mean: float = None, rng: List[float] = None) -> None:
		"""
		Sets the likelihood for the sample.

		:param distribution: The likelihood distribution for the sample. It should be in the format `np.array([p, ...])`, where p is the probability of the calendar year.
		:type distribution: np.ndarray
		:param mean: The mean of the likelihood distribution. Defaults to None.
		:type mean: float, optional
		:param rng: The 2-sigma (95.45%) range of the likelihood distribution in calendar years BP. It should be in the format `[from, to]`. Defaults to None.
		:type rng: List[float], optional
		"""
		
		self._data['likelihood'] = distribution
		self._data['likelihood_mean'] = mean
		if rng is None:
			self._data['likelihood_range'] = None
		else:
			self._data['likelihood_range'] = copy.copy(rng)
	
	def set_posterior(self, distribution: np.ndarray, mean: float = None, rng: List[float] = None,
					  agreement: float = 0) -> None:
		"""
		Sets the posterior for the sample.

		:param distribution: The posterior distribution for the sample. It should be in the format `np.array([p, ...])`, where p is the probability of the calendar year.
		:type distribution: np.ndarray
		:param mean: The mean of the posterior distribution. Defaults to None.
		:type mean: float, optional
		:param rng: The 2-sigma (95.45%) range of the posterior distribution in calendar years BP. It should be in the format `[from, to]`. Defaults to None.
		:type rng: List[float], optional
		:param agreement: The agreement index generated by OxCal modeling. Defaults to 0.
		:type agreement: float, optional
		"""
		
		self._data['posterior'] = distribution
		self._data['posterior_mean'] = mean
		self._data['posterior_agreement'] = agreement
		if rng is None:
			self._data['posterior_range'] = None
		else:
			self._data['posterior_range'] = copy.copy(rng)
	
	def get_range(self) -> [float, float] or [None, None]:
		"""
		Gets the 2-sigma (95.45%) range of the dating of the sample in calendar years BP.

		If the sample has been modeled, it returns the posterior range.
		If the sample is not calibrated, it returns [None, None].
		If the sample is long-lived and the date type is 'R', it keeps the lower boundary and doubles the range.
		Otherwise, it returns the likelihood range.

		:return: A list containing the lower and upper bounds of the 2-sigma range, or [None, None] if it's not set.
		:rtype: list
		"""
		
		if self.is_modeled:
			return self.posterior_range
		if not self.is_calibrated:
			return [None, None]
		if self.long_lived and (self.date_type == 'R'):
			# If the sample is long-lived, keep the lower boundary and double the range
			lower, upper = self.likelihood_range
			if lower is None:
				return [None, None]
			return [lower, lower - 2 * (lower - upper)]
		return self.likelihood_range
	
	def to_oxcal(self) -> str or None:
		"""
		Converts the sample to OxCal model format.

		:return: The sample in OxCal model format or None if the necessary sample data (name, age, uncertainty, date_type, long_lived, outlier) is not set.
		:rtype: str or None
		"""
		
		return oxcal_date(self.name, self.age, self.uncertainty, self.date_type, self.long_lived, self.outlier)
	
	def to_dict(self) -> dict:
		"""
		Converts the sample data to a dictionary.

		The method creates a deep copy of the sample data and converts numpy arrays to lists for JSON serialization.

		:return: A dictionary containing the sample data.
		:rtype: dict
		"""
		
		data = copy.deepcopy(self._data)
		return dict_np_to_list(data)
	
	def from_dict(self, data: dict) -> None:
		"""
		Loads the sample data from a dictionary.

		The method updates the sample data with the provided dictionary. It also converts keys to integers and numpy arrays to lists for JSON serialization.

		:param data: A dictionary containing the sample data.
		:type data: dict
		"""
		
		self._data = dict(
			name=None,
			age=None,
			uncertainty=None,
			date_type=None,
			long_lived=None,
			redeposited=None,
			outlier=None,
			context=None,
			area=None,
			excavation_area_phase=None,
			earlier_than=[],
			
			group=None,
			phase=None,
			years=None,
			likelihood=None,
			posterior=None,
			likelihood_range=None,
			likelihood_mean=None,
			posterior_range=None,
			posterior_mean=None,
			posterior_agreement=None,
		)
		self._data.update(dict_keys_to_int(data))
		for key in ['years', 'likelihood', 'posterior']:
			if self._data[key] is not None:
				self._data[key] = np.array(self._data[key], dtype=np.float64)
	
	def copy(self) -> 'Sample':
		"""
		Creates a copy of the sample instance.

		:return: A new instance of the Sample class with the same data as the current instance.
		:rtype: Sample
		"""
		
		return Sample(self.to_dict())
	
	def __repr__(self) -> str:
		repr_str = f"<Sample '{self.name}': age={self.age}, uncertainty={self.uncertainty}, date_type={self.date_type}, long_lived={self.long_lived}, redeposited={self.redeposited}, outlier={self.outlier}, context={self.context}, area={self.area}, excavation_area_phase={self.excavation_area_phase}, earlier_than={self.earlier_than}, group={self.group}, phase={self.phase}"
		
		if self.is_calibrated:
			repr_str += f", likelihood_range={self.likelihood_range}, likelihood_mean={self.likelihood_mean}"
		
		if self.is_modeled:
			repr_str += f", posterior_range={self.posterior_range}, posterior_mean={self.posterior_mean}, posterior_agreement={self.posterior_agreement}"
		
		repr_str += ">"
		
		return repr_str
