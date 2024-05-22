from sitesyncro.utils.fnc_oxcal import (get_distributions)

from typing import List

import numpy as np
import copy

class Phase(object):
	"""
	A class representing a chronological phase.
	Parameters:
	:param group: Stratigraphic group number.
	:type group: int
	:param phase: Phase number.
	:type phase: int
	
	:param phase: Phase nr. (higher = later/younger)
	:type phase: int
	"""
	
	def __init__(self, group: int, phase: int):
		
		self._data = dict(
			group=group,
			phase=phase,
			
			samples=[],
			years=None,
			
			start=None,
			start_range=None,
			start_mean=None,
			
			end=None,
			end_range=None,
			end_mean=None,
		)
	
	def populate(self, model: 'Model'):
		"""
		Populates the phase with probability distributions for start and end.
		
		:param model: Model object
		:type model: Model
		:return: True if successfully populated
		:rtype: bool
		"""
		self._data['years'] = model.years
		_, posteriors = get_distributions(model.oxcal_data, model.curve)
		# posteriors = {name: [distribution, mean, rng, agreement], ...}
		found = []
		for name in posteriors:
			grph = None
			typ = None
			if name.startswith("Start Gr."):
				grph = name[9:].split('-')
				typ = 'start'
			elif name.startswith("End Gr."):
				grph = name[7:].split('-')
				typ = 'end'
			if (grph is None) or (len(grph) != 2):
				continue
			gr, ph = grph[0], grph[1]
			if (not gr.isnumeric()) or (not ph.isnumeric()):
				continue
			gr, ph = int(gr), int(ph)
			if (gr!=self.group) or (ph!=self.phase):
				continue
			if typ is None:
				continue
			found.append(typ)
			self._data[typ] = posteriors[name][0]
			self._data[typ+'_mean'] = posteriors[name][1]
			self._data[typ+'_range'] = posteriors[name][2]
		
		self._data['samples'] = []
		for name in model.samples:
			if model.samples[name].group == self.group and model.samples[name].phase == self.phase:
				self._data['samples'].append(name)
		
		return (len(found) == 2)			
	
	
	@property
	def group(self) -> int:
		"""
		The stratigraphic group the phase belongs to
		:return: The number of the group
		:rtype: int
		"""
		return self._data['group']
	
	@property
	def phase(self) -> int:
		"""
		The phase number.
		:return: The number of the phase
		:rtype: int
		"""
		return self._data['phase']
	
	@property
	def start(self) -> np.ndarray | None:
		"""
		The probability distribution of the dating of the start of the phase.
		"""
		return self._data['start']
	
	@property
	def start_mean(self) -> float | None:
		"""
		The mean of the probability distribution of the dating of the start of the phase.
		"""
		return self._data['start_mean']
	
	@property
	def start_range(self) -> float | None:
		"""
		The range of the probability distribution of the dating of the start of the phase.
		"""
		return self._data['start_range']
	
	@property
	def end(self) -> np.ndarray | None:
		"""
		The probability distribution of the dating of the end of the phase.
		"""
		return self._data['end']
	
	@property
	def end_mean(self) -> float | None:
		"""
		The mean of the probability distribution of the dating of the end of the phase.
		"""
		return self._data['end_mean']
	
	@property
	def end_range(self) -> float | None:
		"""
		The range of the probability distribution of the dating of the end of the phase.
		"""
		return self._data['end_range']
	
	@property
	def samples(self) -> List[str]:
		"""
		The list of samples that belong to the phase.
		"""
		return copy.copy(self._data['samples'])
	
	@property
	def years(self) -> np.ndarray | None:
		"""
		The years of the model.
		"""
		return self._data['years']
	
	def __repr__(self):
		return f"Phase {self.group}-{self.phase}"
	
	