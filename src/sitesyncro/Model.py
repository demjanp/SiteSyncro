import copy
import gzip
import json
import os
import subprocess
from collections import defaultdict
from typing import List, Dict, Any

import numpy as np

from sitesyncro.Sample import Sample
from sitesyncro.utils.fnc_cluster import (proc_clustering)
from sitesyncro.utils.fnc_data import (dict_keys_to_int, dict_np_to_list)
from sitesyncro.utils.fnc_load import (load_data)
from sitesyncro.utils.fnc_oxcal import (download_oxcal, gen_oxcal_model, load_oxcal_data, get_distributions)
from sitesyncro.utils.fnc_phase import (create_earlier_than_matrix, get_groups_and_phases, find_dating_outliers,
										update_earlier_than_by_clustering, update_earlier_than_by_dating)
from sitesyncro.utils.fnc_plot import (plot_randomized, plot_clusters, save_results_csv, save_outliers)
from sitesyncro.utils.fnc_radiocarbon import (get_curve)
from sitesyncro.utils.fnc_simulate import (test_distributions)


class Model(object):
	"""
	A class representing a Bayesian model of dated samples interconnected by stratigraphic relationships.

	:param directory: Working directory for model data (default is "model").
	:type directory: str

	:param samples: List of samples as instances of the class [Sample](#sample_class)
	:type samples: list

	:param curve_name: The name of the calibration curve to use (default is "intcal20.14c").
	:type curve_name: str

	:param phase_model: OxCal phase model type. Can be 'sequence', 'contiguous', 'overlapping', or 'none' (default is "sequence").
	:type phase_model: str

	:param cluster_n: Number of clusters to form (-1 = automatic; default is -1).
	:type cluster_n: int

	:param min_years_per_cluster: The minimum number of years per cluster.
	:type min_years_per_cluster: int

	:param uniform: Flag indicating whether to use uniform randomization (default is False).
	:type uniform: bool

	:param p_value: The P-value for statistical tests (default is 0.05).
	:type p_value: float

	:param uncertainty_base: The base uncertainty for randomization (default is 15).
	:type uncertainty_base: float

	:param npass: Minimum number of passes for the randomization tests (default is 100).
	:type npass: int

	:param convergence: Convergence threshold for the randomization tests (default is 0.99).
	:type convergence: float

	:param oxcal_url: Url to download the OxCal program (default is "https://c14.arch.ox.ac.uk/OxCalDistribution.zip").
	:type oxcal_url: str

	:param overwrite: Flag indicating whether to overwrite existing data in the model directory (default is False).
	:type overwrite: bool
	"""
	
	def __init__(self, **kwargs):

		defaults = dict(
			directory='model',
			samples=[],
			curve_name='intcal20.14c',
			phase_model='sequence',
			cluster_n=-1,
			min_years_per_cluster=25,
			uniform=False,
			p_value=0.05,
			uncertainty_base=15,
			npass=100,
			convergence=0.99,
			oxcal_url='https://c14.arch.ox.ac.uk/OxCalDistribution.zip',
			overwrite=False,
		)
		
		# Check arguments
		for key in kwargs:
			if key not in defaults:
				raise Exception("Invalid argument: %s" % key)
			if type(kwargs[key]) != type(defaults[key]):
				raise Exception("Invalid argument type for %s: %s" % (key, type(kwargs[key].__name__)))
		
		overwrite = kwargs.pop('overwrite', False)
		
		for sample in kwargs.get('samples', []):
			if not isinstance(sample, Sample):
				raise Exception("Invalid sample format: %s. Sample expected." % (type(sample).__name__))
		
		if kwargs.get('phase_model', 'sequence') not in ['sequence', 'contiguous', 'overlapping', 'none']:
			raise Exception("Invalid phase model specified")
		
		if not download_oxcal(url=kwargs.get('oxcal_url', defaults['oxcal_url'])):
			raise Exception("OxCal not found")
		
		self._data = self._assigned()
		self._data.update(self._calculated())
		
		# Attempt to load data from directory
		if ('directory' in kwargs) and (not overwrite) and ('samples' not in kwargs):
			if self.load(kwargs['directory']):
				del kwargs['directory']
				self.update_params(**kwargs)
				return
		
		# Update missing keys in kwargs with defaults
		kwargs = dict([(key, kwargs[key] if key in kwargs else defaults[key]) for key in defaults])
		
		# Convert samples to dict
		kwargs['samples'] = dict([(sample.name, sample) for sample in kwargs['samples']])
		
		self._data.update(kwargs)
		
		# Create model directory if needed
		self._data['directory'] = self._create_dir(self.directory, overwrite)
	
	def _assigned(self) -> Dict[str, Any]:
		return dict(
			directory=None,
			samples=None,
			curve_name=None,
			phase_model=None,
			cluster_n=None,
			min_years_per_cluster=None,
			uniform=None,
			p_value=None,
			uncertainty_base=None,
			npass=None,
			convergence=None,
			oxcal_url=None,
		)
	
	def _calculated(self) -> Dict[str, Any]:
		return dict(
			years=None,
			curve=None,
			oxcal_data=None,
			
			outlier_candidates=None,
			
			summed=None,
			random_p=None,
			random_lower=None,
			random_upper=None,
			
			clusters=None,
			cluster_means=None,
			cluster_sils=None,
			cluster_ps=None,
			cluster_opt_n=None,
		)
	
	def _calculated_randomization(self) -> List[str]:
		
		return ['summed', 'random_p', 'random_lower', 'random_upper']
	
	def _calculated_clustering(self) -> List[str]:
		
		return ['clusters', 'cluster_means', 'cluster_sils', 'cluster_ps', 'cluster_opt_n']
	
	def _create_dir(self, directory: str, overwrite: bool) -> str:
		if os.path.isdir(directory) and (overwrite or (not os.listdir(directory))):
			return directory
		
		parent_dir = os.path.dirname(directory)
		if not parent_dir:
			parent_dir = "."
		n = None
		for d in os.listdir(parent_dir):
			if d == os.path.basename(directory):
				n = 0
			elif d.startswith(os.path.basename(directory)):
				try:
					n = max(n, int(d.split("_")[-1]))
				except:
					pass
		if n is not None:
			n += 1
			directory = os.path.join(parent_dir, os.path.basename(directory) + "_" + str(n))
		os.makedirs(directory)
		return directory
	
	# Assigned properties
	
	@property
	def directory(self) -> str:
		"""
		The directory where the model data is stored.
		
		Returns:
			str: The directory path as a string.
		"""
		return self._data['directory']
	@property
	def directory(self) -> str:
		"""
		Represents the directory where the model data is stored.

		:return: The directory where the model data is stored.
		:rtype: str
		"""
		
		return self._data['directory']
	@property
	def samples(self) -> Dict[str, Sample]:
		"""
		A dictionary of samples associated with the model.
		
		:return: A dictionary where the keys are the sample names and the values are Sample objects.
		:rtype: Dict[str, Sample]
		"""
		if self._data['samples'] is None:
			return {}
		return self._data['samples']
	
	@property
	def curve_name(self) -> str:
		"""
		File name of the radiocarbon age calibration curve (see OxCal/bin directory).

		:return: The name of the calibration curve.
		:rtype: str
		"""
		return self._data['curve_name']
	
	@property
	def phase_model(self) -> str:
		"""
		OxCal phase model type.

		:return: The type of the phase model. Can be 'sequence', 'contiguous', 'overlapping', or 'none'.
		:rtype: str
		"""
		return self._data['phase_model']
	
	@property
	def cluster_n(self) -> int:
		"""
		Number of clusters to form (-1 = automatic).
		
		:return: The number of clusters to form.
		:rtype: int
		"""
		return self._data['cluster_n']
	
	@property
	def min_years_per_cluster(self) -> int:
		"""
		The minimum number of years per cluster.
		
		:return: The minimum number of years per cluster.
		:rtype: int
		"""
		return self._data['min_years_per_cluster']
	
	@property
	def uniform(self) -> bool:
		"""
		Flag indicating whether the model tests for a uniform distribution of the calendar ages.

		:return: True if a uniform distribution is used, False otherwise.
		:rtype: bool
		"""
		return self._data['uniform']
	
	@property
	def p_value(self) -> float:
		"""
		The p-value used for the randomization test.

		:return: The P-value for statistical tests.
		:rtype: float
		"""
		return self._data['p_value']
	
	@property
	def uncertainty_base(self) -> float:
		"""
		The base uncertainty for the radiocarbon dates.

		:return: The base uncertainty for the radiocarbon dates.
		:rtype: float
		"""
		return self._data['uncertainty_base']
	
	@property
	def npass(self) -> int:
		"""
		The minimum number of passes for the randomization test.
		
		:return: The minimum number of passes.
		:rtype: int
		"""
		return self._data['npass']
	
	@property
	def convergence(self) -> float:
		"""
		The convergence threshold for the randomization test.

		:return: The convergence threshold for the randomization tests.
		:rtype: float
		"""
		return self._data['convergence']
	
	@property
	def oxcal_url(self) -> str:
		"""
		The URL from where the OxCal software can be downloaded.

		:return: The URL of the OxCal software.
		:rtype: str
		"""
		return self._data['oxcal_url']
	
	# Calculated properties
	
	@property
	def years(self) -> np.ndarray:
		"""
		Calendar years BP corresponding to the probability distributions.

		:return: An array of calendar years.
		:rtype: np.ndarray
		"""
		if self._data['years'] is None:
			if self.curve is not None:
				self._data['years'] = self.curve[:, 0]
			else:
				return None
		return self._data['years'].copy()
	
	@property
	def curve(self) -> np.ndarray:
		"""
		Radiocarbon calibration curve.

		2D array containing the calibration curve data. Each row represents a calendar year BP, C-14 year, and uncertainty.
		
		:return: An array of the calibration curve.
		:rtype: np.ndarray
		"""
		
		if self._data['curve'] is None:
			if self._data['curve_name'] is None:
				return None
			else:
				self._data['curve'] = get_curve(self._data['curve_name'])
		return self._data['curve'].copy()
	
	@property
	def uncertainties(self) -> List[float]:
		"""
		List of uncertainties from C-14 dates of samples.

		:return: A list of uncertainties from C-14 dates of samples.
		:rtype: List[float]
		"""
		uncertainties = []
		for name in self.samples:
			if self.samples[name].date_type == 'R':
				uncertainties.append(self.samples[name].uncertainty)
		return uncertainties
	
	@property
	def oxcal_data(self) -> dict or None:
		"""
		The OxCal data associated with the model.

		:return: A dictionary containing the OxCal data if it exists, otherwise None.
		:rtype: dict or None
		"""
		if self._data['oxcal_data'] is None:
			return None
		return copy.deepcopy(self._data['oxcal_data'])
	
	@property
	def outliers(self) -> List[str]:
		"""
		List of outliers among samples which need to be removed for the model to be valid.

		:return: A list of sample names that are considered outliers.
		:rtype: List[str]
		"""
		return [name for name in self.samples if self.samples[name].outlier]
	
	@property
	def outlier_candidates(self) -> List[str]:
		"""
		List of candidates for outliers, from which the final outliers to be eliminated were picked.
		These samples have conflicts between dating ranges and stratigraphic relationships with other samples.

		:return: A list of sample names that are considered candidates for outliers.
		:rtype: List[str]
		"""
		if self._data['outlier_candidates'] is None:
			return []
		return copy.deepcopy(self._data['outlier_candidates'])
	
	@property
	def summed(self) -> np.ndarray or None:
		"""
		Summed probability distribution of the dating of all samples.
		The summed probability is represented as an array where each element is the probability of the calendar year.

		:return: An array of the summed probability distribution if it exists, otherwise None.
		:rtype: np.ndarray or None
		"""
		if self._data['summed'] is None:
			return None
		return self._data['summed'].copy()
	
	@property
	def random_p(self) -> float or None:
		"""
		The calculated p-value from the randomization test.

		:return: The p-value if the randomization test has been performed, otherwise None.
		:rtype: float or None
		"""
		return self._data['random_p']
	
	@property
	def random_lower(self) -> np.ndarray or None:
		"""
		Lower bound of the randomization test.
		The lower bound is represented as an array where each element is the probability of the calendar year.

		:return: An array of the lower bound of the randomization test if it exists, otherwise None.
		:rtype: np.ndarray or None
		"""
		if self._data['random_lower'] is None:
			return None
		return self._data['random_lower'].copy()
	
	@property
	def random_upper(self) -> np.ndarray or None:
		"""
		Upper bound of the randomization test.
		The Upper bound is represented as an array where each element is the probability of the calendar year.

		:return: An array of the upper bound of the randomization test if it exists, otherwise None.
		:rtype: np.ndarray or None
		"""
		if self._data['random_upper'] is None:
			return None
		return self._data['random_upper'].copy()
	
	@property
	def areas(self) -> List[str]:
		"""
		Unique area names extracted from the samples.

		:return: A sorted list of unique area names.
		:rtype: List[str]
		"""
		areas = set()
		for name in self.samples:
			areas.add(self.samples[name].area)
		return sorted(list(areas))
	
	@property
	def contexts(self) -> List[str]:
		"""
		Unique context names extracted from the samples.

		:return: A sorted list of unique context names.
		:rtype: List[str]
		"""
		
		contexts = set()
		for name in self.samples:
			contexts.add(self.samples[name].context)
		return sorted(list(contexts))
	
	@property
	def groups(self) -> Dict[str, List[str]]:
		"""
		Groups that the samples belong to based on stratigraphic interconnection with other samples.
		The groups are represented as a dictionary where the keys are the group names and the values are lists of sample names.

		:return: A dictionary where the keys are the group names and the values are lists of sample names.
		:rtype: Dict[str, List[str]]
		"""
		# groups = {group_name: [sample name, ...], ...}
		
		groups = defaultdict(list)
		for name in self.samples:
			groups[self.samples[name].group].append(name)
		return dict(groups)
	
	@property
	def clusters(self) -> Dict[int, Dict[int, List[str]]]:
		"""
		Clusters of samples based on the similarity of their probability distributions.

		:return: {clusters_n: [cluster: [sample name, ...], ...}, ...}; clusters_n = number of clusters
		:rtype: Dict[int, Dict[int, List[str]]]
		"""
		
		if self._data['clusters'] is None:
			return None
		return copy.deepcopy(self._data['clusters'])
	
	@property
	def cluster_means(self) -> Dict[int, Dict[int, float]]:
		"""
		Mean date of the samples in each cluster in calendar years BP.
		
		:return: {clusters_n: {cluster: year, ...}, ...}; clusters_n = number of clusters
		:rtype: Dict[int, Dict[int, float]]
		"""
		
		if self._data['cluster_means'] is None:
			return None
		return copy.deepcopy(self._data['cluster_means'])
	
	@property
	def cluster_sils(self) -> Dict[int, float]:
		"""
		Silhouette scores of each clustering solution.

		The silhouette score is a measure of how similar an object is to its own cluster compared to other clusters.

		:return: {clusters_n: silhouette, ...}; clusters_n = number of clusters
		:rtype: Dict[int, float]
		"""
		
		if self._data['cluster_sils'] is None:
			return None
		return copy.deepcopy(self._data['cluster_sils'])
	
	@property
	def cluster_ps(self) -> Dict[int, float]:
		"""
		P-values of the clustering solutions.

		:return: {clusters_n: p, ...}; clusters_n = number of clusters
		:rtype: Dict[int, float]
		"""
		
		if self._data['cluster_ps'] is None:
			return None
		return copy.deepcopy(self._data['cluster_ps'])
	
	@property
	def cluster_opt_n(self) -> int or None:
		"""
		Optimal number of clusters based on the silhouette scores.

		The optimal number of clusters is the one that maximizes the average silhouette score for clustering solutions for which the p-value is lower than Model.p_value.

		:return: The optimal number of clusters if the clustering has been performed, None otherwise.
		:rtype: int or None
		"""
		
		return self._data['cluster_opt_n']
	
	@property
	def has_data(self) -> bool:
		"""
		Checks if the model has any associated sample data.

		:return: True if the model has sample data, False otherwise.
		:rtype: bool
		"""
		
		return (len(self.samples) > 0)
	
	@property
	def is_modeled(self) -> bool:
		"""
		Checks if Bayesian modeling of sample dates has been performed for all samples.

		:return: True if Bayesian modeling has been performed for all samples, False otherwise.
		:rtype: bool
		"""
		
		if not self.has_data:
			return False
		for name in self.samples:
			if not self.samples[name].is_modeled:
				return False
		return True
	
	@property
	def is_randomized(self) -> bool:
		"""
		Checks if the randomization test has been performed.
		
		:return: True if the randomization test has been performed, False otherwise.
		:rtype: bool
		"""
		return self._data['random_p'] is not None
	
	@property
	def is_clustered(self) -> bool:
		"""
		Checks if the model has been clustered.
		
		:return: True if the model has been clustered, False otherwise.
		:rtype: bool
		"""
		
		if self._data['clusters']:
			return True
		return False
	
	# Methods
	
	def add_sample(self, *args, **kwargs) -> None:
		"""
		Adds a sample to the model.

		Accepts either a single argument of type Sample or a set of arguments and keyword arguments to create a new Sample instance.

		If a single argument of type Sample is provided, it is added to the model's samples directly.
		If multiple arguments are provided, they are used to create a new Sample instance which is then added to the model's samples.

		:param args: Either a single argument of type Sample or multiple arguments to create a new Sample instance.
		:param kwargs: Keyword arguments to create a new Sample instance.
		:return: None
		"""
		
		def _from_arguments(name: str, age: float, uncertainty: float,
							date_type: str = 'R', long_lived: bool = False, redeposited: bool = False,
							outlier: bool = False, context: bool = None,
							area: str = None,
							excavation_area_phase: float = None,
							earlier_than: List[str] = [],
							):
			self._data['samples'][name] = Sample(
				name, age, uncertainty, date_type, long_lived, redeposited, outlier,
				context, area, excavation_area_phase, earlier_than, self.curve
			)
			self._data['samples'][name].calibrate(self.curve)
		
		self.reset_model()
		if len(args) == 1 and isinstance(args[0], Sample):
			sample = args[0].copy()
			if not sample.is_calibrated:
				sample.calibrate(self.curve)
			self._data['samples'][sample.name] = sample
		else:
			_from_arguments(*args, **kwargs)
	
	def del_sample(self, name: str) -> None:
		"""
		Deletes a sample from the model.

		Removes a sample from the model's samples based on the provided sample name.

		:param name: The name of the sample to be deleted.
		:type name: str
		:return: None
		"""
		
		if name in self._data['samples']:
			del self._data['samples'][name]
	
	def save(self, zipped: bool = False) -> None:
		"""
		Saves the model to a JSON file.

		The file is saved in the model's directory.
		
		:param zipped: If True, the model is saved as a zipped JSON file. Defaults to False.
		:type zipped: bool
		:return: None
		"""
		
		fname = os.path.join(self.directory, 'model.json.gz' if zipped else 'model.json')
		data = dict_np_to_list(copy.deepcopy(self._data))
		for name in self.samples:
			data['samples'][name] = self.samples[name].to_dict()
		
		if zipped:
			with gzip.open(fname, 'wt') as file:
				json.dump(data, file)
		else:
			with open(fname, 'w') as file:
				json.dump(data, file)
	
	def copy(self, directory: str) -> object:
		"""
		Creates a copy of the current model with a new directory.

		Creates a new instance of the Model class with the same parameters and data as the current model,
		but with a different directory. The new directory is provided as an argument.

		:param directory: The directory for the new model.
		:type directory: str
		:return: A new instance of the Model class with the same data as the current model but a different directory.
		:rtype: Model

		"""
		
		samples = dict([(name, self.samples[name].copy()) for name in self.samples])
		model = Model(directory, samples, self.curve_name, self.phase_model,
					  self.cluster_n, self.min_years_per_cluster, self.uniform, self.p_value, self.uncertainty_base,
					  self.oxcal_url
					  )
		for key in self._calculated():
			model._data[key] = getattr(self, key)
		return model
	
	def load(self, directory: str = None) -> bool:
		"""
		Loads the model from a JSON file.

		Attempts to load the model data from a JSON file located in the provided directory.
		If no directory is provided, it uses the model's current directory.
		Supports both regular and zipped JSON files.

		:param directory: The directory from where the model data should be loaded. If None, the model's current directory is used.
		:type directory: str, optional
		:return: True if the model data was successfully loaded, False otherwise.
		:rtype: bool
		"""
		
		if directory is None:
			directory = self.directory
		
		# Check if model file exists
		zipped = False
		fname = os.path.join(directory, 'model.json')
		if not os.path.isfile(fname):
			zipped = True
			fname = os.path.join(directory, 'model.json.gz')
		if not os.path.isfile(fname):
			return False
		
		# Load data
		if zipped:
			with gzip.open(fname, 'rt') as file:
				data = json.load(file)
		else:
			with open(fname, 'r') as file:
				data = json.load(file)
		data = dict_keys_to_int(data)
		
		# Check data structure
		for key in self._assigned():
			if key not in data:
				return False
		for key in self._calculated():
			if key not in data:
				return False
		
		self.update_params(**data)
		
		# Initialize samples
		for name in data['samples']:
			data['samples'][name] = Sample(data['samples'][name])
			if not data['samples'][name].is_calibrated:
				data['samples'][name].calibrate(self.curve)
		
		# Convert lists to numpy arrays
		for key in ['years', 'curve', 'summed', 'random_lower', 'random_upper']:
			if data[key] is not None:
				data[key] = np.array(data[key], dtype=np.float64)
		
		self._data = data
		self._data['directory'] = directory
		return True
	
	def import_csv(self, fname: str) -> None:
		"""
		Loads sample data from a CSV file.

		The file should be in the following format:
		- Each line represents a data record.
		- Data fields are separated by semicolons.
		- The first line is a header and is skipped.
		- Each line should have 10 fields: Sample, Context, Area, C14 Age, Uncertainty, Phase, Earlier-Than, Long-Lived, Redeposited, Outlier.
		- The fields are processed as follows:
			- Sample, Context, and Area are stripped of leading and trailing whitespace.
			- C14 Age and Uncertainty are converted to floats.
			- Phase is converted to a float if possible, otherwise it is set to None.
			- Earlier-Than is split on commas and stripped of leading and trailing whitespace. If it is empty, it is set to an empty list.
			- Long-Lived, Redeposited, and Outlier are converted to integers and then to booleans.

		:param fname: The file path of the CSV file to be imported.
		:type fname: str
		:raises ValueError: If the input file does not exist or is not formatted correctly.
		:return: None
		"""
		
		# Check if the input file exists
		if not os.path.isfile(fname):
			raise ValueError("Input file %s not found" % fname)
		
		samples, contexts, context_area, long_lived, redeposited, outlier, r_dates, context_phase, earlier_than = load_data(
			fname)
		self._data['samples'] = {}
		for name in samples:
			age, uncertainty = r_dates[name]
			self.add_sample(
				name, age, uncertainty, 'R', long_lived[name], redeposited[name], outlier[name], contexts[name],
				context_area[contexts[name]], context_phase[contexts[name]], earlier_than[name]
			)
	
	def plot_randomized(self, fname: str = None, show: bool = False) -> str:
		"""
		Plots the results of the randomization test.

		The plot can either be saved to a file or displayed.

		Args:
			fname (str, optional): The file name to save the plot to. If None, a default file name is used. Defaults to None.
			show (bool, optional): If True, the plot is displayed. Defaults to False.

		:param fname: The file name to save the plot to. If None, a default file name is used. Defaults to None.
		:type fname: str, optional
		:param show: If True, the plot is displayed. Defaults to False.
		:type show: bool, optional
		:return: The file name the plot was saved to.
		:rtype: str
		"""
		
		if not self.is_randomized:
			print("\nNo randomization data to plot")
			return None
		
		# Determine the file path for the plot
		if fname is None:
			fname = os.path.join(self.directory, "randomized.pdf")
		
		plot_randomized(self, fname, show)
		return fname
	
	def plot_clusters(self, fname: str = None, show: bool = False) -> str:
		"""
		Plots the clustering results.

		The plot can either be saved to a file or displayed.

		:param fname: The file name to save the plot to. If None, a default file name is used. Defaults to None.
		:type fname: str, optional
		:param show: If True, the plot is displayed. Defaults to False.
		:type show: bool, optional
		:return: The file name the plot was saved to.
		:rtype: str
		"""
		
		if not self.is_clustered:
			print("\nNo clustering data to plot")
			return None
		
		# Determine the file path for the plot
		if fname is None:
			fname = os.path.join(self.directory, "silhouette.pdf")
		
		# Plot the clustering data
		plot_clusters(self, fname, show)
		return fname
	
	def save_csv(self, fcsv: str = None) -> str:
		"""
		Saves the results to a CSV file.

		:param fcsv: The file path for the CSV file. If None, a default file name and path are used. Defaults to None.
		:type fcsv: str, optional
		:return: The file path the results were saved to.
		:rtype: str
		"""
		
		# Check if there are any results
		if not self.has_data:
			print("\nNo data available")
			return None
		
		# Determine the file path for the CSV file
		if fcsv is None:
			fcsv = os.path.join(self.directory, "results.csv")
		
		# Save the results to the CSV file
		save_results_csv(self, fcsv)
		return fcsv
	
	def save_outliers(self, fname: str = None) -> str:
		"""
		Saves a list of outliers to a text file.

		:param fname: The file name to save the outliers to. If None, a default file name is used. Defaults to None.
		:type fname: str, optional
		:return: The file name the outliers were saved to.
		:rtype: str
		"""
		
		if fname is None:
			fname = os.path.join(self.directory, "outliers.txt")
		
		save_outliers(self, fname)
		return fname
	
	def to_oxcal(self, fname: str = None) -> str:
		"""
		Exports the model to an OxCal file.

		Generates an OxCal file from the current model. The OxCal file can be used for further analysis in the OxCal software.
		If a file name is provided, the OxCal file is saved to that file. If no file name is provided, a default file name is used.

		:param fname: The file name to save the OxCal file to. If None, a default file name is used. Defaults to None.
		:type fname: str, optional
		:return: The file name the OxCal file was saved to.
		:rtype: str
		"""
		
		# Check if there is any data
		if not self.has_data:
			print("\nNo data available")
			return None
		
		if fname is None:
			fname = os.path.join(self.directory, "model.oxcal")
		
		txt = gen_oxcal_model(self)
		with open(fname, "w") as file:
			file.write(txt)
		return fname
	
	def load_oxcal_data(self) -> None:
		"""
		Loads the OxCal data associated with the model.

		Attempts to load the OxCal data from a file named 'model.js' located in the model's directory.
		The loaded data is then stored in the model's 'oxcal_data' attribute.

		:return: None
		"""
		
		fname = os.path.join(self.directory, "model.js")
		data = load_oxcal_data(fname)
		self._data['oxcal_data'] = data
		
		likelihoods, posteriors = get_distributions(data, self.curve)
		# likelihoods = {name: [distribution, mean, rng], ...}
		# posteriors = {name: [distribution, mean, rng, agreement], ...}
		
		for name in self.samples:
			if name in likelihoods:
				self.samples[name].set_likelihood(*likelihoods[name])
			else:
				self.samples[name].set_likelihood(None)
			if name in posteriors:
				self.samples[name].set_posterior(*posteriors[name])
			else:
				self.samples[name].set_posterior(None)
	
	def reset_model(self) -> None:
		"""
		Resets the model.

		Resets all calculated properties of the model to their initial state.
		It also sets the posterior distribution of each sample in the model to None.

		:return: None
		"""
		self._data.update(self._calculated())
		for name in self.samples:
			self.samples[name].set_posterior(None)
	
	def update_params(self, **kwargs) -> (Dict[str, List], set):
		"""
		Updates the model parameters and resets calculated attributes if necessary.

		Accepts keyword arguments that correspond to the model parameters.
		If a parameter is provided that differs from the current value, the model parameter is updated and all related calculated attributes are reset.

		:param kwargs: Keyword arguments corresponding to the model parameters.
		:type kwargs: dict
		:return: (reset_assigned, reset_calculated);
			reset_assigned: {parameter: [old value, new value], ...}; parameters and their values that have been updated;
			reset_calculated: {attribute, ...}; calculated attributes that have been reset
		:rtype: (Dict[str, List], set)
		"""
		
		assigned_full = ['samples', 'curve_name', 'phase_model']
		assigned_randomization = ['uniform', 'p_value', 'uncertainty_base', 'npass', 'convergence']
		assigned_clustering = ['cluster_n', 'min_years_per_cluster']
		
		def _get_n_samples(value):
			if isinstance(value, dict):
				return len(value)
			return 0
		
		reset_assigned = {}
		for key in self._assigned():
			if key not in kwargs:
				continue
			if kwargs[key] is None:
				continue
			if kwargs[key] != self._data[key]:
				if self._data[key] is not None:
					if key == 'samples':
						reset_assigned[key] = ["N: %d" % (_get_n_samples(self._data[key])),
											   "N: %d" % (_get_n_samples(kwargs[key]))]
					else:
						reset_assigned[key] = [self._data[key], kwargs[key]]
				self._data[key] = kwargs[key]
		if reset_assigned:
			print("\nThe following model parameters have changed:")
			for key in self._assigned():
				if key not in reset_assigned:
					continue
				old_val, new_val = reset_assigned[key]
				print("\t%s: %s -> %s" % (key, old_val, new_val))
			print()
		
		reset_assigned_ = set(reset_assigned.keys())
		
		nulled_data = self._calculated()
		reset_calculated = set()
		if reset_assigned_.intersection(assigned_full):
			# Reset all calculated attributes
			for key in nulled_data:
				reset_calculated.add(key)
				self._data[key] = nulled_data[key]
		elif reset_assigned_.intersection(assigned_randomization):
			# Reset calculated attributes associated with randomization and clustering
			for key in self._calculated_randomization() + self._calculated_clustering():
				reset_calculated.add(key)
				self._data[key] = nulled_data[key]
		elif reset_assigned_.intersection(assigned_clustering):
			# Reset calculated attributes associated with clustering
			for key in self._calculated_clustering():
				reset_calculated.add(key)
				self._data[key] = nulled_data[key]
		
		if 'curve_name' in reset_calculated:
			for name in self.samples:
				self.samples[name].calibrate(self.curve)
		
		if reset_calculated:
			print("\nThe following model attributes have been reset:")
			for key in self._calculated():
				if key in reset_calculated:
					print("\t%s" % (key))
			print()
		return reset_assigned, reset_calculated
	
	def process_phasing(self, by_clusters: bool = False, by_dates: bool = False) -> bool:
		"""
		Updates the phasing of samples based on stratigraphic relations.

		Updates the groups and phases of samples based on their stratigraphic relations.

		:param by_clusters: If True, update the phasing by clustering sample dates. Defaults to False.
		:type by_clusters: bool, optional
		:param by_dates: If True, update the phasing by comparing sample dates. Defaults to False.
		:type by_dates: bool, optional
		:return: True if phasing has changed, False otherwise.
		:rtype: bool
		"""
		
		earlier_than, samples = create_earlier_than_matrix(self)
		if by_clusters and self.is_clustered:
			earlier_than = update_earlier_than_by_clustering(self, earlier_than, samples)
		
		if by_dates and self.is_modeled:
			earlier_than = update_earlier_than_by_dating(self, earlier_than, samples)
		
		groups_phases = get_groups_and_phases(earlier_than, samples)
		# groups_phases = {sample: [group, phase], ...}
		for name in self.samples:
			if name in groups_phases:
				group, phase = groups_phases[name]
				self.samples[name].set_group(group)
				self.samples[name].set_phase(phase)
			else:
				self.samples[name].set_group(None)
				self.samples[name].set_phase(None)
		
		earlier_than_1, samples_1 = create_earlier_than_matrix(self)
		if not ((earlier_than.shape == earlier_than.shape) and np.allclose(earlier_than, earlier_than_1) and (
				samples == samples_1)):
			return True
		return False
	
	def process_outliers(self) -> None:
		"""
		Identifies and marks dating outliers among the samples in the model.

		Finds dating outliers among the samples which need to be removed for the model to be valid.
		The outliers are identified based on conflicts between their dating ranges and stratigraphic relationships with other samples.
		The identified outliers can be retrieved via the attributes Model.outliers and Model.outlier_candidates.

		:return: None
		"""
		
		outliers, self._data['outlier_candidates'] = find_dating_outliers(self)
		for name in outliers:
			self.samples[name].set_outlier(True)
	
	def process_dates(self) -> None:
		"""
		Calculates the posterior probabilities of sample dates based on phasing using Bayesian modeling in OxCal.

		Generates an OxCal file from the current model and runs the OxCal software on it.
		The results of the OxCal modeling are then loaded back into the model.
		The calculated posterior probabilities can be retrieved via the samples' attributes.

		Note: This method resets all calculated attributes that depend on the dating posteriors.

		:return: None
		"""
		
		self.to_oxcal()
		r = subprocess.call("OxCal\\bin\\OxCalWin.exe %s" % (os.path.join(self.directory, "model.oxcal")))
		self.load_oxcal_data()
		# Reset calculated attributes depending on dating posteriors
		nulled_data = self._calculated()
		for key in self._calculated_randomization() + self._calculated_clustering():
			self._data[key] = nulled_data[key]
	
	def process_randomization(self, max_cpus: int = -1, max_queue_size: int = -1) -> None:
		"""
		Performs a randomization test on the sample dates.

		Tests if the sample dates represent a uniform or normal distribution in time, depending on the Model.uniform parameter.

		:param max_cpus: Maximum number of CPUs to use for parallel processing. If -1, all available CPUs are used. Defaults to -1.
		:type max_cpus: int, optional
		:param max_queue_size: Maximum queue size for parallel processing. If -1, the queue size is unlimited. Defaults to -1.
		:type max_queue_size: int, optional
		:return: None
		"""
		
		self._data['summed'], self._data['random_lower'], self._data['random_upper'], self._data[
			'random_p'] = test_distributions(self, max_cpus=max_cpus, max_queue_size=max_queue_size)
	
	def process_clustering(self, max_cpus=-1, max_queue_size=-1) -> None:
		"""
		Performs clustering on the sample dates.

		Clusters the sample dates and uses randomization testing to find the optimal clustering solution.
		
		:param max_cpus: Maximum number of CPUs to use for parallel processing. If -1, all available CPUs are used. Defaults to -1.
		:type max_cpus: int, optional
		:param max_queue_size: Maximum queue size for parallel processing. If -1, the queue size is unlimited. Defaults to -1.
		:type max_queue_size: int, optional
		:return: None
		"""
		
		self._data['clusters'], self._data['cluster_means'], self._data['cluster_sils'], self._data['cluster_ps'], \
			self._data['cluster_opt_n'] = proc_clustering(self, max_cpus=max_cpus, max_queue_size=max_queue_size)
	
	def process(self, by_clusters: bool = False, by_dates: bool = False, max_cpus: int = -1, max_queue_size: int = -1,
				save: bool = False) -> None:
		# Process the complete model
		# by_clusters: if True, update the phasing by clustering sample dates
		# by_dates: if True, update the phasing by comparing sample dates
		"""
		Processes the complete model.

		Performs the following steps:
		1. Modeling stratigraphy to determine phasing
		2. Finding outliers
		3. Bayesian modeling of C-14 dates
		4. Testing the distribution of dates for randomness
		5. Clustering temporal distributions
		
		:param by_clusters: If True, update the phasing by clustering sample dates. Defaults to False.
		:type by_clusters: bool, optional
		:param by_dates: If True, update the phasing by comparing sample dates. Defaults to False.
		:type by_dates: bool, optional
		:param max_cpus: Maximum number of CPUs to use for parallel processing. If -1, all available CPUs are used. Defaults to -1.
		:type max_cpus: int, optional
		:param max_queue_size: Maximum queue size for parallel processing. If -1, the queue size is unlimited. Defaults to -1.
		:type max_queue_size: int, optional
		:return: None
		"""
		
		if not self.is_modeled:
			print("\nModeling stratigraphy\n")
			self.process_phasing()
			print("\nFinding outliers\n")
			self.process_outliers()
			print("\nModeling C-14 dates\n")
			self.process_dates()
			if save:
				self.save(zipped=True)
		if by_dates:
			if self.process_phasing(by_dates=True):
				print("\nUpdating phasing by comparing sample dates\n")
				self.process_dates()
				if save:
					self.save(zipped=True)
		if not self.is_randomized:
			print("\nTesting the distribution of dates for randomness\n")
			self.process_randomization(max_cpus=max_cpus, max_queue_size=max_queue_size)
			if save:
				self.save(zipped=True)
		if not self.is_clustered:
			print("\nClustering temporal distributions\n")
			self.process_clustering(max_cpus=max_cpus, max_queue_size=max_queue_size)
			if save:
				self.save(zipped=True)
		if by_clusters:
			if self.process_phasing(by_clusters=True):
				print("\nUpdating phasing by clustering\n")
				self.process_dates()
				if save:
					self.save(zipped=True)
