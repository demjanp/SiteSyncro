from sitesyncro.utils.fnc_oxcal import (download_oxcal, gen_oxcal_model, load_oxcal_data, get_distributions)
from sitesyncro.utils.fnc_radiocarbon import (get_curve)
from sitesyncro.utils.fnc_data import (dict_keys_to_int, dict_np_to_list)
from sitesyncro.utils.fnc_load import (load_data)
from sitesyncro.utils.fnc_plot import (plot_randomized, plot_clusters, save_results_csv, save_outliers)
from sitesyncro.utils.fnc_phase import (create_earlier_than_matrix, get_groups_and_phases, find_dating_outliers, update_earlier_than_by_clustering)
from sitesyncro.utils.fnc_simulate import (test_distributions)
from sitesyncro.utils.fnc_cluster import (proc_clustering)
from sitesyncro.Sample import Sample

import os
import copy
import gzip
import json
import shutil
import subprocess
import numpy as np
from collections import defaultdict

class Model(object):
	
	def __init__(self, 
			directory = 'model',
			samples = [],
			curve_name = 'intcal20.14c',
			phase_model = 'sequence',
			cluster_n = -1,
			uniform = False,
			p_value = 0.05,
			uncertainty_base = 15,
			npass = 100,
			convergence = 0.99,
			oxcal_url = 'https://c14.arch.ox.ac.uk/OxCalDistribution.zip',
			overwrite = False,
		):
		
		for sample in samples:
			if not isinstance(sample, Sample):
				raise Exception("Invalid sample format: %s. Sample expected." % (type(sample).__name__))
		
		if phase_model not in ['sequence', 'contiguous', 'overlapping', 'none']:
			raise Exception("Invalid phase model specified: %s" % (phase_model))
		
		if not download_oxcal(url = oxcal_url):
			raise Exception("OxCal not found")
		
		self._data = self._assigned()
		self._data.update(self._calculated())
		
		self._data.update(dict(
			directory = directory,
			samples = dict([(sample.name, sample) for sample in samples]),
			curve_name = curve_name,
			phase_model = phase_model,
			cluster_n = cluster_n,
			uniform = uniform,
			p_value = p_value,
			uncertainty_base = uncertainty_base,
			npass = npass,
			convergence = convergence,
			oxcal_url = oxcal_url,
		))
		
		if overwrite and os.path.isdir(self._data['directory']):
			shutil.rmtree(self._data['directory'], ignore_errors = True)
		
		# Attempt to load model if directory exists, otherwise create it
		if not self.load():
			self._data['directory'] = self._create_dir(directory)
	
	def _assigned(self):
		return dict(
			directory = None,
			samples = {},
			curve_name = None,
			phase_model = None,
			cluster_n = None,
			uniform = None,
			p_value = None,
			uncertainty_base = None,
			npass = None,
			convergence = None,
			oxcal_url = None,
		)
	
	def _calculated(self):
		return dict(
			years = None,
			curve = None,
			oxcal_data = None,
			
			outlier_candidates = None,
			
			summed = None,
			random_p = None,
			random_lower = None,
			random_upper = None,
			
			clusters = None,
			cluster_means = None,
			cluster_sils = None,
			cluster_ps = None,
			cluster_opt_n = None,
		)
	
	def _create_dir(self, directory):
		if os.path.isdir(directory) and not os.listdir(directory):
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
	def directory(self):
		return self._data['directory']
	
	@property
	def samples(self):
		# samples = {name: Sample, ...}
		return self._data['samples']
	
	@property
	def curve_name(self):
		return self._data['curve_name']
	
	@property
	def phase_model(self):
		return self._data['phase_model']
	
	@property
	def cluster_n(self):
		return self._data['cluster_n']
	
	@property
	def uniform(self):
		return self._data['uniform']
	
	@property
	def p_value(self):
		return self._data['p_value']
	
	@property
	def uncertainty_base(self):
		return self._data['uncertainty_base']
	
	@property
	def npass(self):
		return self._data['npass']
	
	@property
	def convergence(self):
		return self._data['convergence']
	
	@property
	def oxcal_url(self):
		return self._data['oxcal_url']
	
	
	# Calculated properties
	
	@property
	def years(self):
		# Calendar years BP corresponding to the probability distributions
		# years =  np.array([year, ...])
		if self._data['years'] is None:
			if self.curve is not None:
				self._data['years'] = self.curve[:,0]
			else:
				return None
		return self._data['years'].copy()
	
	@property
	def curve(self):
		# Radiocarbon calibration curve
		# curve = np.array([[calendar year BP, C-14 year, uncertainty], ...]), sorted by calendar years
		if self._data['curve'] is None:
			if self._data['curve_name'] is None:
				return None
			else:
				self._data['curve'] = get_curve(self._data['curve_name'])
		return self._data['curve'].copy()
	
	@property
	def uncertainties(self):
		# List of uncertainties from C-14 dates of samples
		uncertainties = []
		for name in self.samples:
			if self.samples[name].date_type == 'R':
				uncertainties.append(self.samples[name].uncertainty)
		return uncertainties
	
	@property
	def oxcal_data(self):
		# OxCal data
		# oxcal_data = {key: data, ...}
		if self._data['oxcal_data'] is None:
			return None
		return copy.deepcopy(self._data['oxcal_data'])
	
	@property
	def outliers(self):
		# Returns dating outliers among samples which need to be removed for the model to be valid
		# outliers = [sample name, ...]
		return [name for name in self.samples if self.samples[name].outlier]
	
	@property
	def outlier_candidates(self):
		# Returns a candidates for outliers, from which the final outliers to be eliminated were picked
		# These samples have conflicts between dating ranges and stratigraphic relationships with other samples
		# outlier_candidates = [sample name, ...]
		if self._data['outlier_candidates'] is None:
			return []
		return copy.deepcopy(self._data['outlier_candidates'])
	
	@property
	def summed(self):
		# Summed probability distribution of the dating of all samples
		# summed = np.array([p, ...]), where p is the probability of the calendar year
		if self._data['summed'] is None:
			return None
		return self._data['summed'].copy()
	
	@property
	def random_p(self):
		# Calculated p-value for the randomization test
		return self._data['random_p']
	
	@property
	def random_lower(self):
		# Lower bound of the randomization test
		# random_lower = np.array([p, ...]), where p is the probability of the calendar year
		if self._data['random_lower'] is None:
			return None
		return self._data['random_lower'].copy()
	
	@property
	def random_upper(self):
		# Upper bound of the randomization test
		# random_upper = np.array([p, ...]), where p is the probability of the calendar year
		if self._data['random_upper'] is None:
			return None
		return self._data['random_upper'].copy()
	
	@property
	def areas(self):
		# Excavation areas that the samples belong to
		# areas = [area_name, ...]
		
		areas = set()
		for name in self.samples:
			areas.add(self.samples[name].area)
		return sorted(list(areas))
	
	@property
	def contexts(self):
		# Contexts that the samples belong to
		# contexts = [context_name, ...]
		
		contexts = set()
		for name in self.samples:
			contexts.add(self.samples[name].context)
		return sorted(list(contexts))
	
	@property
	def groups(self):
		# Groups that the samples belong to based on stratigraphic interconnection with other samples
		# groups = {group_name: [sample name, ...], ...}
		
		groups = defaultdict(list)
		for name in self.samples:
			groups[self.samples[name].group].append(name)
		return dict(groups)
	
	@property
	def clusters(self):
		# Clusters of samples based on the similarity of their probability distributions
		# clusters = {clusters_n: [cluster: [sample name, ...], ...}, ...}
		if self._data['clusters'] is None:
			return None
		return copy.deepcopy(self._data['clusters'])
	
	@property
	def cluster_means(self):
		# Mean date of the samples in each cluster in calendar years BP
		# cluster_means = {clusters_n: {cluster: year, ...}, ...}
		if self._data['cluster_means'] is None:
			return None
		return copy.deepcopy(self._data['cluster_means'])
	
	@property
	def cluster_sils(self):
		# Silhouette score of each clustering solution
		# cluster_sils = {clusters_n: silhouette, ...}
		if self._data['cluster_sils'] is None:
			return None
		return copy.deepcopy(self._data['cluster_sils'])
	
	@property
	def cluster_ps(self):
		# P-value of the clustering solutions
		# cluster_ps = {clusters_n: p, ...}
		if self._data['cluster_ps'] is None:
			return None
		return copy.deepcopy(self._data['cluster_ps'])
	
	@property
	def cluster_opt_n(self):
		# Optimal number of clusters
		return self._data['cluster_opt_n']
	
	@property
	def has_data(self):
		return (len(self.samples) > 0)
	
	@property
	def is_modeled(self):
		if not self.has_data:
			return False
		for name in self.samples:
			if not self.samples[name].is_modeled:
				return False
		return True
	
	@property
	def is_randomized(self):
		return self._data['random_p'] is not None
	
	@property
	def is_clustered(self):
		if self._data['clusters']:
			return True
		return False
	
	
	# Methods
	
	def add_sample(self, *args, **kwargs):
		# positional arguments:
		# 	name: unique name (ID) of the sample
		# 	age: C-14 age (years BP) for date_type 'R'; mean calendar age (years BP) for date_type 'U'
		# 	uncertainty: uncertainty (years BP) for date_type 'R'; 1/2 range (years BP) for date_type 'U'
		# keyword arguments:
		# 	date_type: 'R' for radiocarbon date; 'U' for calendar date as a uniform distribution
		# 	long_lived: True if sample could be older than the examined deposition event due to e.g. old wood effect
		# 	redeposited: True if sample could be redeposited from a different context
		#	outlier: True if sample is an outlier and should not be used for modeling
		# 	context: name of the context where sample was found
		# 	area: excavation area
		# 	area_excavation_phase: chronological phase of the context within the excavation area (integer, lower = earlier (older) phase)
		# 	earlier_than: list of names of samples which are stratigraphically later (younger) than this sample
		
		def _from_arguments(name, age, uncertainty, 
				date_type='R', long_lived=False, redeposited=False, outlier=False, context=None,
				area = None,
				area_excavation_phase = None,
				earlier_than = [],
			):
			self._data['samples'][name] = Sample(
				name, age, uncertainty, date_type, long_lived, redeposited, outlier,
				context, area, area_excavation_phase, earlier_than, self.curve
			)
			self._data['samples'][name].calibrate(self.curve)
		
		if len(args) == 1 and isinstance(args[0], Sample):
			sample = args[0].copy()
			if not sample.is_calibrated:
				sample.calibrate(self.curve)
			self._data['samples'][sample.name] = sample
		else:
			_from_arguments(*args, **kwargs)
	
	def del_sample(self, name):
		if name in self._data['samples']:
			del self._data['samples'][name]
	
	def reset_model(self):
		# Reset calculated properties
		self._data.update(self._calculated())
		for name in self.samples:
			self.samples[name].set_posterior(None)
	
	def save(self, zipped = False):
		# Save the model to a JSON file
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
	
	def copy(self, directory):
		# Create a copy of the model with a new directory
		samples = dict([(name, self.samples[name].copy()) for name in self.samples])
		model = Model(directory, samples, self.curve_name, self.phase_model, 
			self.cluster_n, self.uniform, self.p_value, self.uncertainty_base, self.oxcal_url
		)
		for key in self._calculated():
			model._data[key] = getattr(self, key)
		return model
	
	def load(self):
		# Load the model from a JSON file
		
		# Check if model file exists
		zipped = False
		fname = os.path.join(self.directory, 'model.json')
		if not os.path.isfile(fname):
			zipped = True
			fname = os.path.join(self.directory, 'model.json.gz')
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
		
		# Check if model parameters have changed
		changed = False
		for key in self._assigned():
			if key in ['directory','samples']:
				continue
			if data[key] != self._data[key]:
				changed = True
				break
		if changed:
			self.reset_model()
			print("\nModel parameters have changed. Resetting model.")
		
		# Initialize samples
		for name in data['samples']:
			data['samples'][name] = Sample(data['samples'][name])
			if not data['samples'][name].is_calibrated:
				data['samples'][name].calibrate(self.curve)
		
		# Convert lists to numpy arrays
		for key in ['years', 'curve', 'summed', 'random_lower', 'random_upper']:
			if data[key] is not None:
				data[key] = np.array(data[key], dtype = np.float64)
		
		oxcal_url = self.oxcal_url
		directory = self.directory
		
		self._data = data
		self._data['oxcal_url'] = oxcal_url
		self._data['directory'] = directory
		return True
	
	def import_csv(self, fname):
		# Import model input data from a CSV file
		
		# Check if the input file exists
		if not os.path.isfile(fname):
			raise ValueError("Input file %s not found" % fname)
		
		samples, contexts, context_area, long_lived, redeposited, outlier, r_dates, context_phase, earlier_than = load_data(fname)
		# samples = [sample, ...]
		# contexts = {sample: context, ...}
		# context_area = {context: area, ...}
		# long_lived = {sample: True/False, ...}
		# redeposited = {sample: True/False, ...}
		# outlier = {sample: True/False, ...}
		# r_dates = {sample: (age, uncertainty), ...}
		# context_phase = {context: phase, ...}
		# earlier_than = {sample: [sample, ...], ...}
		
		self.reset_model()
		self._data['samples'] = {}
		for name in samples:
			age, uncertainty = r_dates[name]
			self.add_sample(
				name, age, uncertainty, 'R', long_lived[name], redeposited[name], outlier[name], contexts[name], 
				context_area[contexts[name]], context_phase[contexts[name]], earlier_than[name]
			)
	
	def plot_randomized(self, show = False):
		# Plot the randomization test results
		# If `show` is True, the plot is shown. If False, the plot is saved to a file.
		
		if not self.is_randomized:
			print("\nNo randomization data to plot")
			return False
		
		# Determine the file path for the plot
		fplot = None
		if not show:
			fplot = os.path.join(self.directory, "randomized.pdf")
		
		plot_randomized(self, fplot)
		return True
	
	def plot_clusters(self, show = False):
		# Plot the clustering results
		# If `show` is True, the plot is shown. If False, the plot is saved to a file.
		
		if not self.is_clustered:
			print("\nNo clustering data to plot")
			return False
		
		# Determine the file path for the plot
		fplot = None
		if not show:
			fplot = os.path.join(self.directory, "silhouette.pdf")

		# Plot the clustering data
		plot_clusters(self, fplot)
	
	def save_csv(self, fcsv = None):
		# Saves the results to a CSV file
		# If fcsv is None, a default file path and name is used
		
		# Check if there are any results
		if not self.has_data:
			print("\nNo data available")
			return
		
		# Determine the file path for the CSV file
		if fcsv is None:
			fcsv = os.path.join(self.directory, "results.csv")
		
		# Save the results to the CSV file
		save_results_csv(self, fcsv)
	
	def save_outliers(self, fname = None):
		# Saves a list of outliers to fname
		# If fname is None, a default file path and name is used
		
		if fname is None:
			fname = os.path.join(self.directory, "outliers.txt")
		
		save_outliers(self, fname)
	
	def to_oxcal(self):
		# Export the model to an OxCal file
		
		# Check if there is any data
		if not self.has_data:
			print("\nNo data available")
			return None
		
		txt = gen_oxcal_model(self)
		fmodel = os.path.join(self.directory, "model.oxcal")
		with open(fmodel, "w") as file:
			file.write(txt)
	
	def load_oxcal_data(self):
		
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
	
	def update_params(self, **kwargs):
		if 'directory' in kwargs:
			raise Exception("Cannot update directory parameter. Use Model.copy(directory).")
		if 'samples' in kwargs:
			raise Exception("Cannot update samples parameter. Use add_sample and del_sample methods.")
		
		# Update model parameters
		for key in self._assigned():
			if (key in kwargs) and (kwargs[key] is not None):
				self._data[key] = kwargs[key]
		
		self.reset_model()
		
		if 'curve_name' in kwargs:
			for name in self.samples:
				self.samples[name].calibrate(self.curve)
	
	def process_phasing(self, by_clusters = False):
		# Update groups and phases of samples based on stratigraphic relations
		# by_clusters: if True, update the phasing by clustering sample dates
		earlier_than, samples = create_earlier_than_matrix(self)
		if by_clusters and self.is_clustered:
			earlier_than = update_earlier_than_by_clustering(self, earlier_than, samples)
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
	
	def process_outliers(self, max_cpus = -1, max_queue_size = 10000):
		outliers, self._data['outlier_candidates'] = find_dating_outliers(self)
		for name in outliers:
			self.samples[name].set_outlier(True)
		
	def process_dates(self):
		# Calculate posteriors of sample dates based on phasing using bayesian modeling in OxCal
		self.to_oxcal()
		r = subprocess.call("OxCal\\bin\\OxCalWin.exe %s" % (os.path.join(self.directory, "model.oxcal")))
		self.load_oxcal_data()
	
	def process_randomization(self, max_cpus = -1, max_queue_size = 10000):
		# Test if sample dates represent a uniform / normal (depending on Model.uniform parameter) distribution in time
		self._data['summed'], self._data['random_lower'], self._data['random_upper'], self._data['random_p'] = test_distributions(self, max_cpus = max_cpus, max_queue_size = max_queue_size)
	
	def process_clustering(self, max_cpus = -1, max_queue_size = 10000):
		# Cluster dates and using randomization testing find optimal clustering solution
		self._data['clusters'], self._data['cluster_means'], self._data['cluster_sils'], self._data['cluster_ps'], self._data['cluster_opt_n'] = proc_clustering(self, max_cpus = max_cpus, max_queue_size = max_queue_size)
	
	def process(self, by_clusters = False, max_cpus = -1, max_queue_size = 10000):
		# Process the complete model
		# by_clusters: if True, update the phasing by clustering sample dates
		self.reset_model()
		print("\nModeling stratigraphy\n")
		self.process_phasing()
		print("\nFinding outliers\n")
		self.process_outliers(max_cpus = max_cpus, max_queue_size = max_queue_size)
		print("\nModeling C-14 dates\n")
		self.process_dates()
		print("\nTesting the distribution of dates for randomness\n")
		self.process_randomization(max_cpus = max_cpus, max_queue_size = max_queue_size)
		print("\nClustering temporal distributions\n")
		self.process_clustering(max_cpus = max_cpus, max_queue_size = max_queue_size)
		if by_clusters:
			print("\nUpdating phasing by clustering\n")
			self.process_phasing(by_clusters = True)
			print("\nModeling C-14 dates\n")
			self.process_dates()

