#!/usr/bin/env python
# -*- coding: utf-8 -*-

'''
**SiteSyncro**

Site-specific chronological modeling and synchronization

Created on 10.4.2024

Author:	Peter Demján, Institute of Archaeology of the Czech Academy of Sciences <peter.demjan@gmail.com>
Home:	https://github.com/demjanp/SiteSyncro

'''

from sitesyncro.lib.fnc_process import (process)

from sitesyncro.lib.fnc_load import (create_result_path, load_data)
from sitesyncro.lib.fnc_data import (OxCalData, RandomizeData, ClusterData)
from sitesyncro.lib.fnc_plot import (plot_randomized, plot_clusters, save_results_csv)

import shutil
import os

class SiteSyncro(object):
	"""
	The SiteSyncro class is used for site-specific chronological modeling and synchronization.

	Attributes:
		oc_data (OxCalData): An instance of the OxCalData class.
		rnd_data (RandomizeData): An instance of the RandomizeData class.
		clu_data (ClusterData): An instance of the ClusterData class.
		oc_clu_data (OxCalData): Another instance of the OxCalData class.
	"""

	def __init__(self, **kwargs):
		
		# Initialize the _data attribute with default values
		self._data = dict(
			samples = [],  # List of samples; [sample, ...]
			context_samples = {},  # Dictionary mapping contexts to samples; {context: [sample, ...], ...}
			context_area = {},  # Dictionary mapping contexts to areas; {context: area, ...}
			areas = [],  # List of areas; [area, ...]
			long_lived = {},  # Dictionary mapping samples to boolean values indicating whether they are long-lived; {sample: True/False, ...}
			r_dates = {},  # Dictionary mapping samples to tuples of age and uncertainty; {sample: (age, uncertainty), ...}
			context_phase = {},  # Dictionary mapping contexts to phases; {context: phase, ...}
			earlier_than_rel = {},  # Dictionary mapping samples to lists of other samples that they are earlier than; {sample: [sample, ...], ...}
			input = "",  # Input file path
			result = "result",  # Result directory name
			existing = False,  # Boolean indicating whether to use existing results
			curve = "intcal20.14c",  # Name of the calibration curve to use
			model = "sequence",  # Name of the model to use. Can be 'sequence', 'contiguous', 'overlapping', or 'none'.
			n = -1,  # Number of clusters (-1 = automatic)
			uniform = False,  # Boolean indicating whether to use uniform randomization
			p_value = 0.05,  # P-value for statistical tests
			uncert_base = 15,  # Base uncertainty for randomization
			npass = 100,  # Number of passes for the clustering algorithm
			convergence = 0.99,  # Convergence criterion for the clustering algorithm
			max_cpus = -1,  # Maximum number of CPUs to use
			max_queue_size = 100  # Maximum size of the queue for parallel processing
		)

		# Initialize the other attributes
		self.oc_data = OxCalData()
		self.rnd_data = RandomizeData()
		self.clu_data = ClusterData()
		self.oc_clu_data = OxCalData()

		# Update the _data attribute with the provided keyword arguments
		self._data.update(kwargs)

	def __getitem__(self, key):
		"""
		Overloads the [] operator for getting values from the _data attribute.

		Parameters:
			key (str): The key to get the value for.

		Returns:
			The value associated with the provided key.
		"""
		return self._data[key]

	def __setitem__(self, key, value):
		"""
		Overloads the [] operator for setting values in the _data attribute.

		Parameters:
			key (str): The key to set the value for.
			value: The value to set.
		"""
		self._data[key] = value

	def update_data(self, data):
		"""
		Updates the _data attribute with the provided data.

		Parameters:
			data (dict): The data to update the _data attribute with.
		"""
		self._data.update(data)

	def load_data(self, **kwargs):
		"""
		Loads data from the input file and updates the _data attribute.

		Parameters:
			kwargs (dict): A dictionary of keyword arguments that are used to update the _data attribute before loading the data.
		"""
		# Update the _data attribute with the provided keyword arguments
		for key in kwargs:
			if key in self._data:
				self._data[key] = kwargs[key]

		# Check if the input file exists
		if not os.path.isfile(self._data['input']):
			raise ValueError("Input file %s not found" % self._data['input'])

		# Load the input data
		self._data['samples'], self._data['context_samples'], \
		self._data['context_area'], self._data['areas'], \
		self._data['long_lived'], self._data['r_dates'], \
		self._data['context_phase'], self._data['earlier_than_rel'] = load_data(self._data['input'])

		# Load the models
		fname = os.path.join(self._data['result'], "model.json")
		if os.path.isfile(fname):
			self.oc_data.load_json(fname)

		fname = os.path.join(self._data['result'], "randomized.json.gz")
		if os.path.isfile(fname):
			self.rnd_data.load_json(fname)

		fname = os.path.join(self._data['result'], "clusters.json")
		if os.path.isfile(fname):
			self.clu_data.load_json(fname)

		fname = os.path.join(self._data['result'], "clu_model.json")
		if os.path.isfile(fname):
			self.oc_clu_data.load_json(fname)

	def has_data(self):
		"""
		Checks if there is any data loaded.

		Returns:
			True if there is data loaded, False otherwise.
		"""
		return (len(self._data['samples']) > 0)

	def process(self, **kwargs):
		"""
		Processes the data.
		"""
		# Update the _data attribute with the provided keyword arguments
		for key in kwargs:
			if key in self._data:
				self._data[key] = kwargs[key]

		# Create the result directory if applicable
		self._data['result'] = create_result_path(self._data['result'], self._data['existing'])

		# Copy the input file to the result directory
		shutil.copy(self._data['input'], os.path.join(self._data['result'], "data.csv"))

		# Check if there is data loaded
		if not self.has_data():
			raise Exception("Data missing.")

		# Process the data
		self.oc_data, self.rnd_data, self.clu_data, self.oc_clu_data = process(
			self._data['samples'], self._data['context_samples'], self._data['context_area'],
			self._data['areas'], self._data['context_phase'], self._data['earlier_than_rel'],
			self._data['long_lived'], self._data['r_dates'], self._data['model'],
			curve_name = self._data['curve'], result_path = self._data['result'],
			existing = self._data['existing'], clusters_n = self._data['n'],
			p_value = self._data['p_value'], uncert_base = self._data['uncert_base'], uniform = self._data['uniform'],
			npass = self._data['npass'], convergence = self._data['convergence'],
			max_cpus = self._data['max_cpus'], max_queue_size = self._data['max_queue_size'])

	def plot_randomized(self, show = False):
		"""
		Plots the randomized data.

		Parameters:
			show (bool): If True, the plot is shown. If False, the plot is saved to a file.
		"""
		# Check if there is randomized data
		if not self.rnd_data.has_data():
			print("\nNo randomization data to plot")
			return

		# Determine the file path for the plot
		fplot = None
		if not show:
			fplot = os.path.join(self._data['result'], "randomized.pdf")

		# Plot the randomized data
		plot_randomized(self.rnd_data, fplot)

	def plot_clusters(self, show = False):
		"""
		Plots the clustering data.

		Parameters:
			show (bool): If True, the plot is shown. If False, the plot is saved to a file.
		"""
		# Check if there is any clustering data
		if not self.clu_data.has_data():
			print("\nNo clustering data to plot")
			return

		# Determine the file path for the plot
		fplot = None
		if not show:
			fplot = os.path.join(self._data['result'], "silhouette.pdf")

		# Plot the clustering data
		plot_clusters(self.clu_data, fplot)

	def save_results_csv(self, fcsv = None):
		"""
		Saves the results to a CSV file.

		Parameters:
			fcsv (str): The file path for the CSV file. If None, a default file path is used.
		"""
		# Check if there are any results
		if not self.oc_data.has_data():
			print("\nNo results available")
			return

		# Determine the file path for the CSV file
		if fcsv is None:
			fcsv = os.path.join(self._data['result'], "results.csv")

		# Save the results to the CSV file
		save_results_csv(self.oc_data, self.oc_clu_data, fcsv)
