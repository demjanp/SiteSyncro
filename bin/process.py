#!/usr/bin/env python
# -*- coding: utf-8 -*-

'''
**SiteSyncro**

Site-specific chronological modeling and synchronization

Created on 10.4.2024

Author:	Peter Demján, Institute of Archaeology of the Czech Academy of Sciences <peter.demjan@gmail.com>
Home:	https://github.com/demjanp/SiteSyncro

'''

from sitesyncro import Model
from sitesyncro import __version__

import multiprocessing
import argparse
import sys
import os

DESCRIPTION = "SiteSyncro v%s - Site-specific chronological modeling and synchronization (https://github.com/demjanp/SiteSyncro)" % (__version__)

def parse_arguments(args):
	
	parser = argparse.ArgumentParser(description=DESCRIPTION, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	
	parser.add_argument('-input', type=str, required=False,
		help="File path to load data in semicolon-separated CSV format (8 columns; see README.md)")
	parser.add_argument('-directory', type=str, default="model", required=False, 
		help="Working directory for model data")
	parser.add_argument('-curve_name', type=str, default="intcal20.14c", required=False,
		help="File name of the radiocarbon age calibration curve (see OxCal/bin directory)")
	parser.add_argument('-phase_model', type=str, default="sequence", required=False, 
		help="OxCal phase model type (can be 'sequence', 'contiguous', 'overlapping', or 'none')")
	parser.add_argument('-cluster_n', type=int, default=-1, required=False,
		help="Number of clusters to form (-1 = automatic)")
	parser.add_argument('-by_clusters', type=int, default=0, required=False,
		help="Flag indicating whether to update the phasing by clustering sample dates")
	parser.add_argument('-uniform', type=int, default=0, required=False,
		help="Flag indicating whether to use a uniform distribution for the calendar ages")
	parser.add_argument('-p_value', type=float, default=0.05, required=False,
		help="P-value for the randomization test")
	parser.add_argument('-uncertainty_base', type=float, default=15, required=False,
		help="Base uncertainty for the radiocarbon dates")
	parser.add_argument('-npass', type=int, default=100, required=False,
		help="Minimum number of passes for the randomization test")
	parser.add_argument('-convergence', type=float, default=0.99, required=False,
		help="Convergence threshold for the randomization test")
	parser.add_argument('-max_cpus', type=int, default=-1, required=False,
		help="Maximum number of CPUs to use for parallel processing (-1 = all available)")
	parser.add_argument('-max_queue_size', type=int, default=10000, required=False,
		help="Maximum queue size for parallel processing")
	
	parsed_args = parser.parse_args(args)
	return vars(parsed_args)  # Directly return parsed arguments as a dictionary

if __name__ == '__main__':
	multiprocessing.freeze_support()  # Needed for PyInstaller
	
	arguments = parse_arguments(sys.argv[1:])
	
	arguments['uniform'] = bool(arguments['uniform'])
	
	finput = arguments.pop('input', None)
	by_clusters = bool(arguments.pop('by_clusters', False))
	max_cpus = arguments.pop('max_cpus')
	max_queue_size = arguments.pop('max_queue_size')
	
	print()
	print(DESCRIPTION)
	print()
	print("use option -h or --help to show help message")
	print()
	
	model = Model(**arguments)
	if (finput is not None) and os.path.isfile(finput):
		model.import_csv(finput)	
		print("Loaded data from:", finput)
		print()
	
	print("Working directory:", model.directory)
	print()
	print("Total unique samples:", len(model.samples))
	print("Total unique contexts:", len(model.contexts))
	print("Total unique areas:", len(model.areas))
	print()
	print("Uniform distribution:", model.uniform)
	print("P-value:", model.p_value)
	print("Number of clusters:", "all possible" if model.cluster_n == -1 else model.cluster_n)
	print("Base uncertainty:", model.uncertainty_base)
	print()
	print("Minimum number of passes:", model.npass)
	print("Convergence threshold:", model.convergence)
	print()
	
	model.process(by_clusters = by_clusters, max_cpus = max_cpus, max_queue_size = max_queue_size)
	model.save(zipped = True)
	model.plot_randomized()
	model.plot_clusters()
	model.save_outliers()
	model.save_csv()
	