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

import multiprocessing
import argparse
import sys
import os

DESCRIPTION = "SiteSyncro - Site-specific chronological modeling and synchronization (https://github.com/demjanp/SiteSyncro)"

def parse_arguments(args):
	
	parser = argparse.ArgumentParser(description=DESCRIPTION, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	
	parser.add_argument('input', type=str,
		help="File path to load data in semicolon-separated CSV format (8 columns; see README.md)")
	parser.add_argument('-directory', type=str, default="model", required=False, 
		help="Directory path to store the results")
	parser.add_argument('-curve_name', type=str, default="intcal20.14c", required=False,
		help="File name of the radiocarbon age calibration curve (see OxCal/bin directory)")
	parser.add_argument('-phase_model', type=str, default="sequence", required=False, 
		help="OxCal phase model type (can be 'sequence', 'contiguous', 'overlapping', or 'none')")
	parser.add_argument('-cluster_n', type=int, default=-1, required=False,
		help="Number of clusters to form (-1 = automatic)")
	parser.add_argument('-uniform', type=int, default=0, required=False,
		help="Flag indicating whether to use a uniform distribution for the calendar ages")
	parser.add_argument('-p_value', type=float, default=0.05, required=False,
		help="P-value for the randomization test")
	parser.add_argument('-uncertainty_base', type=float, default=15, required=False,
		help="Base uncertainty for the radiocarbon dates")
	parser.add_argument('-npass', type=int, default=10, required=False,
		help="Minimum number of passes for the randomization test")
	parser.add_argument('-convergence', type=float, default=0.9, required=False,
		help="Convergence threshold for the randomization test")
	
	parsed_args = parser.parse_args(args)
	return vars(parsed_args)  # Directly return parsed arguments as a dictionary

if __name__ == '__main__':
	multiprocessing.freeze_support()  # Needed for PyInstaller
	
	arguments = parse_arguments(sys.argv[1:])
	
	finput = arguments.pop('input')
	directory = arguments.pop('directory')
	
	print()
	print(DESCRIPTION)
	
	if not os.path.isfile(finput):
		raise Exception("Input file not found:", finput)
	
	model = Model(directory, **arguments)
	if not model.is_processed:
		model.import_csv(finput)
	
	print()
	print("Loaded data from:", finput)
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
	
#	model.process()
#	model.process_randomization()
#	model.process_clustering()
	
	model.plot_randomized()
	model.plot_clusters()
	model.save_csv()
	