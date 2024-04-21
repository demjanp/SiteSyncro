#!/usr/bin/env python
# -*- coding: utf-8 -*-

'''
**SiteSync**

Site-specific chronological modeling and synchronization

Created on 10.4.2024

Author:	Peter Demján, Institute of Archaeology of the Czech Academy of Sciences <peter.demjan@gmail.com>
Home:	https://github.com/demjanp/sitesync

'''

from sitesync import SiteSync

import multiprocessing
import argparse
import sys

def parse_arguments(args):
	
	parser = argparse.ArgumentParser(description="SiteSync - Modeling", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	
	parser.add_argument('input', type=str,
		help="File path to load data in semicolon-separated CSV format (8 columns; see README.md)")
	parser.add_argument('-result', type=str, default="result", required=False, 
		help="Directory path to store the results")
	parser.add_argument('-existing', type=int, default=0, required=False,
		help="Flag indicating whether to use existing results")
	parser.add_argument('-curve', type=str, default="intcal20.14c", required=False,
		help="File name of the radiocarbon age calibration curve (see OxCal/bin directory)")
	parser.add_argument('-model', type=str, default="sequence", required=False, 
		help="OxCal model type (can be 'sequence', 'contiguous', 'overlapping', or 'none')")
	parser.add_argument('-n', type=int, default=-1, required=False,
		help="Number of clusters to form (-1 = automatic)")
	
	parser.add_argument('-uniform', type=int, default=0, required=False,
		help="Flag indicating whether to use a uniform distribution for the calendar ages")
	parser.add_argument('-p_value', type=float, default=0.05, required=False,
		help="P-value for the randomization test")
	parser.add_argument('-uncert_base', type=float, default=15, required=False,
		help="Base uncertainty for the radiocarbon dates")
	parser.add_argument('-npass', type=int, default=100, required=False,
		help="Minimum number of passes for the randomization test")
	parser.add_argument('-convergence', type=float, default=0.99, required=False,
		help="Convergence threshold for the randomization test")
	parser.add_argument('-max_cpus', type=int, default=-1, required=False,
		help="Maximum number of CPUs to use for parallel processing (-1 = all available)")
	parser.add_argument('-max_queue_size', type=int, default=100, required=False,
		help="Maximum queue size for parallel processing")
	
	parsed_args = parser.parse_args(args)
	return vars(parsed_args)  # Directly return parsed arguments as a dictionary

if __name__ == '__main__':
	multiprocessing.freeze_support()  # Needed for PyInstaller
	
	arguments = parse_arguments(sys.argv[1:])
	
	sitesync = SiteSync(**arguments)
	sitesync.load_data()
	
	# Report on input data
	print("Total unique samples:", len(sitesync['samples']))
	print("Total unique contexts:", len(sitesync['context_samples']))
	print("Total unique areas:", len(sitesync['areas']))
	print("Total earlier-than relations:", sum([len(sitesync['earlier_than_rel'][sample]) for sample in sitesync['earlier_than_rel']]))
	print()
	# Report on modeling parameters
	print("Uniform distribution:", sitesync['uniform'])
	print("P-value:", sitesync['p_value'])
	print("Number of clusters:", "all possible" if sitesync['n'] == -1 else sitesync['n'])
	print("Base uncertainty:", sitesync['uncert_base'])
	print()
	# Report on MCMC settings
	print("Minimum number of passes:", sitesync['npass'])
	print("Convergence threshold:", sitesync['convergence'])
	print("Maximum number of CPUs:", "all available" if sitesync['max_cpus'] == -1 else sitesync['max_cpus'])
	print("Maximum queue size:", sitesync['max_queue_size'])
	print()
	
	sitesync.process()
	
	sitesync.plot_randomized()
	sitesync.plot_clusters()
	sitesync.save_results_csv()
	