#!/usr/bin/env python
# -*- coding: utf-8 -*-

'''
**SiteSyncro**

Site-specific chronological modeling and synchronization

Created on 10.4.2024

Author:	Peter Demján, Institute of Archaeology of the Czech Academy of Sciences <peter.demjan@gmail.com>
Home:	https://github.com/demjanp/SiteSyncro

'''

from sitesyncro.lib.fnc_data import (OxCalData, RandomizeData, ClusterData)
from sitesyncro.lib.fnc_plot import (plot_randomized, plot_clusters, save_results_csv)

import argparse
import sys
import os

def parse_arguments(args):
	
	parser = argparse.ArgumentParser(description="SiteSyncro - Ploting", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument('-input', type=str, default="result\\model.json", required=False,
	                    help="File path to load model data in JSON format")
	parsed_args = parser.parse_args(args)
	
	return vars(parsed_args)  # Directly return parsed arguments as a dictionary

if __name__ == '__main__':
	
	arguments = parse_arguments(sys.argv[1:])
	input_path = arguments["input"]
	
	
	if not os.path.isfile(input_path):
		raise ValueError("Input file not found")
	
	result_path = os.path.dirname(input_path)
	
	oc_data = OxCalData()
	oc_data.load_json(input_path)
	
	# Plot randomization
	frandomized = os.path.join(result_path, "randomized.json.gz")
	if os.path.isfile(frandomized):
		rnd_data = RandomizeData()
		rnd_data.load_json(frandomized)
		plot_randomized(oc_data, rnd_data)
	
	# Plot clustering
	fclusters = os.path.join(result_path, "clusters.json")
	clu_data = None
	if os.path.isfile(fclusters):
		clu_data = ClusterData()
		clu_data.load_json(fclusters)
		plot_clusters(oc_data, clu_data)
	
	fclu_model = os.path.join(result_path, "clu_model.json")
	clu_model = None
	if os.path.isfile(fclu_model):
		clu_model = OxCalData()
		clu_model.load_json(fclu_model)
	
	save_results_csv(oc_data, clu_model, result_path)
