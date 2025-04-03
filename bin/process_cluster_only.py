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

if __name__ == '__main__':
	
	
	finput = 'data_20250220_no_strat.csv'
	
	model = Model(
		directory = 'stage_0',
		uniform = True,
		cluster_selection = 'silhouette',
		phase_model = 'overlapping',
	)
	
	model.import_csv(finput)	
	model.save(zipped=True)
	model.save_csv_samples()
	
	print()
	print("   Total unique samples:     %s" % len(model.samples))
	print("   Total unique contexts:    %s" % len(model.contexts))
	print("   Total unique areas:       %s" % len(model.areas))
	print()
	print("Model parameters:")
	print("   Calibration curve:      %s" % model.curve_name)
	print("   Phase model:            %s" % model.phase_model)
	print("   Number of clusters:     %s" % ("all possible" if model.cluster_n == -1 else model.cluster_n))
	print("   Uniform distribution:   %s" % model.uniform)
	print("   P-value threshold:      %s" % model.p_value)
	print("   Base uncertainty:       %s" % model.uncertainty_base)
	print("   Min. number of passes:  %s" % model.npass)
	print("   Convergence threshold:  %s" % model.convergence)
	print()
	
	print("Stage 1\n")
	print("\nModelling 1\n")
	model = model.copy('stage_1')
	model.process_phasing()
	model.process_outliers()
	model.process_dates()
	model.save(zipped=True)
	model.save_outliers()
	print("\nRandomization\n")
	model.process_randomization()
	model.save(zipped=True)
	model.plot_randomized()
	print("\nClustering 1\n")
	model.process_clustering()
	model.save(zipped=True)
	model.save_csv_samples()
	model.save_csv_phases()
	model.plot_clusters()

	print("\nStage 2\n")
	model = model.copy('stage_2')
	print("\nModelling 2 - by clusters\n")
	model.process_phasing(by_clusters=True)
	model.process_dates()
	model.save(zipped=True)
	model.save_csv_samples()
	model.save_csv_phases()