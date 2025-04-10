from sitesyncro import Model
from sitesyncro.utils.fnc_stat import calc_sum, calc_range

import numpy as np
import codecs

def save_ranges(data, fname):
	
	with codecs.open(fname, "w", encoding="utf-8-sig") as file:
		file.write("Name;From;To\n")
		for name, d1, d2 in data:
			file.write("%s;%f;%f\n" % (name, d1, d2))

if __name__ == '__main__':
	
	model = Model(directory="model_kap/stage_2")
	data = []
	for gr, ph in sorted(model.phases.keys()):
		phase = model.phases[(gr, ph)]
		dist_sum = calc_sum([model.samples[sample_id].posterior for sample_id in phase.samples])
		ph_from, ph_to = calc_range(phase.years, dist_sum)
		data.append([str(ph), ph_from - 1950, ph_to - 1950])
	save_ranges(data, "phase_ranges.csv")
	
	'''
	model = Model(directory="model_kap/stage_2")
	clusters = model.clusters[5]
	data = []
	for clu in [5,1,3,2,4]:
		dist_sum = calc_sum([model.samples[sample_id].posterior for sample_id in clusters[clu]])
		clu_from, clu_to = calc_range(model.years, dist_sum)
		data.append([str(clu), clu_from - 1950, clu_to - 1950])
	save_ranges(data, "cluster_ranges.csv")
	'''
	
