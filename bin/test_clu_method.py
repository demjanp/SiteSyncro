from sitesyncro import Model
from sitesyncro.utils.fnc_radiocarbon import (get_curve)

import os
import json
import numpy as np
import matplotlib.pyplot as plt

UNCERTAINTY_BASE = 30
DURATION = 200
GAP = 100
CENTER_DATE = 3500
DATES_N = {'small': 10, 'large': 50}
ITERATIONS = 100
#METHODS = ['silhouette', 'mcst']
METHODS = ['silhouette']

FRESULTS = "test_results.json"
DIRECTORY = "test_clu"

def sim_event(cal_age_bp, curve):
	# Simulate a C-14 dated event
	
	# Find the closest index in the calibration curve
	idx = np.argmin(np.abs(curve[:, 0] - cal_age_bp))
	# Get the corresponding radiocarbon age
	return curve[idx, 1]

def sim_phased_events(cal_bp_mean, duration, gap, phases_n, events_n, normal=True):
    # Simulate events_n events in phases_n phases with gaps of `gap` years and `duration` years per phase
    # Ensures non-overlapping phases in decreasing chronological order (older to younger)

    if phases_n < 2:
        raise Exception(f"Invalid number of phases: {phases_n}. Required: >2.")
    
    # Randomly split events_n among phases
    ns = np.random.multinomial(events_n, [1 / phases_n] * phases_n)

    # Compute total range from oldest to youngest phase center
    total_range = (duration + gap) * (phases_n - 1)
    start = cal_bp_mean + total_range / 2  # oldest phase center
    means = [start - i * (duration + gap) for i in range(phases_n)]  # phase centers

    events = []
    for i, (mean, n) in enumerate(zip(means, ns)):
        if normal:
            cal_ages = np.random.normal(loc=mean, scale=duration / 4, size=n)
        else:
            cal_ages = np.random.uniform(low=mean - duration / 2, high=mean + duration / 2, size=n)

        cal_ages = sorted(cal_ages, reverse=True)
        for cal_age_bp in cal_ages:
            events.append([i + 1, cal_age_bp])
    
    return events

def sim_dates(events, uncertainty_base, curve):
	# Simulate C-14 dates based on supplied events
	# returns [[phase, c14age, uncert], ...]
	
	dates = []
	for phase, cal_age_bp in events:
		c14age = sim_event(cal_age_bp, curve)
		uncert = uncertainty_base * np.exp(c14age / (2 * 8033))
		dates.append([phase, c14age, np.random.normal(uncert, 2)])
	return dates

def vis_events(events, fname):
	"""
	Visualize events grouped by phase as box plots.
	
	Parameters:
		events (list): List of [phase_number, cal_age_bp] pairs.
	"""
	# Organize events into phases
	from collections import defaultdict
	phase_data = defaultdict(list)
	for phase, cal_age_bp in events:
		phase_data[phase].append(cal_age_bp)

	# Sort by phase number
	sorted_phases = sorted(phase_data.items())
	labels = [f"Phase {phase}" for phase, _ in sorted_phases]
	data = [cal_ages for _, cal_ages in sorted_phases]

	# Plot
	plt.figure(figsize=(10, 6))
	plt.boxplot(data, labels=labels, vert=False)
	plt.gca().invert_yaxis()  # Older (higher BP) at top
	plt.xlabel("Calibrated Age BP")
	plt.title("Simulated Phased Events")
	plt.grid(axis='x')
	plt.tight_layout()
	plt.savefig(fname)
	plt.close()

if __name__ == '__main__':
	
	curve = get_curve()
	results = {}
	for phases_n in [2, 4, 10]:
		results[phases_n] = {}
		for sample_size in DATES_N:
			results[phases_n][sample_size] = {}
			for method in METHODS:
				results[phases_n][sample_size][method] = {}
				for use_wd in [False, True]:
					results[phases_n][sample_size][method][use_wd] = []
			for i in range(ITERATIONS):
				print("\nPhases: %d, Sample size: %s, Iter. %d/%d" % (phases_n, sample_size, i + 1, ITERATIONS))
				events = sim_phased_events(CENTER_DATE, DURATION, GAP, phases_n, DATES_N[sample_size], normal = False)
				
#				vis_events(events, "%d_%s_%03d.png" % (phases_n, sample_size, i))
#				tst = np.array([evt[1] - 1950 for evt in events])  # DEBUG
#				print(tst.min(), tst.max())  # DEBUG
#				continue  # DEBUG
				
				dates = sim_dates(events, UNCERTAINTY_BASE, curve)
				
				for method in METHODS:
					for use_wd in [False, True]:
						model = Model(directory = DIRECTORY, cluster_selection = method, use_wasserstein = use_wd, uniform = True, overwrite = True)
						n = 1
						for phase, age, uncert in dates:
							name = "%d_%d" % (phase, n)
							n += 1
							model.add_sample(name, age, uncert)
						model.process_clustering(max_clusters = 2*phases_n)
						model.process_phasing(by_clusters = True)
						ph_good = 0
						for name in model.samples:
							if model.samples[name].phase == int(name.split("_")[0]):
								ph_good += 1
						phasing_check = ph_good / len(model.samples)
						results[phases_n][sample_size][method][use_wd].append([ph_good, len(model.samples), model.cluster_opt_n])
						with open(FRESULTS, 'w') as file:
							json.dump(results, file)
						print("Phasing check: %0.4f" % (phasing_check))
