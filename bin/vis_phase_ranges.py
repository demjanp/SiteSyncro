from sitesyncro import Model
#from sitesyncro.utils.fnc_stat import calc_mean, calc_sum, calc_range

import matplotlib.ticker as ticker
import matplotlib.pyplot as plt
import numpy as np
import codecs

def save_ranges(data, fname):
	
	with codecs.open(fname, "w", encoding="utf-8-sig") as file:
		file.write("Name;Start From;Start To;Start Mean;End From;End To;End Mean\n")
		for row in data:
			file.write("%s;%f;%f;%f;%f;%f;%f\n" % tuple(row))

def vis_ranges(data, fname):
	# Sort data by phase descending (top to bottom in the plot)
	data = sorted(data, key=lambda x: x[0])
	
	phases = [row[0] for row in data]
	start_means = [row[3] for row in data]
	end_means = [row[6] for row in data]
	start_froms = [row[1] for row in data]
	start_tos = [row[2] for row in data]
	end_froms = [row[4] for row in data]
	end_tos = [row[5] for row in data]

	y_pos = list(range(len(data)))

	fig, ax = plt.subplots(figsize=(10, 0.6 * len(data)))
	
	# Draw horizontal bars for phase ranges
	for i, (y, sm, em) in enumerate(zip(y_pos, start_means, end_means)):
		
		label = "Range between means of Start and End" if i==0 else ""
		ax.plot([sm, em], [y, y], color='grey', linewidth=10, label=label)

		# Start error bar
		label = "Start (95.45% range)" if i==0 else ""
		ax.plot([start_froms[i], start_tos[i]], [y - 0.02, y - 0.02], color='blue', linewidth=1, label=label)
		ax.plot([start_froms[i], start_froms[i]], [y - 0.2, y + 0.2], color='blue', linewidth=1)
		ax.plot([start_tos[i], start_tos[i]], [y - 0.2, y + 0.2], color='blue', linewidth=1)

		# End error bar
		label = "End (95.45% range)" if i==0 else ""
		ax.plot([end_froms[i], end_tos[i]], [y + 0.02, y + 0.02], color='red', linewidth=1, label=label)
		ax.plot([end_froms[i], end_froms[i]], [y - 0.2, y + 0.2], color='red', linewidth=1)
		ax.plot([end_tos[i], end_tos[i]], [y - 0.2, y + 0.2], color='red', linewidth=1)
	
	# Y axis
	ax.set_yticks(y_pos)
	ax.set_yticklabels(phases)
	ax.set_ylabel("Phase")

	# X axis (Years BCE, reversed)
	ax.set_xlabel("Years BCE")
	ax.invert_xaxis()
	ax.xaxis.set_major_locator(ticker.MultipleLocator(50))
	ax.xaxis.set_minor_locator(ticker.MultipleLocator(10))
	ax.tick_params(axis='x', which='minor', length=4)
	ax.tick_params(axis='x', which='major', length=8)
	
	plt.legend()
	plt.title("Phase Dating Intervals")
	
	# Tight layout and save
	plt.tight_layout()
	plt.savefig(fname, format='pdf')
	plt.close()

ROOT = "model_kap_no_cer"

if __name__ == '__main__':
	
	model = Model(directory="%s/stage_2" % (ROOT))
	data = []
	for gr, ph in sorted(model.phases.keys()):
		phase = model.phases[(gr, ph)]
		start_from, start_to = phase.start_range
		start_mean = phase.start_mean
		end_from, end_to = phase.end_range
		end_mean = phase.end_mean
		
		start_from, start_to, start_mean, end_from, end_to, end_mean = [val - 1950 for val in [start_from, start_to, start_mean, end_from, end_to, end_mean]]
		data.append([ph, start_from, start_to, start_mean, end_from, end_to, end_mean])
	save_ranges(data, "%s/phases.csv" % (ROOT))
	vis_ranges(data, "%s/phases.pdf" % (ROOT))
