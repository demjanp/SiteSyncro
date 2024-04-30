from matplotlib import pyplot
import numpy as np

def plot_randomized(model: object, fname: str, show: bool = False):
	# Plot the randomization test results
	
	if model.random_p < model.p_value:
		null_hypothesis_txt = "Dates are not %s distributed." % ("uniformly" if model.uniform else "normally")
	else:
		null_hypothesis_txt = "Dates are %s distributed." % ("uniformly" if model.uniform else "normally")
	
	perc_lower = (model.p_value * 100) / 2
	perc_upper = 100 - perc_lower
	
	fig = pyplot.figure(figsize = (15, 4))
	pyplot.fill_between(model.years - 1950, model.random_lower, model.random_upper, color = "lightgrey", label = "%0.2f%% of randomized dates" % (perc_upper - perc_lower))
	pyplot.plot(model.years - 1950, model.summed, color = "k", label = "Observed dates")
	idxs = np.where(model.random_upper > 0)[0]
	idx1, idx2 = idxs.min(), idxs.max()
	pyplot.xlim(model.years[int(idx1)] - 1950, model.years[int(idx2)] - 1950)
	pyplot.gca().invert_xaxis()
	pyplot.xlabel("Calendar age (yrs BC)")
	pyplot.ylabel("Summed p")
	pyplot.annotate("p: %0.5f\n%s" % (model.random_p, null_hypothesis_txt), xy = (0.05, 0.95), xycoords = "axes fraction", fontsize = 12, horizontalalignment = "left", verticalalignment = "top")
	pyplot.legend()
	pyplot.tight_layout()
	pyplot.savefig(fname)
	if show:
		pyplot.show()
	pyplot.close()

def plot_clusters(model: object, fname: str, show: bool = False):
	# Plot Silhouette and p-value for solutions with different numbers of clusters
	
	clu_ns = np.array(sorted(list(model.clusters.keys())), dtype = int)
	ps_plot = np.array([model.cluster_ps[clu_n] for clu_n in clu_ns])
	sils_plot = np.array([model.cluster_sils[clu_n] for clu_n in clu_ns])
	
	color1 = "#414487"
	color2 = "#7AD151"
	color3 = "#FDE725"
	
	fig, ax1 = pyplot.subplots()
	
	ax1.set_xlabel("Clusters")
	ax1.set_ylabel("Mean Silhouette Coefficient", color = color1)
	ax1.plot(clu_ns, sils_plot, color = color1)
	ax1.plot(clu_ns, sils_plot, ".", color = color1)
	ax1.tick_params(axis = "y", labelcolor = color1)
	
	ax2 = ax1.twinx()
	ax2.set_ylabel("p", color = color2)
	ax2.plot(clu_ns, ps_plot, color = color2)
	ax2.plot(clu_ns, ps_plot, ".", color = color2)
	if model.cluster_opt_n is not None:
		ax2.plot([model.cluster_opt_n, model.cluster_opt_n], [0, max(ps_plot.max(), sils_plot.max())], color = color3, linewidth = 0.7, label = "Optimal no. of clusters")
	if model.p_value is not None:
		ax2.plot([clu_ns[0], clu_ns[-1]], [model.p_value, model.p_value], "--", color = color2, linewidth = 0.7)
		ax2.annotate("p = %0.3f" % (model.p_value), xy = (clu_ns.mean(), model.p_value), xytext = (0, -3), textcoords = "offset pixels", va = "top", ha = "center", color = color2)
	ax2.tick_params(axis = "y", labelcolor = color2)
	
	pyplot.xticks(clu_ns, clu_ns)
	pyplot.legend()
	
	fig.tight_layout()
	pyplot.savefig(fname)
	if show:
		pyplot.show()
	pyplot.close()

def save_outliers(model: object, fname: str):
	
	txt = "Eliminated outliers:\n"
	outliers = model.outliers
	if outliers:
		txt += "%s\n" % (", ".join(outliers))
	else:
		txt += "None\n"
	txt += "\nOutlier candidates:\n"
	candidates = model.outlier_candidates
	if candidates:
		txt += "%s\n" % (", ".join(candidates))
	else:
		txt += "None\n"
	
	with open(fname, "w") as file:
		file.write(txt)

def save_results_csv(model: object, fcsv: str):
	# Save results to a CSV file
	
	def _format_year(value):
		
		if value is None:
			return "None"
		return "%0.2f" % (-(value - 1950))
	
	samples = list(model.samples.keys())
	
	def _sum_range(rng):
		if None in rng:
			return -1
		return sum(rng)
	
	samples = sorted(samples, key = lambda name: [
		model.samples[name].group, 
		model.samples[name].phase, 
		_sum_range(model.samples[name].likelihood_range)
	])
	
	cluster = dict([(name, None) for name in samples])
	if model.is_clustered:
		m_clusters = model.clusters[model.cluster_opt_n]
		for clu in m_clusters:
			for name in m_clusters[clu]:
				cluster[name] = clu
	
	with open(fcsv, "w") as file:
		file.write("Name;Context;Area;C-14 Date;C-14 Uncertainty;Long-Lived;Redeposited;Outlier;Group;Phase;Cluster;Unmodeled From (CE);Unmodeled To (CE);Modeled From (CE);Modeled To (CE)\n")
		for name in samples:
			likelihood_min, likelihood_max = model.samples[name].likelihood_range
			posterior_min, posterior_max = model.samples[name].posterior_range
			file.write('''"%s";"%s";"%s";%0.2f;%0.2f;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s\n''' % (
				name, model.samples[name].context, model.samples[name].area, 
				model.samples[name].age, model.samples[name].uncertainty, 
				int(model.samples[name].long_lived), int(model.samples[name].redeposited), int(model.samples[name].outlier), 
				model.samples[name].group, model.samples[name].phase, cluster[name], 
				_format_year(likelihood_min), _format_year(likelihood_max),
				_format_year(posterior_min), _format_year(posterior_max)
			))

