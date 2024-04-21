from matplotlib import pyplot
import numpy as np
import os

def plot_randomized(rnd_data, fplot = None):
	
	years = rnd_data.get_years()
	sum_obs = rnd_data.get_sum_obs()
	uniform = rnd_data.get_uniform()
	p = rnd_data.get_p()
	p_value = rnd_data.get_p_value()
	sums_rnd_lower = rnd_data.get_sums_rnd_lower()
	sums_rnd_upper = rnd_data.get_sums_rnd_upper()
	
	if p < p_value:
		null_hypothesis_txt = "Dates are not %s distributed." % ("uniformly" if uniform else "normally")
	else:
		null_hypothesis_txt = "Dates are %s distributed." % ("uniformly" if uniform else "normally")
	
	perc_lower = (p_value * 100) / 2
	perc_upper = 100 - perc_lower
	
	fig = pyplot.figure(figsize = (15, 4))
	pyplot.fill_between(years - 1950, sums_rnd_lower, sums_rnd_upper, color = "lightgreen", label = "%0.2f%% of randomized dates" % (perc_upper - perc_lower))
	pyplot.plot(years - 1950, sum_obs, color = "k", label = "Observed dates")
	idxs = np.where(sums_rnd_upper > 0)[0]
	idx1, idx2 = idxs.min(), idxs.max()
	pyplot.xlim(years[int(idx1)] - 1950, years[int(idx2)] - 1950)
	pyplot.gca().invert_xaxis()
	pyplot.xlabel("Calendar age (yrs BC)")
	pyplot.ylabel("Summed p")
	pyplot.annotate("p: %0.5f\n%s" % (p, null_hypothesis_txt), xy = (0.05, 0.95), xycoords = "axes fraction", fontsize = 12, horizontalalignment = "left", verticalalignment = "top")
	pyplot.legend()
	pyplot.tight_layout()
	if fplot is None:
		pyplot.show()
	else:
		pyplot.savefig(fplot)
	pyplot.close()

def plot_clusters(clu_data, fplot = None):
	
	clusters = clu_data.get_clusters()
	means = clu_data.get_means()
	sils = clu_data.get_sils()
	ps = clu_data.get_ps()
	p_value = clu_data.get_p_value()
	opt_n = clu_data.get_opt_n()
	
	# Plot Silhouette and p-value for solutions with different numbers of clusters
	
	clu_ns = np.array(sorted(list(clusters.keys())), dtype = int)
	ps_plot = np.array([ps[clu_n] for clu_n in clu_ns])
	sils_plot = np.array([sils[clu_n] for clu_n in clu_ns])
	
	color1 = "blue"
	color2 = "green"
	
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
	if opt_n is not None:
		ax2.plot([opt_n, opt_n], [0, max(ps_plot.max(), sils_plot.max())], color = "red", linewidth = 0.7, label = "Optimal no. of clusters")
	if p_value is not None:
		ax2.plot([clu_ns[0], clu_ns[-1]], [p_value, p_value], "--", color = color2, linewidth = 0.7)
		ax2.annotate("p = %0.3f" % (p_value), xy = (clu_ns.mean(), p_value), xytext = (0, -3), textcoords = "offset pixels", va = "top", ha = "center", color = color2)
	ax2.tick_params(axis = "y", labelcolor = color2)
	
	pyplot.xticks(clu_ns, clu_ns)
	pyplot.legend()
	
	fig.tight_layout()
	if fplot is None:
		pyplot.show()
	else:
		pyplot.savefig(fplot)
	pyplot.close()

def save_results_csv(oc_data, oc_clu_data, fcsv):
	
	samples = oc_data.get_samples()
	
	r_dates = oc_data.get_r_dates()
	# r_dates = {sample: (age, uncertainty), ...}
	
	context_samples = oc_data.get_context_samples()
	# context_samples = {context: [sample, ...], ...}
	sample_contexts = {}
	for context in context_samples:
		for sample in context_samples[context]:
			sample_contexts[sample] = context
	
	context_area = oc_data.get_context_area()
	# context_area = {context: area, ...}
	
	phases = oc_data.get_phases()
	# phases = {group: {sample: phase, ...}, ...}
	phases_groups = {}
	for group in phases:
		for sample in phases[group]:
			phases_groups[sample] = (group, phases[group][sample])
	
	long_lived = oc_data.get_long_lived()
	# long_lived = {sample: True/False, ...}
	
	likelihoods = oc_data.get_likelihoods()
	ranges_prior = {}
	for sample in likelihoods:
		ranges_prior[sample] = likelihoods[sample]['range']
	# ranges_prior = {sample: [min, max], ...}
	
	if not oc_clu_data.has_data():
		posteriors = oc_data.get_posteriors()
	else:
		posteriors = oc_clu_data.get_posteriors()
	
	ranges = {}
	for sample in posteriors:
		ranges[sample] = posteriors[sample]['range']
	# ranges = {sample: [min, max], ...}
	
	phases_clustered = {}
	if not oc_clu_data.has_data():
		samples = sorted(samples, key = lambda sample: sum(ranges[sample]))
	else:
		phases_clustered = oc_clu_data.get_phases()[1]
		# phases_clustered = {sample: phase, ...}
		samples = sorted(samples, key = lambda sample: [phases_clustered[sample], sum(ranges[sample])])
	
	# Save results to a CSV file with the follwoing columns:
	# sample, context, area, age, uncertainty, long_lived, group, phase, phase_clustered, range_min, range_max
	with open(fcsv, "w") as file:
		file.write("Name;Context;Area;C-14 Date;C-14 Uncertainty;Long-Lived;Group (Stratigraphic);Phase (Stratigraphic);Phase (C-14 Clustering);Unmodelled From (CE);Unmodelled To (CE);Modelled From (CE);Modelled To (CE)\n")
		for sample in samples:
			context = sample_contexts[sample]
			area = context_area[context]
			age, uncertainty = r_dates[sample]
			long_lived_sample = int(long_lived[sample])
			group, phase = phases_groups[sample]
			phase_clustered = phases_clustered.get(sample, -1)
			range_min, range_max = ranges[sample]
			range_prior_min, range_prior_max = ranges_prior[sample]
			file.write('''"%s";"%s";"%s";%0.2f;%0.2f;%d;%d;%d;%d;%0.2f;%0.2f;%0.2f;%0.2f\n''' % (sample, context, area, age, uncertainty, long_lived_sample, group, phase, phase_clustered, range_prior_min, range_prior_max, range_min, range_max))

