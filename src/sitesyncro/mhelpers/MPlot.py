import codecs
import numpy as np
from matplotlib import pyplot

class MPlot(object):
	"""
	A class implementing results plotting functionality
	
	:param model: The parent Model object
	:type model: Model
	"""
	def __init__(self, model: 'Model'):
		
		self.model = model
	
	def plot_randomized(self, fname: str, show: bool = False) -> None:
		"""
		Plots the randomization test results.

		Parameters:
		fname (str): The filename where the plot will be saved.
		show (bool): If True, the plot will be displayed. Default is False.

		Returns:
		None
		"""
		
		if self.model.random_p < self.model.p_value:
			null_hypothesis_txt = "Dates are not %s distributed." % ("uniformly" if self.model.uniform else "normally")
		else:
			null_hypothesis_txt = "Dates are %s distributed." % ("uniformly" if self.model.uniform else "normally")
		
		perc_lower = (self.model.p_value * 100) / 2
		perc_upper = 100 - perc_lower
		
		fig = pyplot.figure(figsize=(15, 4))
		pyplot.fill_between(self.model.years - 1950, self.model.random_lower, self.model.random_upper, color="lightgrey",
							label="%0.2f%% of randomized dates" % (perc_upper - perc_lower))
		pyplot.plot(self.model.years - 1950, self.model.summed, color="k", label="Observed dates")
		idxs = np.where(self.model.random_upper > 0)[0]
		idx1, idx2 = idxs.min(), idxs.max()
		pyplot.xlim(self.model.years[int(idx1)] - 1950, self.model.years[int(idx2)] - 1950)
		pyplot.gca().invert_xaxis()
		pyplot.xlabel("Calendar age (yrs BC)")
		pyplot.ylabel("Summed p")
		pyplot.annotate("p: %0.5f\n%s" % (self.model.random_p, null_hypothesis_txt), xy=(0.05, 0.95), xycoords="axes fraction",
						fontsize=12, horizontalalignment="left", verticalalignment="top")
		pyplot.legend()
		pyplot.tight_layout()
		pyplot.savefig(fname)
		if show:
			pyplot.show()
		pyplot.close()
	
	def plot_clusters(self, fname: str, show: bool = False) -> None:
		"""
		Plots Silhouette and p-value for solutions with different numbers of clusters.

		Parameters:
		fname (str): The filename where the plot will be saved.
		show (bool): If True, the plot will be displayed. Default is False.

		Returns:
		None
		"""
		
		clu_ns = np.array(sorted(list(self.model.clusters.keys())), dtype=int)
		ps_plot = np.array([self.model.cluster_ps[clu_n] for clu_n in clu_ns])
		sils_plot = np.array([self.model.cluster_sils[clu_n] for clu_n in clu_ns])
		
		color1 = "#414487"
		color2 = "#7AD151"
		color3 = "#FDE725"
		
		fig, ax1 = pyplot.subplots()
		
		ax1.set_xlabel("Clusters")
		ax1.set_ylabel("Mean Silhouette Coefficient", color=color1)
		ax1.plot(clu_ns, sils_plot, color=color1)
		ax1.plot(clu_ns, sils_plot, ".", color=color1)
		ax1.tick_params(axis="y", labelcolor=color1)
		
		ax2 = ax1.twinx()
		ax2.set_ylabel("p", color=color2)
		ax2.plot(clu_ns, ps_plot, color=color2)
		ax2.plot(clu_ns, ps_plot, ".", color=color2)
		if self.model.cluster_opt_n is not None:
			ax2.plot([self.model.cluster_opt_n, self.model.cluster_opt_n], [0, max(ps_plot.max(), sils_plot.max())], color=color3,
					 linewidth=0.7, label="Optimal no. of clusters")
		if self.model.p_value is not None:
			ax2.plot([clu_ns[0], clu_ns[-1]], [self.model.p_value, self.model.p_value], "--", color=color2, linewidth=0.7)
			ax2.annotate("p = %0.3f" % (self.model.p_value), xy=(clu_ns.mean(), self.model.p_value), xytext=(0, -3),
						 textcoords="offset pixels", va="top", ha="center", color=color2)
		ax2.tick_params(axis="y", labelcolor=color2)
		
		pyplot.xticks(clu_ns, clu_ns)
		pyplot.legend()
		
		fig.tight_layout()
		pyplot.savefig(fname)
		if show:
			pyplot.show()
		pyplot.close()
	
	def save_outliers(self, fname: str) -> None:
		"""
		Saves the outliers and outlier candidates to a file.

		Parameters:
		fname (str): The filename where the outliers will be saved.

		Returns:
		None
		"""
		txt = "Eliminated outliers:\n"
		outliers = self.model.outliers
		if outliers:
			txt += "%s\n" % (", ".join(outliers))
		else:
			txt += "None\n"
		txt += "\nOutlier candidates:\n"
		candidates = self.model.outlier_candidates
		if candidates:
			txt += "%s\n" % (", ".join(candidates))
		else:
			txt += "None\n"
		
		with open(fname, "w", encoding="utf-8") as file:
			file.write(txt)
	
	def save_results_csv(self, fcsv: str) -> None:
		"""
		Saves the results to a CSV file.

		Parameters:
		fcsv (str): The filename where the results will be saved.

		Returns:
		None
		"""
		
		def _format_year(value):
			
			if value is None:
				return "None"
			return "%0.2f" % (-(value - 1950))
		
		samples = list(self.model.samples.keys())
		
		def _sum_range(rng):
			if None in rng:
				return -1
			return sum(rng)
		
		samples = sorted(samples, key=lambda name: [
			self.model.samples[name].group,
			self.model.samples[name].phase,
			_sum_range(self.model.samples[name].likelihood_range)
		])
		
		cluster = dict([(name, None) for name in samples])
		if self.model.is_clustered:
			m_clusters = self.model.clusters[self.model.cluster_opt_n]
			for clu in m_clusters:
				for name in m_clusters[clu]:
					cluster[name] = clu
		
		with codecs.open(fcsv, "w", encoding="utf-8-sig") as file:
			file.write(
				"Name;Context;Area;C-14 Date;C-14 Uncertainty;Long-Lived;Redeposited;Outlier;EAP;Group;Phase;Cluster;Unmodeled From (CE);Unmodeled To (CE);Modeled From (CE);Modeled To (CE)\n")
			for name in samples:
				likelihood_min, likelihood_max = self.model.samples[name].likelihood_range
				posterior_min, posterior_max = self.model.samples[name].posterior_range
				file.write('''"%s";"%s";"%s";%0.2f;%0.2f;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s\n''' % (
					name, self.model.samples[name].context, self.model.samples[name].area,
					self.model.samples[name].age, self.model.samples[name].uncertainty,
					int(self.model.samples[name].long_lived), int(self.model.samples[name].redeposited),
					int(self.model.samples[name].outlier), self.model.samples[name].excavation_area_phase,
					self.model.samples[name].group, self.model.samples[name].phase, cluster[name],
					_format_year(likelihood_min), _format_year(likelihood_max),
					_format_year(posterior_min), _format_year(posterior_max)
				))

