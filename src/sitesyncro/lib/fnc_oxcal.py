from sitesyncro.lib.fnc_data import (OxCalData)
from sitesyncro.lib.fnc_phase import (update_earlier_than_matrix, get_groups_and_phases)

import os
import requests
import zipfile
import subprocess
from tqdm import tqdm
from collections import defaultdict

def download_oxcal(url = None):
	
	if os.path.isfile("OxCal\\bin\\OxCalWin.exe"):
		return True
	# Download and unzip OxCal from url to the current directory
	
	if url is None:
		url = "https://c14.arch.ox.ac.uk/OxCalDistribution.zip"
	
	print("Downloading OxCal from %s" % (url))
	
	# get filename from url
	local_zip = os.path.basename(url)
	
	# Download the file from `url` and save it locally under `local_zip`:
	try:
		response = requests.get(url, stream=True)
	except requests.exceptions.RequestException as e:
		print (f"Error: Unable to connect to {url}. Please check the URL or your internet connection.")
		return False
	
	total_size_in_bytes = int(response.headers.get('content-length', 0))
	progress_bar = tqdm(total=total_size_in_bytes, unit='iB', unit_scale=True)
	
	with open(local_zip, 'wb') as f:
		for chunk in response.iter_content(chunk_size=8192):
			progress_bar.update(len(chunk))
			f.write(chunk)
	
	progress_bar.close()
	
	# Unzip the file
	with zipfile.ZipFile(local_zip, 'r') as zip_ref:
		zip_ref.extractall()
	
	# Remove the zip file
	os.remove(local_zip)
	
	return True

def run_oxcal(fmodel, url = None):
	
	if not download_oxcal(url):
		raise Exception("OxCal not found")
	r = subprocess.call("OxCal\\bin\\OxCalWin.exe %s" % (fmodel))

def gen_dates(dates):
	
	txt = ""
	for lab_code, d1, d2, typ in dates:
		if typ == "R":
			txt += '''
					R_Date("%s", %f, %f);
			''' % (lab_code, d1, d2)
		elif typ == "U":
			txt += '''
				Date("%s", U(CE(%f), CE(%f)));
			''' % (lab_code, d1, d2)
		else:
			raise Exception("Invalid date type specified: %s (must be 'R' or 'U')" % (typ))
	return txt

def gen_oxcal_calib(dates, curve_name = 'intcal20.14c'):
	# Generate OxCal calibration model
	# dates = [(name, c14age, uncert, 'R'), ...]
	
	return '''
Curve("%s","%s");
Plot()
{
	%s
};
	''' % (curve_name, curve_name, gen_dates(dates))

def gen_sequence(name, data):
	
	txt = ""
	for phase in sorted(list(data.keys())):
		txt += '''
		Boundary("Start %(name)s-%(phase)d");
		Phase("%(name)s-%(phase)d")
		{
			%(dates)s
		};
		Boundary("End %(name)s-%(phase)d");
		''' % dict(name = name, phase = phase, dates = gen_dates(data[phase]))
	
	return '''
	Sequence(%s)
	{
		%s
	};
	''' % (name, txt)
	
def gen_contiguous(name, data):
	
	txt = ""
	last_phase = None
	for phase in sorted(list(data.keys())):
		if last_phase is None:
			txt += '''
		Boundary("Start %s-%d");
			''' % (name, phase)
		else:
			txt += '''
		Boundary("Transition %s-%d/%s-%d");
			''' % (name, last_phase, name, phase)
		txt += '''
		Phase("%s-%d")
		{
			%s
		};
		''' % (name, phase, gen_dates(data[phase]))
		last_phase = phase
	txt += '''
		Boundary("End %s-%d");
	''' % (name, last_phase)
	
	return '''
	Sequence(%s)
	{
		%s
	};
	''' % (name, txt)

def gen_overlapping(name, data):
	
	txt = ""
	for phase in sorted(list(data.keys())):
		txt += '''
		Sequence()
		{
			Boundary("Start %(name)s-%(phase)d");
			Phase("%(name)s-%(phase)d")
			{
				%(dates)s
			};
			Boundary("End %(name)s-%(phase)d");
		};
		''' % dict(name = name, phase = phase, dates = gen_dates(data[phase]))
	return '''
	Phase(%s)
	{
		%s
	};
	''' % (name, txt)

def gen_none(name, data):
	
	txt = ""
	for phase in sorted(list(data.keys())):
		txt += '''
		Label("%s-%d");
		%s
		''' % (name, phase, gen_dates(data[phase]))
	return txt

def gen_oxcal_model(dates, phases, model, curve_name = 'intcal20.14c'):
	# dates = {sample: (age, uncertainty, 'R'), sample: (CE from, CE to, 'U'), ...}
	# phases = {group: {sample: phase, ...}, ...}
	
	model_fncs = {
		'sequence': gen_sequence,
		'contiguous': gen_contiguous,
		'overlapping': gen_overlapping,
		'none': gen_none,
	}
	if model not in model_fncs:
		raise Exception("Invalid model specified: %s" % (model))
	
	txt = ''
	for group in phases:
		data = defaultdict(list)
		for sample in phases[group]:
			data[phases[group][sample]].append([sample] + list(dates[sample]))
		if len(data) > 1:
			txt += model_fncs[model]("Gr.%d" % (group), data)
		else:
			txt += model_fncs['none']("Gr.%d" % (group), data)
	
	return '''
Curve("%s","%s");
Plot()
{
	%s
};
	''' % (curve_name, curve_name, txt)


def r_dates_to_dates(long_lived, r_dates, curve_name, result_path):
	# Convert radiocarbon dates to dates for modeling
	
	# Calibrate radiocarbon dates of long-lived samples and calculate their 2-sigma ranges
	ranges_2sig = {}
	if any(long_lived.values()):
		txt = gen_oxcal_calib([[sample] + list(r_dates[sample]) + ['R'] for sample in long_lived], curve_name)
		fcalib = os.path.join(result_path, "calib.oxcal")
		with open(fcalib, "w") as file:
			file.write(txt)
		run_oxcal(fcalib)
		data = OxCalData(os.path.join(result_path, "calib.js"))
		likelihoods = data.get_likelihoods()
		for name in likelihoods:
			ranges_2sig[name] = likelihoods[name]['range']
	# ranges_2sig = {sample: (lower, upper), ...}
	
	# Generate dates for modeling
	dates = {}
	for sample in r_dates:
		if long_lived[sample]:
			# If the sample is long-lived, use the oldest date in the 2-sigma range as lower boundary and 1950 CE as upper boundary
			dates[sample] = (ranges_2sig[sample][0], 1950, 'U')
		else:
			dates[sample] = (r_dates[sample][0], r_dates[sample][1], 'R')
	# dates = {sample: (age, uncertainty, 'R'), sample: (CE from, CE to, 'U'), ...}
	
	return dates

def run_model(dates, phases, model, curve_name, result_path):
	txt = gen_oxcal_model(dates, phases, model, curve_name)
	fmodel = os.path.join(result_path, "model.oxcal")
	with open(fmodel, "w") as file:
		file.write(txt)
	run_oxcal(fmodel)
	data = OxCalData(os.path.join(result_path, "model.js"))
	return data

def update_model_by_clustering(oc_data, clu_data, result_path):
	
	samples = oc_data.get_samples()
	context_samples = oc_data.get_context_samples()
	context_area = oc_data.get_context_area()
	areas = oc_data.get_areas()
	earlier_than = oc_data.get_earlier_than()
	long_lived = oc_data.get_long_lived()
	r_dates = oc_data.get_r_dates()
	dates = oc_data.get_dates()
	model = oc_data.get_model()
	curve_name = oc_data.get_curve()
	
	clusters = clu_data.get_clusters()
	means = clu_data.get_means()
	opt_n = clu_data.get_opt_n()
	
	if opt_n is None:
		raise ValueError("Optimal number of clusters not found")
	
	clusters = clusters[opt_n]  # clusters = {label: [sample, ...], ...}
	means = means[opt_n]  # means = {label: mean, ...}
	
	# Update earlier_than and phasing based on temporal clustering of the samples
	earlier_than = update_earlier_than_matrix(earlier_than, samples, clusters, means)
	groups, phases = get_groups_and_phases(earlier_than, samples)
	
	# Generate OxCal model
	txt = gen_oxcal_model(dates, phases, model, curve_name)
	fmodel = os.path.join(result_path, "clu_model.oxcal")
	with open(fmodel, "w") as file:
		file.write(txt)

	# Run OxCal model
	run_oxcal(fmodel)
	
	# Save results
	data = OxCalData(os.path.join(result_path, "clu_model.js"))
	data.set_priors(samples, context_samples, context_area, areas, groups, phases, earlier_than, long_lived, r_dates, dates, model, curve_name)
	
	return data
