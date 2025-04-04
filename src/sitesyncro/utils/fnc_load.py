import os
from collections import defaultdict
from typing import List, Dict, Any

from sitesyncro.utils.fnc_phase import (eap_to_int)


def load_input(fname: str) -> List[Dict[str, Any]]:
	"""
	Load input data from a file.

	This function reads a file line by line, parses each line into a dictionary, and appends the dictionary to a list.
	The file should be in the following format:
	- Each line represents a data record.
	- Data fields are separated by semicolons.
	- The first line is a header and is skipped.
	- Each line should have 10 fields: Sample, Context, Area, C14 Age, Uncertainty, Excavation Area Phase, Earlier-Than, Site Phase, Long-Lived, Redeposited, Outlier.
	
	Parameters:
	fname (str): The name of the file to load.

	Returns:
	A list of dictionaries representing the data records in the file.

	Raises:
	Exception: If a line does not have exactly 10 fields, or if a field cannot be converted to the required type, an exception is raised.
	"""
	data = []
	with open(fname, "r", encoding="utf-8") as file:
		next(file)  # Skip header
		for line in file:
			line = line.strip()
			if line:
				elements = line.split(";")
				if len(elements) != 11:
					raise Exception(f"Incorrect data format in line: {line}")
				elements = [val.strip() for val in elements]
				sample, context, area, age, uncertainty, eap, earlier_than, site_phase, long_lived, redeposited, outlier = elements
				
				if age and uncertainty:
					try:
						age = float(age)
						uncertainty = float(uncertainty)
					except:
						raise Exception(f"Incorrect data format in line: {line}")
				else:
					age, uncertainty = None, None
				try:
					long_lived = int(long_lived) if long_lived else 0
					redeposited = int(redeposited) if redeposited else 0
					outlier = int(outlier) if outlier else 0
				except ValueError:
					raise Exception(f"Incorrect data format in line: {line}")
				# EAP should be in format "1" or "1a" or "1-2" or "1a-b" or "1a-2b"
				if not eap:
					eap = None
				elif eap_to_int(eap) is None:
					raise Exception(f"Incorrect EAP format in line: {line}")
				
				# Site Phase should be in format "1" or "1a" or "1-2" or "1a-b" or "1a-2b"
				if not site_phase:
					site_phase = None
				elif eap_to_int(site_phase) is None:
					raise Exception(f"Incorrect Site Phase format in line: {line}")
				
				# split comma-separated values from earlier_than
				if earlier_than:
					earlier_than = [val.strip() for val in earlier_than.split(",")]
				else:
					earlier_than = []
				
				if not context:
					context = None
				if not area:
					area = None
				data.append({
					"Sample": sample,
					"Context": context,
					"Area": area,
					"C14 Age": age,
					"Uncertainty": uncertainty,
					"EAP": eap,
					"Earlier-Than": earlier_than,
					"Site Phase": site_phase,
					"Long-Lived": long_lived,
					"Redeposited": redeposited,
					"Outlier": outlier,
				})
	return data


def get_samples_contexts_and_areas(data: List[Dict[str, Any]]) -> (
		List[str], Dict[str, List[str]], Dict[str, Any], List[str]):
	
	# Create a set of samples
	samples = set([line["Sample"] for line in data])
	
	# Check if sample names are unique
	done = set()
	found = False
	for line in data:
		if line["Sample"] in done:
			print("Sample %s is duplicate" % (line["Sample"]))
			found = True
		done.add(line["Sample"])
	if found:
		raise Exception("Duplicate sample names found")
	
	# Create a list of samples for each context
	context_samples = defaultdict(list)
	for line in data:
		context_samples[line["Context"]].append(line["Sample"])
	
	# Create a dictionary of contexts and their areas
	context_area = {}
	areas = set()
	for line in data:
		context_area[line["Context"]] = line["Area"]
		areas.add(line["Area"])
	areas = sorted(list(areas))
	
	# Collect sample-less contexts from the Earlier-Than lists which stratigraphically connect contexts with samples
	connected_contexts = defaultdict(set)
	for line in data:
		for context in line["Earlier-Than"]:
			connected_contexts[context].add(line["Context"])
	
	for line in data:
		for context in line["Earlier-Than"]:
			if (context not in context_samples) and (context in connected_contexts) and (len(connected_contexts[context]) >= 2):
				sample = "%s_no-sample" % (context)
				samples.add(sample)
				context_samples[context] = [sample]
				context_area[context] = line["Area"]
	
	samples = sorted(list(samples))
	
	return samples, context_samples, context_area, areas


def get_long_lived(data: List[Dict[str, Any]]) -> Dict[str, bool]:
	# Create a list of long-lived contexts
	long_lived = {}
	for line in data:
		long_lived[line["Sample"]] = bool(line["Long-Lived"])
	return long_lived


def get_redeposited(data: List[Dict[str, Any]]) -> Dict[str, bool]:
	# Create a list of redeposited contexts
	redeposited = {}
	for line in data:
		redeposited[line["Sample"]] = bool(line["Redeposited"])
	return redeposited


def get_outlier(data: List[Dict[str, Any]]) -> Dict[str, bool]:
	# Create a list of outliers
	outlier = {}
	for line in data:
		outlier[line["Sample"]] = bool(line["Outlier"])
	return outlier


def get_c14_dates(data: List[Dict[str, Any]]) -> Dict[str, Any]:
	r_dates = {}
	for line in data:
		age, uncertainty = line["C14 Age"], line["Uncertainty"]
		if None not in [age, uncertainty]:
			r_dates[line["Sample"]] = (age, uncertainty)
	return r_dates


def get_context_eap(data: List[Dict[str, Any]]) -> Dict[str, str]:
	# Create a dictionary of contexts and their excavation area phases
	context_eap = {}
	for line in data:
		context = line["Context"]
		if (context not in context_eap) or (context_eap[context] is None):
			context_eap[context] = line["EAP"]
	return context_eap

def get_context_site_phases(data: List[Dict[str, Any]]) -> Dict[str, str]:
	# Create a dictionary of contexts and their site phases
	context_site_phases = {}
	for line in data:
		context = line["Context"]
		if (context not in context_site_phases) or (context_site_phases[context] is None):
			context_site_phases[context] = line["Site Phase"]
	return context_site_phases

def get_earlier_than_relations(data: List[Dict[str, Any]], context_samples: Dict[str, List[str]]) -> Dict[
	str, List[str]]:
	earlier_than_rel = defaultdict(list)
	for line in data:
		sample1 = line["Sample"]
		for et in line["Earlier-Than"]:
			if et in context_samples:
				for sample2 in context_samples[et]:
					earlier_than_rel[sample1].append(sample2)
	return dict(earlier_than_rel)


def create_result_path(result_path: str, existing: bool = False) -> str:
	"""
	Create a directory for storing results.

	This function creates a new directory for storing results. If a directory with the same name already exists, a new directory with a suffix is created.

	Parameters:
	result_path (str): The path of the directory to create.
	existing (bool, optional): If True, the function will return the path if the directory already exists. If False, the function will create a new directory with a suffix. Default is False.

	Returns:
	str: The path of the created directory.
	"""
	if existing and os.path.isdir(result_path):
		return result_path
	
	parent_dir = os.path.dirname(result_path)
	if not parent_dir:
		parent_dir = "."
	n = None
	for d in os.listdir(parent_dir):
		if d == os.path.basename(result_path):
			n = 0
		elif d.startswith(os.path.basename(result_path)):
			try:
				n = max(n, int(d.split("_")[-1]))
			except:
				pass
	if n is not None:
		n += 1
		result_path = os.path.join(parent_dir, os.path.basename(result_path) + "_" + str(n))
	os.makedirs(result_path)
	return result_path


def load_data(fname: str) -> (
		List[str], Dict[str, str], Dict[str, Any], Dict[str, bool], Dict[str, bool], Dict[str, bool], Dict[str, Any],
		Dict[str, float], Dict[str, float], Dict[str, List[str]]):
	"""
	Load data from a file and process it.

	This function reads a file, processes the data, and returns it in a structured format.

	Parameters:
	fname (str): The name of the file to load.

	Returns:
	(samples, contexts, context_area, long_lived, redeposited, outlier, r_dates, context_eap, context_site_phases, earlier_than_rel):
		- samples = [sample, ...]
		- contexts = {sample: context, ...}
		- context_area = {context: area, ...}
		- long_lived = {sample: True/False, ...}
		- redeposited = {sample: True/False, ...}
		- outlier = {sample: True/False, ...}
		- r_dates = {sample: (age, uncertainty), ...}
		- context_eap = {context: eap, ...}
		- context_site_phases = {context: site_phase, ...}
		- earlier_than_rel = {sample: [sample, ...], ...}
	"""
	
	data = load_input(fname)
	samples, context_samples, context_area, areas = get_samples_contexts_and_areas(data)
	long_lived = get_long_lived(data)
	redeposited = get_redeposited(data)
	outlier = get_outlier(data)
	r_dates = get_c14_dates(data)
	context_eap = get_context_eap(data)
	context_site_phases = get_context_site_phases(data)
	earlier_than_rel = get_earlier_than_relations(data, context_samples)
	
	contexts = {}
	for context in context_samples:
		for sample in context_samples[context]:
			contexts[sample] = context
		if context not in context_area:
			context_area[context] = None
		if context not in context_eap:
			context_eap[context] = None
		if context not in context_site_phases:
			context_site_phases[context] = None
	
	for sample in samples:
		if sample not in contexts:
			contexts[sample] = None
		if sample not in long_lived:
			long_lived[sample] = False
		if sample not in redeposited:
			redeposited[sample] = False
		if sample not in outlier:
			outlier[sample] = False
		if sample not in earlier_than_rel:
			earlier_than_rel[sample] = []
	
	return samples, contexts, context_area, long_lived, redeposited, outlier, r_dates, context_eap, context_site_phases, earlier_than_rel
