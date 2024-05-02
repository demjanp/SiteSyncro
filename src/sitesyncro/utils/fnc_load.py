import os
from collections import defaultdict
from typing import List, Dict, Any


def load_input(fname: str) -> List[Dict[str, Any]]:
	"""
	Load input data from a file.

	This function reads a file line by line, parses each line into a dictionary, and appends the dictionary to a list.
	The file should be in the following format:
	- Each line represents a data record.
	- Data fields are separated by semicolons.
	- The first line is a header and is skipped.
	- Each line should have 10 fields: Sample, Context, Area, C14 Age, Uncertainty, Phase, Earlier-Than, Long-Lived, Redeposited, Outlier.
	- The fields are processed as follows:
		- Sample, Context, and Area are stripped of leading and trailing whitespace.
		- C14 Age and Uncertainty are converted to floats.
		- Phase is converted to a float if possible, otherwise it is set to None.
		- Earlier-Than is split on commas and stripped of leading and trailing whitespace. If it is empty, it is set to an empty list.
		- Long-Lived, Redeposited, and Outlier are converted to integers and then to booleans.

	Parameters:
	fname (str): The name of the file to load.

	Returns:
	A list of dictionaries representing the data records in the file.

	Raises:
	Exception: If a line does not have exactly 10 fields, or if a field cannot be converted to the required type, an exception is raised.
	"""
	data = []
	with open(fname, "r") as file:
		next(file)  # Skip header
		for line in file:
			line = line.strip()
			if line:
				elements = line.split(";")
				if len(elements) != 10:
					raise Exception(f"Incorrect data format in line: {line}")
				sample, context, area, age, uncertainty, phase, earlier_than, long_lived, redeposited, outlier = elements
				try:
					age = float(age.strip())
					uncertainty = float(uncertainty.strip())
					long_lived = int(long_lived.strip())
					redeposited = int(redeposited.strip())
					outlier = int(outlier.strip())
				except ValueError:
					raise Exception(f"Incorrect data format in line: {line}")
				phase = phase.strip()
				try:
					phase = float(phase)
				except:
					phase = None
				# split comma-separated values from earlier_than
				earlier_than = earlier_than.strip()
				if earlier_than:
					earlier_than = [val.strip() for val in earlier_than.split(",")]
				else:
					earlier_than = []
				context = context.strip()
				if not context:
					context = None
				area = area.strip()
				if not area:
					area = None
				data.append({
					"Sample": sample.strip(),
					"Context": context,
					"Area": area,
					"C14 Age": age,
					"Uncertainty": uncertainty,
					"Phase": phase,
					"Earlier-Than": earlier_than,
					"Long-Lived": long_lived,
					"Redeposited": redeposited,
					"Outlier": outlier,
				})
	return data


def get_samples_contexts_and_areas(data: List[Dict[str, Any]]) -> (
		List[str], Dict[str, List[str]], Dict[str, Any], List[str]):
	# Create a list of samples
	samples = sorted(list(set([line["Sample"] for line in data])))
	
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
		r_dates[line["Sample"]] = (line["C14 Age"], line["Uncertainty"])
	return r_dates


def get_context_phase(data: List[Dict[str, Any]]) -> Dict[str, float]:
	# Create a dictionary of contexts and their phases
	context_phase = {}
	for line in data:
		context_phase[line["Context"]] = line["Phase"]
	return context_phase


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
		Dict[str, float], Dict[str, List[str]]):
	"""
	Load data from a file and process it.

	This function reads a file, processes the data, and returns it in a structured format.

	Parameters:
	fname (str): The name of the file to load.

	Returns:
	(samples, contexts, context_area, long_lived, redeposited, outlier, r_dates, context_phase, earlier_than_rel):
		- samples = [sample, ...]
		- contexts = {sample: context, ...}
		- context_area = {context: area, ...}
		- long_lived = {sample: True/False, ...}
		- redeposited = {sample: True/False, ...}
		- outlier = {sample: True/False, ...}
		- r_dates = {sample: (age, uncertainty), ...}
		- context_phase = {context: phase, ...}
		- earlier_than_rel = {sample: [sample, ...], ...}
	"""
	
	data = load_input(fname)
	samples, context_samples, context_area, areas = get_samples_contexts_and_areas(data)
	long_lived = get_long_lived(data)
	redeposited = get_redeposited(data)
	outlier = get_outlier(data)
	r_dates = get_c14_dates(data)
	context_phase = get_context_phase(data)
	earlier_than_rel = get_earlier_than_relations(data, context_samples)
	
	contexts = {}
	for context in context_samples:
		for sample in context_samples[context]:
			contexts[sample] = context
		if context not in context_area:
			context_area[context] = None
		if context not in context_phase:
			context_phase[context] = None
	
	for sample in samples:
		if sample not in contexts:
			contexts[sample] = None
		if sample not in long_lived:
			long_lived[sample] = False
		if sample not in redeposited:
			redeposited[sample] = False
		if sample not in outlier:
			outlier[sample] = False
		if sample not in r_dates:
			r_dates[sample] = (None, None)
		if sample not in earlier_than_rel:
			earlier_than_rel[sample] = []
	
	return samples, contexts, context_area, long_lived, redeposited, outlier, r_dates, context_phase, earlier_than_rel
