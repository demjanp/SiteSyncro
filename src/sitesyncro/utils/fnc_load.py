from collections import defaultdict
import os

def load_input(fname):
	data = []
	with open(fname, "r") as file:
		next(file)  # Skip header
		for line in file:
			line = line.strip()
			if line:
				elements = line.split(";")
				if len(elements) != 8:
					raise Exception(f"Incorrect data format in line: {line}")
				sample, context, area, age, uncertainty, phase, earlier_than, long_lived = elements
				try:
					age = float(age.strip())
					uncertainty = float(uncertainty.strip())
					phase = float(phase.strip())
					long_lived = int(long_lived.strip())
				except ValueError:
					raise Exception(f"Incorrect data format in line: {line}")
				# split comma-separated values from earlier_than
				earlier_than = [val.strip() for val in earlier_than.split(",")]
				data.append({
					"Sample": sample.strip(),
					"Context": context.strip(),
					"Area": area.strip(),
					"C14 Age": age,
					"Uncertainty": uncertainty,
					"Phase": phase,
					"Earlier-Than": earlier_than,
					"Long-Lived": long_lived
				})
	return data

def get_samples_contexts_and_areas(data):
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

def get_long_lived(data):
	
	# Create a list of long-lived or redeposited contexts
	long_lived = {}
	for line in data:
		long_lived[line["Sample"]] = bool(line["Long-Lived"])
	return long_lived

def get_c14_dates(data):
	
	r_dates = {}
	for line in data:
		r_dates[line["Sample"]] = (line["C14 Age"], line["Uncertainty"])
	return r_dates

def get_context_phase(data):
	
	# Create a dictionary of contexts and their phases
	context_phase = {}
	for line in data:
		context_phase[line["Context"]] = line["Phase"]
	return context_phase

def get_earlier_than_relations(data, context_samples):
	
	earlier_than_rel = defaultdict(list)
	for line in data:
		sample1 = line["Sample"]
		for et in line["Earlier-Than"]:
			if et in context_samples:
				for sample2 in context_samples[et]:
					earlier_than_rel[sample1].append(sample2)
	return dict(earlier_than_rel)

def create_result_path(result_path, existing = False):
	
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

def load_data(fname):
	data = load_input(fname)
	samples, context_samples, context_area, areas = get_samples_contexts_and_areas(data)
	long_lived = get_long_lived(data)
	r_dates = get_c14_dates(data)
	context_phase = get_context_phase(data)
	earlier_than_rel = get_earlier_than_relations(data, context_samples)
	
	# samples = [sample, ...]
	# contexts = {sample: context, ...}
	# context_area = {context: area, ...}
	# long_lived = {sample: True/False, ...}
	# r_dates = {sample: (age, uncertainty), ...}
	# context_phase = {context: phase, ...}
	# earlier_than_rel = {sample: [sample, ...], ...}
	
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
		if sample not in r_dates:
			r_dates[sample] = (None, None)
		if sample not in earlier_than_rel:
			earlier_than_rel[sample] = {}
	
	return samples, contexts, context_area, long_lived, r_dates, context_phase, earlier_than_rel
