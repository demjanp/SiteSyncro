import re
import copy
import json
import gzip
import numpy as np

class AbstractData(object):
	
	def __init__(self):
		self._data = {}
	
	def _convert_keys(self, data):
		if isinstance(data, dict):
			new_data = data.copy()  # Create a copy of the dictionary
			for key in new_data:
				new_data[key] = self._convert_keys(new_data[key])
				if isinstance(key, str) and key.isdigit():
					data[int(key)] = data.pop(key)
		elif isinstance(data, list):
			for i, val in enumerate(data):
				data[i] = self._convert_keys(val)
		return data
	
	def save_json(self, fname):
		
		if fname.endswith('.gz'):
			with gzip.open(fname, 'wt') as file:
				json.dump(self._data, file)
		else:
			with open(fname, 'w') as file:
				json.dump(self._data, file)
	
	def load_json(self, fname):
		
		if fname.endswith('.gz'):
			with gzip.open(fname, 'rt') as file:
				data = json.load(file)
		else:
			with open(fname, 'r') as file:
				data = json.load(file)
		
		# Recursively convert all numeric keys in self._data to integers
		self._data = self._convert_keys(data)
	
	def has_data(self):
		
		return len(self._data) > 0
	
	def __repr__(self):
		
		def format_list(value):
			
			if not isinstance(value, list):
				return value
			if len(value) <= 3:
				return value
			return "[%s, ... (%d items)]" % (", ".join([str(val) for val in value[:3]]), len(value))
		
		def print_data(data, chain=""):
			text = ""
			for key in sorted(data.keys()):
				if isinstance(key, int):
					new_chain = chain + "[%s]" % (key)
				else:
					new_chain = chain + "['%s']" % (key)
				if isinstance(data[key], dict):
					text += print_data(data[key], new_chain)
				else:
					text += "%s = %s\n" % (new_chain, format_list(data[key]))
			return text
		
		return print_data(self._data)
	
class OxCalData(AbstractData):
	
	def __init__(self, fname=None):
		
		def replace_colons_in_braces(s):
			# Find the innermost braces
			inner_braces = re.findall(r'\{[^{}]*\}', s)
			if not inner_braces:
				return s
			
			for brace in inner_braces:
				# Replace colons in the innermost braces
				s = s.replace(brace, "dict(%s)" % brace.replace(':', '=')[1:-1])
			
			# Recursively process the modified string
			return replace_colons_in_braces(s)
		
		def read_js_file(file_path):
			
			data = {}
			
			with (open(file_path, 'r') as file):
				for line in file:
					if line.startswith('if('):
						continue
					
					# Replace last character of line with a comma
					line = line[:-2] + ',' + line[-1]
					
					# Put everything before "=" in quotes
					line = line.split('=')
					line = [line[0], '='.join(line[1:])]
					
					key, value = line
					key = key.strip()
					value = value.strip().strip(',')
					
					value = replace_colons_in_braces(value)
					
					# Replace NaN (only if whole word) with None
					value = re.sub(r'\bNaN\b', 'None', value)
					
					# Parse value using python eval
					value = eval(value)
					if value == []:
						value = {}
					
					keys = key.split('.')
					collect = []
					for key in keys:
						index1 = None
						index2 = None
						if '[' in key:
							key, index = key.split('[', 1)
							index1 = int(index.split(']')[0])
							index2 = None
							if '[' in index:
								index2 = int(index.split('[')[1].split(']')[0])
						if index2 is not None:
							collect += [key, index1, index2]
						elif index1 is not None:
							collect += [key, index1]
						else:
							collect.append(key)
					keys = collect
					
					current = data
					for i, key in enumerate(keys):
						if i == len(keys) - 1:
							current[key] = value
						elif key not in current:
							current[key] = {}
						current = current[key]
			
			return data
		
		super().__init__()
		
		self._data['data'] = {}
		if fname is not None:
			self._data['data'] = read_js_file(fname)
		self._data['priors'] = {}
	
	def __getitem__(self, key):
		return self._data['data'][key]
	
	def __setitem__(self, key, value):
		self._data['data'][key] = value
	
	def has_data(self):
		
		return len(self._data['data']) > 0
	
	def set_priors(self, samples, context_samples, context_area, areas, groups, phases, earlier_than, long_lived,
	               r_dates, dates, model, curve):
		# samples = [sample, ...]
		# context_samples = {context: [sample, ...], ...}
		# context_area = {context: area, ...}
		# areas = [area, ...]
		# groups = {group: [sample, ...], ...}
		# phases = {group: {sample: phase, ...}, ...}
		# long_lived = {sample: True/False, ...}
		# r_dates = {sample: (age, uncertainty), ...}
		# dates = {sample: (age, uncertainty, 'R'), sample: (CE from, CE to, 'U'), ...}; R = radiocarbon, U = uniform
		
		self._data['priors']['samples'] = copy.deepcopy(samples)
		self._data['priors']['context_samples'] = copy.deepcopy(context_samples)
		self._data['priors']['context_area'] = copy.deepcopy(context_area)
		self._data['priors']['areas'] = copy.deepcopy(areas)
		self._data['priors']['groups'] = copy.deepcopy(groups)
		self._data['priors']['phases'] = copy.deepcopy(phases)
		self._data['priors']['earlier_than'] = copy.deepcopy(earlier_than.tolist())
		self._data['priors']['long_lived'] = copy.deepcopy(long_lived)
		self._data['priors']['r_dates'] = copy.deepcopy(r_dates)
		self._data['priors']['dates'] = copy.deepcopy(dates)
		self._data['priors']['model'] = model
		self._data['priors']['curve'] = curve
	
	@staticmethod
	def __get_params(data):
		
		mean = None
		if 'mean' in data:
			mean = data['mean']
		prob = data['prob']
		start = data['start']
		resolution = data['resolution']
		years = [start + j * resolution for j in range(len(prob))]
		
		rng = (None, None)
		if ('range' in data) and (2 in data['range']):
			L = np.inf
			U = -np.inf
			for key in data['range'][2]:
				L = min(L, data['range'][2][key][0])
				U = max(U, data['range'][2][key][1])
			rng = (L, U)
		
		return mean, rng, prob, years
	
	def get_posteriors(self):
		result = {}
		
		if 'ocd' not in self._data['data']:
			return result
		
		for i in self._data['data']['ocd']:
			if i == 0:
				continue
			if 'posterior' not in self._data['data']['ocd'][i]:
				continue
			if 0 not in self._data['data']['ocd'][i]['posterior']['comment']:
				continue
			name = self._data['data']['ocd'][i]['posterior']['comment'][0][:-10].strip()
			
			agreement = None
			if 'agreement' in self._data['data']['ocd'][i]['posterior']:
				agreement = self._data['data']['ocd'][i]['posterior']['agreement']
			
			mean, rng, prob, years = self.__get_params(self._data['data']['ocd'][i]['posterior'])
			
			result[name] = {'prob': prob, 'years': years, 'agreement': agreement, 'mean': mean, 'range': rng}
		
		return result
	
	def get_likelihoods(self):
		
		result = {}
		
		if 'ocd' not in self._data['data']:
			return result
		
		for i in self._data['data']['ocd']:
			if i == 0:
				continue
			if 'likelihood' not in self._data['data']['ocd'][i]:
				continue
			if 'prob' not in self._data['data']['ocd'][i]['likelihood']:
				continue
			
			name = self._data['data']['ocd'][i]['name'].strip()
			
			mean, rng, prob, years = self.__get_params(self._data['data']['ocd'][i]['likelihood'])
			
			result[name] = {'prob': prob, 'years': years, 'mean': mean, 'range': rng}
		
		return result
	
	def get_samples(self):
		# samples = [sample, ...]
		
		if 'samples' not in self._data['priors']:
			return None
		
		return copy.deepcopy(self._data['priors']['samples'])
	
	def get_context_samples(self):
		# context_samples = {context: [sample, ...], ...}
		
		if 'context_samples' not in self._data['priors']:
			return None
		
		return copy.deepcopy(self._data['priors']['context_samples'])
	
	def get_context_area(self):
		# context_area = {context: area, ...}
		
		if 'context_area' not in self._data['priors']:
			return None
		
		return copy.deepcopy(self._data['priors']['context_area'])
	
	def get_areas(self):
		# areas = [area, ...]
		
		if 'areas' not in self._data['priors']:
			return None
		
		return copy.deepcopy(self._data['priors']['areas'])
	
	def get_groups(self):
		# groups = {group: [sample, ...], ...}
		
		if 'groups' not in self._data['priors']:
			return None
		
		return copy.deepcopy(self._data['priors']['groups'])
	
	def get_phases(self):
		# phases = {group: {sample: phase, ...}, ...}
		
		if 'phases' not in self._data['priors']:
			return None
		
		return copy.deepcopy(self._data['priors']['phases'])
	
	def get_earlier_than(self):
		# earlier_than: matrix[n_samples x n_samples] = [True/False, ...]; sample in row is earlier than sample in column based on stratigraphy
		
		if 'earlier_than' not in self._data['priors']:
			return None
		
		return np.array(self._data['priors']['earlier_than'], dtype=bool)
	
	def get_long_lived(self):
		# long_lived = {sample: True/False, ...}
		
		if 'long_lived' not in self._data['priors']:
			return None
		
		return copy.deepcopy(self._data['priors']['long_lived'])
	
	def get_r_dates(self):
		# r_dates = {sample: (age, uncertainty), ...}
		
		if 'r_dates' not in self._data['priors']:
			return None
		
		return copy.deepcopy(self._data['priors']['r_dates'])
	
	def get_dates(self):
		# dates = {sample: (age, uncertainty, 'R'), sample: (CE from, CE to, 'U'), ...}
		
		if 'dates' not in self._data['priors']:
			return None
		
		return copy.deepcopy(self._data['priors']['dates'])
	
	def get_model(self):
		
		if 'model' not in self._data['priors']:
			return None
		return self._data['priors']['model']
		
	def get_curve(self):
		
		if 'curve' not in self._data['priors']:
			return None
		return self._data['priors']['curve']
	
class ClusterData(AbstractData):
	
	def __init__(self, clusters = {}, means = {}, sils = {}, ps = {}, p_value = None, opt_n = None):
		
		super().__init__()
		
		self._data['clusters'] = copy.deepcopy(clusters)
		self._data['means'] = copy.deepcopy(means)
		self._data['sils'] = copy.deepcopy(sils)
		self._data['ps'] = copy.deepcopy(ps)
		self._data['p_value'] = p_value
		self._data['opt_n'] = opt_n
	
	def has_data(self):
		
		return len(self._data['clusters']) > 0
	
	def get_clusters(self):
		
		return copy.deepcopy(self._data['clusters'])
	
	def get_means(self):
		
		return copy.deepcopy(self._data['means'])
	
	def get_sils(self):
		
		return copy.deepcopy(self._data['sils'])
	
	def get_ps(self):
		
		return copy.deepcopy(self._data['ps'])
	
	def get_p_value(self):
		
		return self._data['p_value']
	
	def get_opt_n(self):
		
		return self._data['opt_n']
	
	def __getitem__(self, key):
		return self._data[key]
	
	def __setitem__(self, key, value):
		self._data[key] = value

class RandomizeData(AbstractData):
	
	def __init__(self, years = None, sum_obs = None, uniform = None, p = None, p_value = None, sums_rnd_lower = None, sums_rnd_upper = None):
		
		super().__init__()
		
		self._data['years'] = years if years is None else years.tolist()
		self._data['sum_obs'] = sum_obs if sum_obs is None else sum_obs.tolist()
		self._data['uniform'] = uniform
		self._data['p'] = p
		self._data['p_value'] = p_value
		self._data['sums_rnd_lower'] = sums_rnd_lower if sums_rnd_lower is None else sums_rnd_lower.tolist()
		self._data['sums_rnd_upper'] = sums_rnd_upper if sums_rnd_upper is None else sums_rnd_upper.tolist()
	
	def has_data(self):
		
		return self._data['years'] is not None
	
	def get_years(self):
		
		return np.array(self._data['years'])
	
	def get_sum_obs(self):
		
		return np.array(self._data['sum_obs'])
	
	def get_uniform(self):
		
		return self._data['uniform']
	
	def get_p(self):
		
		return self._data['p']
	
	def get_p_value(self):
		
		return self._data['p_value']
	
	def get_sums_rnd_lower(self):
		
		return np.array(self._data['sums_rnd_lower'])
	
	def get_sums_rnd_upper(self):
		
		return np.array(self._data['sums_rnd_upper'])
	
	def __getitem__(self, key):
		return self._data[key]
	
	def __setitem__(self, key, value):
		self._data[key] = value