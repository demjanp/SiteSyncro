import os
import re
import zipfile
from collections import defaultdict
from typing import Dict, Any

import numpy as np
import requests
import unicodedata
from scipy.interpolate import interp1d
from tqdm import tqdm


def download_oxcal(url: str = None) -> bool:
	"""
	Download and unzip OxCal from a specified URL.

	This function checks if OxCal is already downloaded. If not, it downloads and unzips OxCal from a specified URL to
	the current directory.

	Parameters:
	url (str, optional): The URL to download OxCal from. If not provided, the default OxCal distribution URL is used.

	Returns:
	bool: True if OxCal is successfully downloaded and unzipped, False otherwise.
	"""
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
		print(f"Error: Unable to connect to {url}. Please check the URL or your internet connection.")
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


def oxcal_date(name: str, age: float, uncertainty: float, date_type: str, long_lived: bool, outlier: bool) -> str:
	"""
	Generate an OxCal date string.

	This function generates an OxCal date string based on the provided parameters.

	Parameters:
	name (str): The name of the date.
	age (float): The age of the date.
	uncertainty (float): The uncertainty of the date.
	date_type (str): The type of the date. It should be 'R' or 'U'.
	long_lived (bool): Whether the date is long-lived.
	outlier (bool): Whether the date is an outlier.

	Returns:
	str: The generated OxCal date string.

	Raises:
	Exception: If an invalid date type is specified.
	"""
	txt = None
	if date_type == 'R':
		txt = '''R_Date("%s", %f, %f)''' % (name, age, uncertainty)
	elif date_type == 'U':
		txt = '''Date("%s", U(CE(%f), CE(%f)))''' % (name, -(age + uncertainty - 1950), -(age - uncertainty - 1950))
	if txt is None:
		raise Exception("Invalid date type specified: %s (must be 'R' or 'U')" % (date_type))
	if outlier:
		txt += '''{Outlier("General", 1);};'''
	elif long_lived:
		txt += '''{Outlier("Charcoal", 1);};'''
	else:
		txt += ";"
	return txt


def gen_sequence(name: str, data: Dict[int, str]) -> str:
	txt = ""
	for phase in sorted(list(data.keys())):
		txt += '''
		Boundary("Start %(name)s-%(phase)d");
		Phase("%(name)s-%(phase)d")
		{
			%(dates)s
		};
		Boundary("End %(name)s-%(phase)d");
		''' % dict(name=name, phase=phase, dates=data[phase])
	
	return '''
	Sequence(%s)
	{
		%s
	};
	''' % (name, txt)


def gen_contiguous(name: str, data: Dict[int, str]) -> str:
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
		''' % (name, phase, data[phase])
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


def gen_overlapping(name: str, data: Dict[int, str]) -> str:
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
		''' % dict(name=name, phase=phase, dates=data[phase])
	return '''
	Phase(%s)
	{
		%s
	};
	''' % (name, txt)


def gen_none(name: str, data: Dict[int, str]) -> str:
	txt = ""
	for phase in sorted(list(data.keys())):
		txt += '''
		Label("%s-%d");
		%s
		''' % (name, phase, data[phase])
	return txt


def load_oxcal_data(fname: str) -> Dict:
	"""
	Load OxCal data from a file.

	This function loads OxCal data from a specified file.

	Parameters:
	fname (str): The name of the file to load the data from.

	Returns:
	Dict: The loaded OxCal data.

	Raises:
	Exception: If the data file is not found.
	"""
	
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
		
		with (open(file_path, 'r', encoding='utf-8') as file):
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
				value = eval(value, dict(true=True, false=False))
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
	
	if not os.path.isfile(fname):
		raise Exception("Data file %s not found" % (fname))
	
	return read_js_file(fname)


def get_distributions(data: Dict, curve: np.ndarray) -> (Dict[str, Any], Dict[str, Any]):
	"""
	Get distributions from OxCal data.

	This function gets likelihoods and posteriors from OxCal data unified to match the range and calendar ages on the calibration curve.

	Parameters:
	data (Dict): The OxCal data (from load_oxcal_data).
	curve: np.array([[calendar year BP, C-14 year, uncertainty], ...]), sorted by calendar years. The radiocarbon calibration curve.

	Returns:
	(likelihoods, posteriors)
		- likelihoods = {name: [distribution, mean, rng], ...}
		- posteriors = {name: [distribution, mean, rng, agreement], ...}
	"""
	
	def _get_params(params_data):
		mean = None
		if 'mean' in params_data:
			mean = 1950 - params_data['mean']
		prob = params_data['prob']
		start = params_data['start']
		resolution = params_data['resolution']
		years = [start + j * resolution for j in range(len(prob))]
		rng = (None, None)
		if ('range' in params_data) and (2 in params_data['range']):
			L = np.inf
			U = -np.inf
			for key in params_data['range'][2]:
				L = min(L, params_data['range'][2][key][0])
				U = max(U, params_data['range'][2][key][1])
			rng = (1950 - L, 1950 - U)
		return mean, rng, prob, years
	
	def _unify(dist, years, curve):
		# years are in CE
		all_years = curve[:, 0]
		years = 1950 - np.array(years)
		# Make sure years are in the correct order
		if np.sign(years[-1] - years[0]) != np.sign(all_years[-1] - all_years[0]):
			years = years[::-1]
			dist = dist[::-1]
		years = np.concatenate(([all_years[0]], years, [all_years[-1]]))
		dist = np.concatenate(([0], dist, [0]))
		s = dist.sum()
		if s > 0:
			dist /= s
		new_prob = interp1d(years, dist, kind="linear")(all_years)
		s = new_prob.sum()
		if s > 0:
			new_prob /= s
		return new_prob
	
	def _load_likelihoods(data):
		result = {}
		if 'ocd' not in data:
			return
		for i in data['ocd']:
			if i == 0:
				continue
			if 'likelihood' not in data['ocd'][i]:
				continue
			if 'prob' not in data['ocd'][i]['likelihood']:
				continue
			name = data['ocd'][i]['name'].strip()
			mean, rng, prob, years = _get_params(data['ocd'][i]['likelihood'])
			result[name] = [_unify(prob, years, curve), mean, rng]
		return result
	
	def _load_posteriors(data):
		result = {}
		if 'ocd' not in data:
			return result
		for i in data['ocd']:
			if i == 0:
				continue
			if 'posterior' not in data['ocd'][i]:
				continue
			if 'prob' not in data['ocd'][i]['posterior']:
				continue
			if 0 not in data['ocd'][i]['posterior']['comment']:
				continue
			name = data['ocd'][i]['posterior']['comment'][0][:-10].strip()
			agreement = None
			if 'agreement' in data['ocd'][i]['posterior']:
				agreement = data['ocd'][i]['posterior']['agreement']
			mean, rng, prob, years = _get_params(data['ocd'][i]['posterior'])
			result[name] = [_unify(prob, years, curve), mean, rng, agreement]
		return result
	
	likelihoods = _load_likelihoods(data)
	posteriors = _load_posteriors(data)
	return likelihoods, posteriors
