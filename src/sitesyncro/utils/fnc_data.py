import os
import numpy as np

from typing import Union

def dict_keys_to_int(data: Union[dict, list]):
	# Convert all numeric keys in a dictionary to int (usefull when loading JSON files)
	if isinstance(data, dict):
		new_data = data.copy()  # Create a copy of the dictionary
		for key in new_data:
			new_data[key] = dict_keys_to_int(new_data[key])
			if isinstance(key, str) and key.isdigit():
				data[int(key)] = data.pop(key)
	elif isinstance(data, list):
		for i, val in enumerate(data):
			data[i] = dict_keys_to_int(val)
	return data

def dict_np_to_list(data: Union[dict, list]):
	# Convert all numpy arrays in a dictionary to lists
	if isinstance(data, dict):
		new_data = data.copy()  # Create a copy of the dictionary
		for key in new_data:
			new_data[key] = dict_np_to_list(new_data[key])
			if isinstance(new_data[key], np.ndarray):
				data[key] = new_data[key].tolist()
	elif isinstance(data, list):
		for i, val in enumerate(data):
			data[i] = dict_np_to_list(val)
	return data
