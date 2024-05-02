from typing import Union

import numpy as np


def dict_keys_to_int(data: Union[dict, list]) -> Union[dict, list]:
	"""
	Convert all numeric keys in a dictionary or list to integers.

	This function is useful when loading JSON files, as JSON objects can only have string keys.
	If a key is a string representation of an integer, it will be converted to an integer.
	The function works recursively, so it will convert keys in any nested dictionaries as well.

	Parameters:
	data: The input data. It can be a dictionary or a list.
	If it's a dictionary, the function will convert its keys.
	If it's a list, the function will iterate over its elements and convert keys in any dictionaries found.

	Returns:
	The input data with all numeric keys converted to integers.
	The type of the returned value will be the same as the type of the input.
	"""
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


def dict_np_to_list(data: Union[dict, list]) -> Union[dict, list]:
	"""
	Convert all numpy arrays in a dictionary or list to lists.

	This function is useful when preparing data for serialization, as numpy arrays cannot be serialized directly.
	If a value is a numpy array, it will be converted to a list.
	The function works recursively, so it will convert values in any nested dictionaries or lists as well.

	Parameters:
	data: The input data. It can be a dictionary or a list.
	If it's a dictionary, the function will convert its values.
	If it's a list, the function will iterate over its elements and convert values in any dictionaries or arrays found.

	Returns:
	The input data with all numpy arrays converted to lists.
	The type of the returned value will be the same as the type of the input.
	"""
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
