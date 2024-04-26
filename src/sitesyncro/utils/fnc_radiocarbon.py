import os
import numpy as np
from scipy.interpolate import interp1d

def get_curve(curve_name = 'intcal20.14c'):
	# Load calibration curve
	# returns: [[calendar year BP, C-14 year, uncertainty], ...], sorted by calendar years
	
	fcurve = os.path.join("OxCal\\bin", curve_name)
	
	if not os.path.isfile(fcurve):
		raise ValueError("Calibration curve not found")
	
	with open(fcurve, "r", encoding="latin1") as f:
		data = f.read()
	data = data.split("\n")
	cal_curve = []
	for line in data:
		line = line.strip()
		if not line:
			continue
		if line.startswith("#"):
			continue
		cal_curve.append([np.float64(value) for value in line.split(",")])
	cal_curve = np.array(cal_curve, dtype = np.float64)
	cal_curve = cal_curve[np.argsort(cal_curve[:,0])]
	
	years = np.arange(np.floor(cal_curve[:,0].min()), np.ceil(cal_curve[:,0].max()) + 1, 1)
	cal_curve = np.vstack((
		years,
		interp1d(cal_curve[:,0], cal_curve[:,1], kind = "quadratic")(years),
		interp1d(cal_curve[:,0], cal_curve[:,2], kind = "linear")(years),
	)).T
	
	return cal_curve.astype(np.float64)

def calibrate(age, uncert, curve, date_type = 'R'):
	# Calibrate a C-14 date
	# Calibration formula as defined by Bronk Ramsey 2008, doi: 10.1111/j.1475-4754.2008.00394.x
	#
	# age: C-14 age (years BP) for date_type 'R'; mean calendar age (years BP) for date_type 'U'
	# uncert: uncertainty (years BP) for date_type 'R'; 1/2 range (years BP) for date_type 'U'
	# curve: [[calendar year BP, C-14 year, uncertainty], ...]
	# date_type: 'R' for radiocarbon date; 'U' for calendar date as a uniform distribution
	#
	# returns np.array([p, ...]), where p is the probability of the calendar year
	
	if date_type == 'R':
		sigma_sum = uncert**2 + curve[:,2]**2
		dist = (np.exp(-(age - curve[:,1])**2 / (2 * sigma_sum)) / np.sqrt(sigma_sum))
	elif date_type == 'U':
		dist = np.ones(curve.shape[0], dtype = np.float64)
		dist[curve[:,0] < age - uncertainty] = 0
		dist[curve[:,0] > age + uncertainty] = 0
	else:
		raise Exception("Invalid date type specified: %s (must be 'R' or 'U')" % (date_type))
	
	s = dist.sum()
	if s > 0:
		dist /= s
	
	return dist
