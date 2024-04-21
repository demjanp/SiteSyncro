import numpy as np
from scipy.interpolate import interp1d

def load_calibration_curve(fcalib, interpolate = False):
	# load calibration curve
	# data from: fcalib 14c file
	# returns: [[CalBP, ConvBP, CalSigma], ...], sorted by CalBP
	
	with open(fcalib, "r", encoding="latin1") as f:
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
	
	if interpolate:
		cal_bp = np.arange(cal_curve[:,0].min(), cal_curve[:,0].max() + 1, 1)
		cal_curve = np.vstack((
			cal_bp,
			interp1d(cal_curve[:,0], cal_curve[:,1], kind = "quadratic")(cal_bp),
			interp1d(cal_curve[:,0], cal_curve[:,2], kind = "linear")(cal_bp),
		)).T
	
	return cal_curve.astype(np.float64)

def calibrate(age, uncert, curve, normalize = False):
	# calibrate a 14C measurement
	# calibration formula as defined by Bronk Ramsey 2008, doi: 10.1111/j.1475-4754.2008.00394.x
	# age: uncalibrated 14C age BP
	# uncert: 1 sigma uncertainty
	# curve: [[CalBP, ConvBP, CalSigma], ...]
	
	sigma_sum = uncert**2 + curve[:,2]**2
	dist = (np.exp(-(age - curve[:,1])**2 / (2 * sigma_sum)) / np.sqrt(sigma_sum))
	if normalize:
		s = dist.sum()
		if s > 0:
			dist /= s
	return dist
