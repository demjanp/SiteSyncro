from collections import defaultdict

from sitesyncro.utils.fnc_oxcal import (gen_sequence, gen_contiguous, gen_overlapping, gen_none, gen_multiphase)

class MOxCal(object):
	"""
	A class implementing OxCal modeling functionality
	
	:param model: The parent Model object
	:type model: Model
	"""
	def __init__(self, model: 'Model'):
		
		self.model = model
	
	def gen_oxcal_model(self) -> str:
		"""
		Generate an OxCal model string.

		This function generates an OxCal model string based on the provided model object.
		
		Returns:
		str: The generated OxCal model string.

		Raises:
		Exception: If an invalid phase model is specified.
		"""
		model_fncs = {
			'sequence': gen_sequence,
			'contiguous': gen_contiguous,
			'overlapping': gen_overlapping,
			'none': gen_none,
		}
		if self.model.phase_model not in model_fncs:
			raise Exception("Invalid model specified: %s" % (self.model.phase_model))
		
		txt = ''
		groups = self.model.groups
		for group in groups:
			data_phase = {}
			data_multiphase = defaultdict(list)
			for name in groups[group]:
				phase_min, phase_max = self.model.samples[name].phasing_range
				if phase_min not in data_phase:
					data_phase[phase_min] = []
				if phase_max not in data_phase:
					data_phase[phase_max] = []
				if phase_min == phase_max:
					data_phase[phase_min].append(self.model.samples[name])
				else:
					data_multiphase[(phase_min, phase_max)].append(self.model.samples[name])
			data_phase = dict(data_phase)
			for phase in data_phase:
				data_phase[phase] = sorted(data_phase[phase], key=lambda sample: sum(sample.get_range()))
				if data_phase[phase]:
					data_phase[phase] = "\n".join([sample.to_oxcal() for sample in data_phase[phase]])
				else:
					data_phase[phase] = ""
			txt += model_fncs[self.model.phase_model]("Gr.%d" % (group), data_phase)
			if data_multiphase:
				for key in data_multiphase:
					phase_min, phase_max = key
					if data_multiphase[key]:
						data_multiphase[key] = "\n".join([sample.to_oxcal() for sample in data_multiphase[key]])
					else:
						data_multiphase[key] = ""
				txt += gen_multiphase("Gr.%d" % (group), data_multiphase)
		
		txt = '''
Curve("%s","%s");
Plot()
{
	Outlier_Model("Charcoal",Exp(1,-10,0),U(0,3),"t");
	Outlier_Model("General",T(5),U(0,4),"t");
	%s
};
		''' % (self.model.curve_name, self.model.curve_name, txt)
		
		# Replace non-ASCII characters with their closest ASCII equivalent
		# txt = ''.join(c for c in unicodedata.normalize('NFKD', txt) if unicodedata.category(c) != 'Mn')
		return txt

