import numpy as np

# Source: https://github.com/milzj/SAA4PDE/blob/semilinear_complexity/base/signif.py
def signif(x, precision=3):
	"""Rounds the input to significant figures.

	Parameters:
	----------
		x : float
			a floating point number

		precision : int (optional)
			number of significant figures

	"""
	y = np.format_float_positional(x, precision=precision, unique=True, trim="k", fractional=False)
	return np.float64(y)

# Source: https://github.com/milzj/SAA4PDE/blob/semilinear_complexity/base/lsqs_label.py
def lsqs_label(constant=0.0, rate=0.0, base=10.0, precision=3):
	constant = signif(constant, precision=precision)
	rate = signif(rate, precision=precision)
	return r"${}\cdot {}^{}$".format(constant, base, "{"+ str(rate)+"}")

def lsqs_label_base(rate=2, base=2):
	return r"${}={}^{}$".format("n_{\mathrm{ref}}", base, "{"+ str(rate)+"}")

