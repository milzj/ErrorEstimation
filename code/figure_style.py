# https://raw.githubusercontent.com/milzj/SAA4PDE/semilinear_complexity/stats/figure_style.py

import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from cycler import cycler

from mpl_toolkits.axes_grid1.inset_locator import inset_axes, mark_inset

fontsize = 12.5
linewidth = 1

color_styles = 2*['k']
marker_styles = ['d', 'o', 'p', 's', '<']
line_styles = ['-', '--', "-.", ':', '--']
line_styles = ['-', "-."]

plt.rcParams.update({
	"text.usetex": True,
	"text.latex.preamble": r"\usepackage{amsfonts}",
	"font.family": "serif",
	"font.serif": "Computer Modern Roman",
	"font.monospace": "Computer Modern Typewriter",
	"lines.linewidth": linewidth,
	"font.size": fontsize,
	"legend.frameon": True,
	"legend.fontsize": "medium",
	"legend.framealpha": 1.0,
	"figure.figsize": [5.0, 5.0],
	"savefig.bbox": "tight",
	"savefig.pad_inches": 0.1})

from matplotlib.colors import LinearSegmentedColormap

def cmap_blue_orange():
	"Color map inspired by cm.coolwarm."
	return LinearSegmentedColormap.from_list(name="cmap_BlueOrange",
           colors =["tab:blue", "lightgrey", "tab:orange"], N=256)


def latex_float(f, signif = 1):
	float_str = "{0:.2g}".format(f)
	if signif == 1:
		float_str = "{0:1.0e}".format(f)
	elif signif == 2:
		float_str = "{0:1.1e}".format(f)
	if "e" in float_str:
		base, exponent = float_str.split("e")
		return r"{0} \cdot 10^{{{1}}}".format(base, int(exponent))
	else:
		return float_str

