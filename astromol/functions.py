import os, sys, argparse, math
import numpy as np
from operator import itemgetter
from math import ceil
from datetime import date
from matplotlib.ticker import AutoMinorLocator
import matplotlib.pyplot as plt
import matplotlib
import matplotlib.ticker as ticker
from matplotlib.ticker import FormatStrFormatter
import matplotlib.gridspec as gridspec
import periodictable as pt
import matplotlib.patches as patches
from colour import Color
import seaborn as sns
from scipy.stats import gaussian_kde as gkde
from scipy.interpolate import make_interp_spline, BSpline

from astromol.molecules import *
from astromol.sources import *
from astromol.telescopes import *

matplotlib.rc("text", usetex=True)
matplotlib.rc(
    "text.latex", preamble=r"\usepackage{cmbright}\usepackage[version=4]{mhchem}"
)

#############################################################
# 						Functions	 						#
#############################################################

def update_plots(mol_list=None):

    """
    A meta function that, when run, will call every plot command and generate new plots based on
    the input list of Molecule objects using default parameters.  Useful for rapidly re-generating all figures.

    Defaults to using all_molecules.
    """

    if mol_list is None:
        mol_list = all_molecules

    # cumu_det_plot(mol_list)
    # cumu_det_natoms_plot(mol_list)
    # det_per_year_per_atoms(mol_list)
    # facility_shares(scopes_list, mol_list)
    # cumu_det_facility(mol_list)
    # periodic_heatmap(mol_list)
    # mass_by_wavelengths(mol_list)
    # mols_waves_by_atoms(mol_list)
    # du_histogram(mol_list)
    # type_pie_chart(mol_list)
    # source_pie_chart(mol_list)
    # mol_type_by_source_type(mol_list)
    # du_by_source_type(mol_list)
    # rel_du_by_source_type(mol_list)
    # mass_by_source_type(mol_list)
    # waves_by_source_type(mol_list)

    return
