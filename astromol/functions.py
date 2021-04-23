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
import matplotlib.patches as patches
import periodictable
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


def make_all_plots():

    """
    A meta function that, when run, will call every plot command and generate new plots based on
    the input list of Molecule objects using default parameters.  Useful for rapidly re-generating all figures.

    Will always use the full database and default filenames.
    """

    cumu_det_plot(mol_list=all_molecules)
    cumu_det_natoms_plot(mol_list=all_molecules)
    det_per_year_per_atom(mol_list=all_molecules)
    facility_shares(telescopes_list=all_telescopes, mol_list=all_molecules)
    # cumu_det_facility(mol_list=all_molecules)
    # periodic_heatmap(mol_list=all_molecules)
    # mass_by_wavelengths(mol_list=all_molecules)
    # mols_waves_by_atoms(mol_list=all_molecules)
    # du_histogram(mol_list=all_molecules)
    # type_pie_chart(mol_list=all_molecules)
    # source_pie_chart(mol_list=all_molecules)
    # mol_type_by_source_type(mol_list=all_molecules)
    # du_by_source_type(mol_list=all_molecules)
    # rel_du_by_source_type(mol_list=all_molecules)
    # mass_by_source_type(mol_list=all_molecules)
    # waves_by_source_type(mol_list=all_molecules)

    return


def change_color(color, amount=1.0):
    """
    Lightens the given color by multiplying (1-luminosity) by the given amount.
    Input can be matplotlib color string, hex string, or RGB tuple.

    Examples:
    >> lighten_color('g', 0.3)
    >> lighten_color('#F034A3', 0.6)
    >> lighten_color((.3,.55,.1), 0.5)
    """
    import matplotlib.colors as mc
    import colorsys

    try:
        c = mc.cnames[color]
    except:
        c = color
    c = colorsys.rgb_to_hls(*mc.to_rgb(c))
    return colorsys.hls_to_rgb(c[0], 1 - amount * (1 - c[1]), c[2])


#############################################################
# 						    Plots	 						#
#############################################################


def cumu_det_plot(mol_list=None, syear=None, eyear=None, filename=None):

    """
    Makes a plot of the cumulative detections by year.
    Defaults to using all molecules in the database, but can be passed a subset of molecules as a list of Molecule objects in mol_list.
    The start year and end year are defined by default to be the earliest year in the list and the current year, by default, but take integers as overrides.
    The filename defaults to 'cumulative_detections.pdf' but can be overriden.
    """

    # If a list wasn't specified, default to all molecules
    if mol_list is None:
        mol_list = all_molecules

    # Close an old figure if it exists and initialize a new figure
    plt.close("Cumulative Detections")
    fig = plt.figure(num="Cumulative Detections", figsize=(10, 8))
    plt.ion()

    # set some font defaults
    fontparams = {"size": 24, "family": "sans-serif", "sans-serif": ["Helvetica"]}
    plt.rc("font", **fontparams)
    plt.rc("mathtext", fontset="stixsans")

    # get the starting and ending years, if they aren't set by the user
    if syear is None:
        # grab the earliest year in the list of molecules
        syear = min([x.year for x in mol_list])

    if eyear is None:
        # grab the current year
        eyear = date.today().year

    # make the x-axis array of years
    years = np.arange(syear, eyear + 1)

    # make an array to hold the detections
    dets = np.copy(years) * 0

    # loop through the years and the list and add everything up.  There's gotta be a better way to do this, but who cares, its fast enough.
    for x in range(len(dets)):
        i = 0
        for mol in mol_list:
            if mol.year < years[x] + 1:
                i += 1
        dets[x] = i

    # get some year indicies for years we care about
    def iyear(x):
        return np.argwhere(years == x)[0][0]

    # do linear fits to the data for the two ranges we care about (1968-Present, 2005-Present)
    # get the slope of the fit to detections since 1968 and since 2005

    trend1968 = (
        np.polynomial.polynomial.Polynomial.fit(
            years[iyear(1968) :], dets[iyear(1968) :], 1
        )
        .convert()
        .coef[1]
    )
    trend2005 = (
        np.polynomial.polynomial.Polynomial.fit(
            years[iyear(2005) :], dets[iyear(2005) :], 1
        )
        .convert()
        .coef[1]
    )

    # load up an axis
    ax = fig.add_subplot(111)

    # label the axes
    plt.xlabel("Year")
    plt.ylabel("Cumulative Number of Detected Molecules")

    # customize tick marks
    ax.tick_params(axis="x", which="both", direction="in", length=15, width=1)
    ax.tick_params(axis="y", which="both", direction="in", length=15, width=1)

    ax.yaxis.set_ticks_position("both")
    ax.xaxis.set_ticks_position("both")

    # plot it
    ax.plot(years, dets, color="dodgerblue", linewidth=4)

    # add annotations; #first, the detections per year
    args = {"ha": "left", "size": "24"}
    det_str = r"\noindent Since 1968: {:.1f} detections/year \\ Since 2005: {:.1f} detections/year".format(
        trend1968, trend2005
    )
    ax.annotate(det_str, xy=(0.05, 0.85), xycoords="axes fraction", **args)

    # Now the total number of detections
    ax.annotate(
        "Total: {}".format(dets[iyear(eyear)]),
        xy=(eyear - 2, dets[iyear(eyear)]),
        xycoords="data",
        va="center",
        ha="right",
        size=16,
    )

    # Now the facilities
    args = {"size": "16"}
    arrowprops = {"arrowstyle": "-|>", "connectionstyle": "arc3", "facecolor": "black"}
    ax.annotate(
        "NRAO 36-foot (1968)",
        xy=(1968, dets[iyear(1968)] + 4),
        xycoords="data",
        xytext=(0, 35),
        textcoords="offset points",
        rotation=90,
        arrowprops=arrowprops,
        va="bottom",
        ha="center",
        **args
    )
    ax.annotate(
        "",
        xy=(1982, dets[iyear(1982)] - 4),
        xycoords="data",
        xytext=(0, -35),
        textcoords="offset points",
        rotation=0,
        arrowprops=arrowprops,
        va="bottom",
        ha="left",
        **args
    )
    ax.annotate(
        "Nobeyama (1982)", xy=(1982, dets[iyear(1982)] - 31), xycoords="data", **args
    )
    ax.annotate(
        "IRAM (1984)",
        xy=(1984, dets[iyear(1984)] + 4),
        xycoords="data",
        xytext=(0, 35),
        textcoords="offset points",
        rotation=90,
        arrowprops=arrowprops,
        va="bottom",
        ha="center",
        **args
    )
    ax.annotate(
        "GBT (2001)",
        xy=(2001, dets[iyear(2001)] + 4),
        xycoords="data",
        xytext=(0, 35),
        textcoords="offset points",
        rotation=90,
        arrowprops=arrowprops,
        va="bottom",
        ha="center",
        **args
    )
    ax.annotate(
        "ALMA (2011)",
        xy=(2011, dets[iyear(2011)] - 4),
        xycoords="data",
        xytext=(0, -35),
        textcoords="offset points",
        rotation=90,
        arrowprops=arrowprops,
        va="top",
        ha="center",
        **args
    )

    # show the plot
    plt.show()

    # write out the figure
    plt.savefig(
        filename if filename is not None else "cumulative_detections.pdf",
        format="pdf",
        transparent=True,
        bbox_inches="tight",
    )


def cumu_det_natoms_plot(mol_list=None, syear=None, eyear=None, filename=None):

    """
    Makes a plot of the cumulative detections (sorted by atoms) by year using the molecules in mol_list.
    Defaults to using all molecules in the database, but can be passed a subset of molecules as a list of Molecule objects in mol_list.
    The start year and end year are defined by default to be the earliest year in the list and the current year + 20 (to give room for labels),
    but these can be overriden by setting syear and eyear to the desired years (integers).

    Default filename is 'cumulative_by_atoms.pdf', but that can also be overriden.

    """

    # If a list wasn't specified, default to all molecules
    if mol_list is None:
        mol_list = all_molecules

    # Close an old figure if it exists and initialize a new figure
    plt.close("Cumulative Detections By Atoms")
    fig = plt.figure(num="Cumulative Detections By Atoms", figsize=(10, 8))
    plt.ion()

    # set some font defaults
    fontparams = {"size": 24, "family": "sans-serif", "sans-serif": ["Helvetica"]}
    plt.rc("font", **fontparams)
    plt.rc("mathtext", fontset="stixsans")

    # get the starting and ending years, if they aren't set by the user

    if syear is None:
        # grab the earliest year in the list of molecules
        syear = min([x.year for x in mol_list])

    if eyear is None:
        # grab the current year
        eyear = date.today().year

    # make the x-axis array of years
    years = np.arange(syear, eyear + 1)

    # now we make a dictionary of the different traces we're gonna want, loop through the list, and populate an array of detections for that number of atoms or special case
    dets_dict = {}
    for natoms in range(2, 14):
        # fill it with an empte. array
        dets_dict[natoms] = np.copy(years) * 0
        # loop through the years and the list and add everything up.  There's gotta be a better way to do this, but who cares, its fast enough.
        for x in range(len(dets_dict[natoms])):
            i = 0
            for mol in mol_list:
                if mol.year < years[x] + 1 and mol.natoms == natoms:
                    i += 1
            dets_dict[natoms][x] = i

    # do the fullerenes and pahs
    dets_dict["fullerenes"] = np.copy(years) * 0
    for x in range(len(dets_dict["fullerenes"])):
        i = 0
        for mol in mol_list:
            if mol.year < years[x] + 1 and mol.fullerene is True:
                i += 1
        dets_dict["fullerenes"][x] = i

    dets_dict["pahs"] = np.copy(years) * 0
    for x in range(len(dets_dict["pahs"])):
        i = 0
        for mol in mol_list:
            if mol.year < years[x] + 1 and mol.pah is True:
                i += 1
        dets_dict["pahs"][x] = i

    # load up an axis
    ax = fig.add_subplot(111)

    # label the axes
    plt.xlabel("Year")
    plt.ylabel("Cumulative Number of Detected Molecules")

    # customize tick marks
    ax.tick_params(axis="x", which="both", direction="in", length=15, width=1)
    ax.tick_params(axis="y", which="both", direction="in", length=15, width=1)
    ax.yaxis.set_ticks_position("both")
    ax.xaxis.set_ticks_position("both")
    ax.set_xticks([1940, 1960, 1980, 2000, 2020])
    ax.set_xlim([syear, eyear])

    # plot the traces
    ax.plot(years, dets_dict[2], color="#000000")
    ax.plot(years, dets_dict[3], color="#800000")
    ax.plot(years, dets_dict[4], color="#f032e6")
    ax.plot(years, dets_dict[5], color="#9A6324")
    ax.plot(years, dets_dict[6], color="dodgerblue")
    ax.plot(years, dets_dict[7], color="#e6194B")
    ax.plot(years, dets_dict[8], color="#469990")
    ax.plot(years, dets_dict[9], color="#f58231")
    ax.plot(years, dets_dict[10], color="#42d4f4")
    ax.plot(years, dets_dict[11], color="#ffe119")
    ax.plot(years, dets_dict[12], color="#3cb44b")
    ax.plot(years, dets_dict[13], color="#e6beff")
    ax.plot(years, dets_dict["fullerenes"], color="#000075")
    ax.plot(years, dets_dict["pahs"], color="#aaffc3")

    # add the annotations

    start_y = 0.92
    step_y = 0.04

    start_x = 0.10

    ax.annotate(
        r"\textbf{2 atoms}",
        xy=(start_x, start_y - (0 * step_y)),
        xycoords="axes fraction",
        size=16,
        color="#000000",
        va="top",
        ha="left",
    )
    ax.annotate(
        r"\textbf{3 atoms}",
        xy=(start_x, start_y - (1 * step_y)),
        xycoords="axes fraction",
        size=16,
        color="#800000",
        va="top",
        ha="left",
    )
    ax.annotate(
        r"\textbf{4 atoms}",
        xy=(start_x, start_y - (2 * step_y)),
        xycoords="axes fraction",
        size=16,
        color="#f032e6",
        va="top",
        ha="left",
    )
    ax.annotate(
        r"\textbf{5 atoms}",
        xy=(start_x, start_y - (3 * step_y)),
        xycoords="axes fraction",
        size=16,
        color="#9A6324",
        va="top",
        ha="left",
    )
    ax.annotate(
        r"\textbf{6 atoms}",
        xy=(start_x, start_y - (4 * step_y)),
        xycoords="axes fraction",
        size=16,
        color="dodgerblue",
        va="top",
        ha="left",
    )
    ax.annotate(
        r"\textbf{7 atoms}",
        xy=(start_x, start_y - (5 * step_y)),
        xycoords="axes fraction",
        size=16,
        color="#e6194B",
        va="top",
        ha="left",
    )
    ax.annotate(
        r"\textbf{8 atoms}",
        xy=(start_x, start_y - (6 * step_y)),
        xycoords="axes fraction",
        size=16,
        color="#469990",
        va="top",
        ha="left",
    )
    ax.annotate(
        r"\textbf{9 atoms}",
        xy=(start_x, start_y - (7 * step_y)),
        xycoords="axes fraction",
        size=16,
        color="#f58231",
        va="top",
        ha="left",
    )
    ax.annotate(
        r"\textbf{10 atoms}",
        xy=(start_x, start_y - (8 * step_y)),
        xycoords="axes fraction",
        size=16,
        color="#42d4f4",
        va="top",
        ha="left",
    )
    ax.annotate(
        r"\textbf{11 atoms}",
        xy=(start_x, start_y - (9 * step_y)),
        xycoords="axes fraction",
        size=16,
        color="#ffe119",
        va="top",
        ha="left",
    )
    ax.annotate(
        r"\textbf{12 atoms}",
        xy=(start_x, start_y - (10 * step_y)),
        xycoords="axes fraction",
        size=16,
        color="#3cb44b",
        va="top",
        ha="left",
    )
    ax.annotate(
        r"\textbf{13+ atoms}",
        xy=(start_x, start_y - (11 * step_y)),
        xycoords="axes fraction",
        size=16,
        color="#e6beff",
        va="top",
        ha="left",
    )

    ax.annotate(
        r"\textbf{Fullerenes}",
        xy=(start_x, start_y - (12 * step_y)),
        xycoords="axes fraction",
        size=16,
        color="#000075",
        va="top",
        ha="left",
    )
    ax.annotate(
        r"\textbf{PAHs}",
        xy=(start_x, start_y - (13 * step_y)),
        xycoords="axes fraction",
        size=16,
        color="#aaffc3",
        va="top",
        ha="left",
    )

    # show the plot
    plt.show()

    # write out the figure
    plt.savefig(
        filename if filename is not None else "cumulative_by_atoms.pdf",
        format="pdf",
        transparent=True,
        bbox_inches="tight",
    )


def det_per_year_per_atom(mol_list=None, filename=None):

    """
    Makes a plot of the average number of detections per year (y) for a molecule with (x) atoms, starting in the year they were first detected.
    Has the ability to plot PAHs and fullerenes, but doesn't.
    Defaults to using all molecules in the database, but can be passed a subset of molecules as a list of Molecule objects in mol_list.
    Default filename is 'rate_by_atoms.pdf', but that can also be overriden.
    """

    # If a list wasn't specified, default to all molecules
    if mol_list is None:
        mol_list = all_molecules

    # Close an old figure if it exists and initialize a new figure
    plt.close("Detects Per Year Per Atom")
    fig = plt.figure(num="Detects Per Year Per Atom", figsize=(10, 8))
    plt.ion()

    # set some font defaults
    fontparams = {"size": 24, "family": "sans-serif", "sans-serif": ["Helvetica"]}
    plt.rc("font", **fontparams)
    plt.rc("mathtext", fontset="stixsans")

    # We're going to cheat later on and lump PAHs and Fullerenes together, but these need fake natoms.
    # They'll always be +2 and +4, respectively, beyond the maximum of 'normal' molecules, and we'll
    # leave +1 blank for visual separation.  For now, though, time to make a list and loop through our molecules.
    natoms = np.arange(2, 18)

    # and an array of the average number of detections
    avg_dets = np.copy(natoms) * 0.0
    n_dets = np.copy(natoms) * 0.0

    # get the years that are being spanned
    eyear = max([x.year for x in mol_list])
    for x in range(len(natoms)):
        i = 0
        years = []
        if natoms[x] < 14:
            for mol in mol_list:
                if mol.natoms == natoms[x]:
                    i += 1
                    years.append(mol.year)
        elif natoms[x] == 15:
            for mol in mol_list:
                if mol.pah is True:
                    i += 1
                    years.append(mol.year)
        elif natoms[x] == 17:
            for mol in mol_list:
                if mol.fullerene is True:
                    i += 1
                    years.append(mol.year)
        if i == 0:
            avg_dets[x] = np.NaN
        else:
            avg_dets[x] = i / (eyear - min(years) + 1)
            n_dets[x] = i
    n_dets[n_dets == 0] = np.nan

    # load up an axis
    ax = fig.add_subplot(111)

    # label the axes
    plt.xlabel("Number of Atoms")
    plt.ylabel("Detections/Year*")

    # customize tick marks
    ax.tick_params(axis="x", which="both", direction="in", length=15, width=1)
    ax.tick_params(axis="y", which="both", direction="in", length=15, width=1)
    ax.yaxis.set_ticks_position("both")
    ax.xaxis.set_ticks_position("both")
    ax.set_xticks([2, 4, 6, 8, 10, 12, 15, 17])
    ax.set_xticklabels(["2", "4", "6", "8", "10", "12", "PAHs", "Fullerenes"])
    ax.set_xlim([1, 12.8])
    ax.set_ylim([0, 1.0])

    # plot the data
    mfc = "dodgerblue"
    mec = change_color(mfc, 1.5)
    sizes = np.copy(n_dets) * 200 / np.nanmin(n_dets)
    ax.scatter(natoms, avg_dets, marker="o", c=mfc, edgecolors=mec, s=sizes, alpha=0.9)

    # add some annotations
    ax.annotate(
        "*Since year of first detection",
        xy=(0.35, 0.9),
        xycoords="axes fraction",
        ha="left",
    )
    ax.annotate(
        r"\noindent Marker size proportional to\\ total \# of detections",
        xy=(0.05, 0.1),
        xycoords="axes fraction",
        ha="left",
    )
    for label in ax.get_xmajorticklabels():
        if label._text == "Fullerenes" or label._text == "PAHs":
            label.set_rotation(-45)
    fig.tight_layout()

    # show the plot
    plt.show()

    # write out the figure
    plt.savefig(
        filename if filename is not None else "rate_by_atoms.pdf",
        format="pdf",
        transparent=True,
        bbox_inches="tight",
    )

    return


def facility_shares(mol_list=None, telescopes_list=None, filename=None):

    """
    Makes a plot of the percentage share of yearly detections that a facility contributed over its operational lifetime for the top 9 facilities.
    Defaults to using all molecules in the database, but can be passed a subset of molecules as a list of Molecule objects in mol_list.
    Defaults to using all telescopes in the database, but can be passed a subset of telescopes as a list of Telescope objects in mol_list.
    Default filename is 'facility_shares.pdf', but that can also be overriden.
    """

    # If a list wasn't specified, default to all molecules
    if mol_list is None:
        mol_list = all_molecules
    # If a list wasn't specified, default to all telescopes
    if telescopes_list is None:
        telescopes_list = all_telescopes

    # Close an old figure if it exists and initialize a new figure
    plt.close("Facility Shares")
    fig, axs = plt.subplots(3, 3, num="Facility Shares")
    plt.ion()

    # set some font defaults
    fontparams = {"size": 14, "family": "sans-serif", "sans-serif": ["Helvetica"]}
    plt.rc("font", **fontparams)
    plt.rc("mathtext", fontset="stixsans")

    # we need to generate the data now, which we'll store in a dictionary for each telescope.
    # Each entry will be [syear,eyear,ndetects,ntotal,shortname] for the start year, end year,
    # number of detections, and total number of detections over those years
    detects = [x.ndetects for x in telescopes_list]
    detects.sort(reverse=True)
    # determine the 9 most productive facilities on an absolute scale
    detects = detects[:9]
    # get the minimum number of detections of these
    min_allowed = min(detects)

    my_dict = {}
    for scope in telescopes_list:
        # number of detections
        ndetects = scope.ndetects
        # determine if this telescope has enough detections to qualify
        # if it does, just pass and do the rest
        if ndetects >= min_allowed:
            pass
        # otherwise, move on to the next
        else:
            continue
        # year it was built
        syear = scope.built
        # year it was decommissioned or this year if it's still in operation
        eyear = (
            scope.decommissioned
            if scope.decommissioned is not None
            else date.today().year
        )

        # now we go get the total number of detections in that time
        ntotal = 0
        for mol in mol_list:
            if mol.year <= eyear and mol.year >= syear:
                ntotal += 1
        my_dict[scope.shortname] = [syear, eyear, ndetects, ntotal, scope.shortname]

    my_list = [my_dict[x] for x in my_dict]
    my_list.sort(key=lambda x: x[2] / x[3], reverse=True)

    idx = 0
    for i in range(3):
        for j in range(3)[::]:
            fracs = [
                my_list[idx][2] / my_list[idx][3],
                1.0 - my_list[idx][2] / my_list[idx][3],
            ]
            label = (
                r"\textbf"
                + "{"
                + "{}".format(int(my_list[idx][2] / my_list[idx][3] * 100))
                + r"\%"
            )

            if my_list[idx][1] == date.today().year:
                color = "dodgerblue"
            else:
                color = "#F87070"

            if idx > 5:
                slices, labels = axs[i, j].pie(
                    fracs,
                    colors=[color, "#F5F6FF"],
                    wedgeprops={"linewidth": 1.0, "edgecolor": "black"},
                )

                axs[i, j].annotate(
                    my_list[idx][-1],
                    xy=(0.5, 1.08),
                    xycoords="axes fraction",
                    ha="center",
                    size=12,
                )
                axs[i, j].annotate(
                    "{} - {}".format(my_list[idx][0], my_list[idx][1]),
                    xy=(0.5, 0.95),
                    xycoords="axes fraction",
                    ha="center",
                    size=10,
                )

                kw = dict(arrowprops=dict(arrowstyle="-"), zorder=0, va="center")
                ang = (slices[0].theta2 - slices[0].theta1) / 2.0 + slices[0].theta1
                y = np.sin(np.deg2rad(ang))
                x = np.cos(np.deg2rad(ang))
                horizontalalignment = {-1: "right", 1: "left"}[int(np.sign(x))]
                connectionstyle = "angle,angleA=0,angleB={}".format(ang)
                kw["arrowprops"].update({"connectionstyle": connectionstyle})
                axs[i, j].annotate(
                    label,
                    xy=(x, y),
                    xytext=(1.1 * np.sign(x), y),
                    horizontalalignment=horizontalalignment,
                    size=10,
                    **kw
                )

            else:
                my_labels = [
                    r"\textbf"
                    + "{"
                    + "{}".format(int(my_list[idx][2] / my_list[idx][3] * 100))
                    + r"\%",
                    "",
                ]

                slices, labels = axs[i, j].pie(
                    fracs,
                    labels=my_labels,
                    colors=[color, "#F5F6FF"],
                    wedgeprops={"linewidth": 1.0, "edgecolor": "black"},
                    labeldistance=0.7,
                )

                axs[i, j].annotate(
                    my_list[idx][-1],
                    xy=(0.5, 1.08),
                    xycoords="axes fraction",
                    ha="center",
                    size=12,
                )
                axs[i, j].annotate(
                    "{} - {}".format(my_list[idx][0], my_list[idx][1]),
                    xy=(0.5, 0.95),
                    xycoords="axes fraction",
                    ha="center",
                    size=10,
                )

                for label in labels:
                    label.set_horizontalalignment("center")
                    label.set_color("white")
                    label.set_fontsize(10)

            # draw a line around the outside of the pie
            # get the center and radius of the pie wedges
            center = slices[0].center
            r = slices[0].r
            # create a new circle with the desired properties
            circle = matplotlib.patches.Circle(
                center, r, fill=False, edgecolor="black", linewidth=1
            )
            # add the circle to the axes
            axs[i, j].add_patch(circle)

            idx += 1

    fig.tight_layout()
    fig.subplots_adjust(wspace=-0.5, hspace=0.15)
    plt.show()

    plt.savefig(
        filename if filename is not None else "facility_shares.pdf",
        format="pdf",
        transparent=True,
        bbox_inches="tight",
    )


def cumu_det_facility(mol_list=None, telescopes_list=None, min_detects=10, filename=None):

    """
    Makes a plot of the cumulative number of detections of a facility with time.
    Defaults to using all molecules in the database, but can be passed a subset of molecules as a list of Molecule objects in mol_list.
    Defaults to using all telescopes in the database, but can be passed a subset of telescopes as a list of Telescope objects in mol_list.
    Will only plot facilities with min_detects or greater, where min_detects defaults to 10 but can be overriden with an integer.
    Default filename is 'facility_shares.pdf', but that can also be overriden.
    """

    #If a list wasn't specified, default to all molecules
    if mol_list is None:
        mol_list = all_molecules
    # If a list wasn't specified, default to all telescopes
    if telescopes_list is None:
        telescopes_list = all_telescopes

    # Close an old figure if it exists and initialize a new figure
    plt.close("Detections Per Facility Over Time")
    fig = plt.figure(num="Detections Per Facility Over Time", figsize=(10, 8))
    plt.ion()

    # set some font defaults
    fontparams = {"size": 24, "family": "sans-serif", "sans-serif": ["Helvetica"]}
    plt.rc("font", **fontparams)
    plt.rc("mathtext", fontset="stixsans")
    
    # We're only going to do facilities with 10 or more total detections.
    # set this thing up to kick off a few years before 1968 and run until today

    years = np.arange(1965, date.today().year + 1)
    scopes = [x for x in telescopes_list if x.ndetects >= min_detects]

    my_dict = {}
    for scope in scopes:
        tmp_years = np.copy(years) * 0
        i = 0
        for x in range(len(years)):
            for mol in mol_list:
                if mol.year == years[x] and scope in mol.telescopes:
                    i += 1
            tmp_years[x] = i
        my_dict[scope.shortname] = tmp_years

    # pull up an axis
    ax = fig.add_subplot(111)

    # label the axes
    plt.xlabel("Year")
    plt.ylabel("Cumulative Number of Detected Molecules")

    # customize tick marks
    ax.tick_params(axis="x", which="both", direction="in", length=15, width=1)
    ax.tick_params(axis="y", which="both", direction="in", length=15, width=1)
    ax.yaxis.set_ticks_position("both")
    ax.xaxis.set_ticks_position("both")
    ax.set_xticks([1970, 1980, 1990, 2000, 2010, 2020])

    ax.plot(years, my_dict["GBT"], color="#000000")
    ax.plot(years, my_dict["IRAM"], color="#800000")
    ax.plot(years, my_dict["NRAO 140-ft"], color="#f032e6")
    ax.plot(years, my_dict["NRAO/ARO 12-m"], color="dodgerblue")
    ax.plot(years, my_dict["NRAO 36-ft"], color="#e6194B")
    ax.plot(years, my_dict["Nobeyama"], color="#469990")

    ax.annotate(
        r"\textbf{GBT}",
        xy=(date.today().year + 1, my_dict["GBT"][-1]),
        xycoords="data",
        size=16,
        color="#000000",
        va="center",
        zorder=100,
    )
    ax.annotate(
        r"\textbf{IRAM}",
        xy=(date.today().year + 1, my_dict["IRAM"][-1]),
        xycoords="data",
        size=16,
        color="#800000",
        va="center",
        zorder=100,
    )
    ax.annotate(
        r"\textbf{NRAO 140-ft}",
        xy=(date.today().year + 1, my_dict["NRAO 140-ft"][-1]),
        xycoords="data",
        size=16,
        color="#f032e6",
        va="top",
        zorder=100,
    )
    ax.annotate(
        r"\textbf{NRAO/ARO 12-m}",
        xy=(date.today().year + 1, my_dict["NRAO/ARO 12-m"][-1]),
        xycoords="data",
        size=16,
        color="dodgerblue",
        va="center",
        zorder=100,
    )
    ax.annotate(
        r"\textbf{NRAO 36-ft}",
        xy=(date.today().year + 1, my_dict["NRAO 36-ft"][-1]),
        xycoords="data",
        size=16,
        color="#e6194B",
        va="center",
        zorder=100,
    )
    ax.annotate(
        r"\textbf{Nobeyama}",
        xy=(date.today().year + 1, my_dict["Nobeyama"][-1]),
        xycoords="data",
        size=16,
        color="#469990",
        va="bottom",
        zorder=100,
    )

    ax.set_xlim([1965, date.today().year + 19])
    ax.set_ylim(0, max([my_dict[x][-1] for x in my_dict]) + 5)

    # do linear fits to the data for the ranges we care about for each facility:
    # get some year indicies for years we care about
    def iyear(x):
        return np.argwhere(years == x)[0][0]

    trendGBT = (
        np.polynomial.polynomial.Polynomial.fit(
            years[iyear(GBT.built) :], my_dict["GBT"][iyear(GBT.built) :], 1
        )
        .convert()
        .coef[1]
    )
    ax.annotate(
        "{:.1f}/yr".format(trendGBT),
        xy=(2014, 7.5),
        xycoords="data",
        size=16,
        color="#000000",
        ha="center",
    )
    ax.annotate(
        "{}{: <7}".format(GBT.built, " - "),
        xy=(2014, 4.5),
        xycoords="data",
        size=16,
        color="#000000",
        ha="right",
    )

    trendIRAMold = (
        np.polynomial.polynomial.Polynomial.fit(
            years[iyear(IRAM30.built) : iyear(2006)],
            my_dict["IRAM"][iyear(IRAM30.built) : iyear(2006)],
            1,
        )
        .convert()
        .coef[1]
    )
    ax.annotate(
        "{:.1f}/yr".format(trendIRAMold),
        xy=(1990, 21),
        xycoords="data",
        size=16,
        color="#800000",
        ha="center",
    )
    ax.annotate(
        "{} - 2006".format(IRAM30.built),
        xy=(1990, 18),
        xycoords="data",
        size=16,
        color="#800000",
        ha="center",
    )

    trendIRAMnew = (
        np.polynomial.polynomial.Polynomial.fit(
            years[iyear(2006) :], my_dict["IRAM"][iyear(2006) :], 1
        )
        .convert()
        .coef[1]
    )
    ax.annotate(
        "{:.1f}/yr".format(trendIRAMnew),
        xy=(2020, 45),
        xycoords="data",
        size=16,
        color="#800000",
        ha="center",
    )
    ax.annotate(
        "2006{: <7}".format(" - "),
        xy=(2020, 42),
        xycoords="data",
        size=16,
        color="#800000",
        ha="right",
    )

    trend140 = (
        np.polynomial.polynomial.Polynomial.fit(
            years[iyear(NRAO140.built) : 1993],
            my_dict["NRAO 140-ft"][iyear(NRAO140.built) : 1993],
            1,
        )
        .convert()
        .coef[1]
    )
    ax.annotate(
        "{:.1f}/yr".format(trend140),
        xy=(1980, 12.5),
        xycoords="data",
        size=16,
        color="#f032e6",
        ha="center",
    )
    ax.annotate(
        "{} - 1993".format(NRAO140.built),
        xy=(1980, 9.5),
        xycoords="data",
        size=16,
        color="#f032e6",
        ha="center",
    )

    trend12 = (
        np.polynomial.polynomial.Polynomial.fit(
            years[iyear(NRAOARO12.built) :],
            my_dict["NRAO/ARO 12-m"][iyear(NRAOARO12.built) :],
            1,
        )
        .convert()
        .coef[1]
    )
    ax.annotate(
        "{:.1f}/yr".format(trend12),
        xy=(2014.8, 30.2),
        xycoords="data",
        size=16,
        color="dodgerblue",
        ha="center",
    )
    ax.annotate(
        "{}{: <7}".format(NRAOARO12.built, " - "),
        xy=(2014.8, 27.2),
        xycoords="data",
        size=16,
        color="dodgerblue",
        ha="right",
    )

    trend36 = (
        np.polynomial.polynomial.Polynomial.fit(
            years[iyear(NRAO36.built) : 1985],
            my_dict["NRAO 36-ft"][iyear(NRAO36.built) : 1985],
            1,
        )
        .convert()
        .coef[1]
    )
    ax.annotate(
        "{:.1f}/yr".format(trend36),
        xy=(1975, 34),
        xycoords="data",
        size=16,
        color="#e6194B",
        ha="center",
    )
    ax.annotate(
        "{} - 1985".format(NRAO36.built),
        xy=(1975, 31),
        xycoords="data",
        size=16,
        color="#e6194B",
        ha="center",
    )

    trendNobeyama = (
        np.polynomial.polynomial.Polynomial.fit(
            years[iyear(Nobeyama45.built) :],
            my_dict["Nobeyama"][iyear(Nobeyama45.built) :],
            1,
        )
        .convert()
        .coef[1]
    )
    ax.annotate(
        "{:.1f}/yr".format(trendNobeyama),
        xy=(2012, 19),
        xycoords="data",
        size=16,
        color="#469990",
        ha="center",
    )
    ax.annotate(
        "{}{: <7}".format(Nobeyama45.built, " - "),
        xy=(2012, 16),
        xycoords="data",
        size=16,
        color="#469990",
        ha="right",
    )

    plt.show()

    plt.savefig(
        filename if filename is not None else "scopes_by_year.pdf",
        format="pdf",
        transparent=True,
        bbox_inches="tight",
    )

# def periodic_heatmap(mol_list=None, filename=None):


#     '''
#     Makes a periodic table heat map.
#     Defaults to using all molecules in the database, but can be passed a subset of molecules as a list of Molecule objects in mol_list.
#     Default filename is 'periodic_heatmap.pdf', but that can also be overriden.
#     '''

#     # If a list wasn't specified, default to all molecules
#     if mol_list is None:
#         mol_list = all_molecules

#     #Close an old figure if it exists and initialize a new figure   
#     plt.close('Periodic Heatmap')    
#     plt.figure(num='Periodic Heatmap',figsize=(20,9.5))    
#     plt.ion()
    
#     #set some font defaults   
#     fontparams = {'size':24, 'family':'sans-serif','sans-serif':['Helvetica']}   
#     plt.rc('font',**fontparams)
#     plt.rc('mathtext', fontset='stixsans')	
#     #test
#     #get a list of elements involved - we'll use the masses list from molecules.py, as that will be update to date
#     els = list(masses.keys())
#     #open a new dictionary to store detections from these
#     census = {}
    
#     #loop through, creating the dictionary entry if needed
#     for mol in mol_list:
#         for el in els:    
#             if el in mol.atoms:
#                 if mol.atoms[el] > 0:
#                     if el in census:    
#                         census[el] += 1
#                     else:
#                         census[el] = 1
                
#     maxdets = max([census[x] for x in census])   
#     map_colors = list(Color("#f1fb53").range_to(Color("#f00707"),maxdets))
#     map_colors = [str(x) for x in map_colors]		
    
#     #Dictionary for the periodic table   
#     elements = {
#         'H'		:	[periodictable.H,1,6.9],
#         'He'	:	[periodictable.He,18,6.9],
#         'Li'	:	[periodictable.Li,1,5.75],
#         'Be'	:	[periodictable.Be,2,5.75],
#         'B'		:	[periodictable.B,13,5.75],
#         'C'		:	[periodictable.C,14,5.75],
#         'N'		:	[periodictable.N,15,5.75],
#         'O'		:	[periodictable.O,16,5.75],
#         'F'		:	[periodictable.F,17,5.75],
#         'Ne'	:	[periodictable.Ne,18,5.75],
#         'Na'	:	[periodictable.Na,1,4.6],
#         'Mg'	:	[periodictable.Mg,2,4.6],
#         'Al'	:	[periodictable.Al,13,4.6],
#         'Si'	:	[periodictable.Si,14,4.6],
#         'P'		:	[periodictable.P,15,4.6],
#         'S'		:	[periodictable.S,16,4.6],
#         'Cl'	:	[periodictable.Cl,17,4.6],
#         'Ar'	:	[periodictable.Ar,18,4.6],
#         'K'		:	[periodictable.K,1,3.45],
#         'Ca'	:	[periodictable.Ca,2,3.45],
#         'Sc'	:	[periodictable.Sc,3,3.45],
#         'Ti'	:	[periodictable.Ti,4,3.45],
#         'V'		:	[periodictable.V,5,3.45],
#         'Cr'	:	[periodictable.Cr,6,3.45],
#         'Mn'	:	[periodictable.Mn,7,3.45],
#         'Fe'	:	[periodictable.Fe,8,3.45],
#         'Co'	:	[periodictable.Co,9,3.45],
#         'Ni'	:	[periodictable.Ni,10,3.45],
#         'Cu'	:	[periodictable.Cu,11,3.45],
#         'Zn'	:	[periodictable.Zn,12,3.45],
#         'Ga'	:	[periodictable.Ga,13,3.45],
#         'Ge'	:	[periodictable.Ge,14,3.45],
#         'As'	:	[periodictable.As,15,3.45],
#         'Se'	:	[periodictable.Se,16,3.45],
#         'Br'	:	[periodictable.Br,17,3.45],
#         'Kr'	:	[periodictable.Kr,18,3.45],
#         'Rb'	:	[periodictable.Rb,1,2.3],
#         'Sr'	:	[periodictable.Sr,2,2.3],
#         'Y'		:	[periodictable.Y,3,2.3],
#         'Zr'	:	[periodictable.Zr,4,2.3],
#         'Nb'	:	[periodictable.Nb,5,2.3],
#         'Mo'	:	[periodictable.Mo,6,2.3],
#         'Tc'	:	[periodictable.Tc,7,2.3],
#         'Ru'	:	[periodictable.Ru,8,2.3],
#         'Rh'	:	[periodictable.Rh,9,2.3],
#         'Pd'	:	[periodictable.Pd,10,2.3],
#         'Ag'	:	[periodictable.Ag,11,2.3],
#         'Cd'	:	[periodictable.Cd,12,2.3],
#         'In'	:	[periodictable.In,13,2.3],
#         'Sn'	:	[periodictable.Sn,14,2.3],
#         'Sb'	:	[periodictable.Sb,15,2.3],
#         'Te'	:	[periodictable.Te,16,2.3],
#         'I'		:	[periodictable.I,17,2.3],
#         'Xe'	:	[periodictable.Xe,18,2.3],
#         'Cs'	:	[periodictable.Cs,1,1.15],
#         'Ba'	:	[periodictable.Ba,2,1.15],
#         'Hf'	:	[periodictable.Hf,4,1.15],
#         'Ta'	:	[periodictable.Ta,5,1.15],
#         'W'		:	[periodictable.W,6,1.15],
#         'Re'	:	[periodictable.Re,7,1.15],
#         'Os'	:	[periodictable.Os,8,1.15],
#         'Ir'	:	[periodictable.Ir,9,1.15],
#         'Pt'	:	[periodictable.Pt,10,1.15],
#         'Au'	:	[periodictable.Au,11,1.15],
#         'Hg'	:	[periodictable.Hg,12,1.15],
#         'Tl'	:	[periodictable.Tl,13,1.15],
#         'Pb'	:	[periodictable.Pb,14,1.15],
#         'Bi'	:	[periodictable.Bi,15,1.15],
#         'Po'	:	[periodictable.Po,16,1.15],
#         'At'	:	[periodictable.At,17,1.15],
#         'Rn'	:	[periodictable.Rn,18,1.15],
#         'Fr'	:	[periodictable.Fr,1,0.],
#         'Ra'	:	[periodictable.Ra,2,0.],
#         'Rf'	:	[periodictable.Rf,4,0.],
#         'Db'	:	[periodictable.Db,5,0.],
#         'Sg'	:	[periodictable.Sg,6,0.],
#         'Bh'	:	[periodictable.Bh,7,0.],
#         'Hs'	:	[periodictable.Hs,8,0.],
#         'Mt'	:	[periodictable.Mt,9,0.],
#         'Ds'	:	[periodictable.Ds,10,0.],
#         'Rg'	:	[periodictable.Rg,11,0.],
#         'Cn'	:	[periodictable.Cn,12,0.],
#         'Nh'	:	[periodictable.Nh,13,0.],
#         'Fl'	:	[periodictable.Fl,14,0.],
#         'Mc'	:	[periodictable.Mc,15,0.],
#         'Lv'	:	[periodictable.Lv,16,0.],
#         'Ts'	:	[periodictable.Ts,17,0.],
#         'Og'	:	[periodictable.Og,18,0.],
#         }

#     #load up an axis 
#     ax = plt.axes([0,0,1,1])  
#     ax.set_xlim([0,18])
#     ax.set_ylim([0,8])
    
#     for el in elements:      
#         x = elements[el][1]-1
#         y = elements[el][2]      
#         sym = r'\textbf{' + elements[el][0].symbol + '}'
#         num = elements[el][0].number
#         mass = elements[el][0].mass
#         name = elements[el][0].name.capitalize()
        
#         this_color = 'white'        
#         if el in census:      
#             this_color = map_colors[census[el]-1] 
#             ndets = r'\textbf{' + str(census[el]) + '}'
#             ax.annotate(ndets,xy=(x+0.8,y+.95),xycoords='data',size=14,color='black',ha='right',va='top')

#         rect = patches.Rectangle((x,y),0.9,1.05,linewidth=1,edgecolor='black',facecolor=this_color,alpha=0.5)
#         ax.add_patch(rect)
#         ax.annotate(num,xy=(x+0.1,y+.95),xycoords='data',size=14,color='black',ha='left',va='top')
#         ax.annotate(sym,xy=(x+0.1,y+.70),xycoords='data',size=20,color='black',ha='left',va='top')
#         ax.annotate(mass,xy=(x+0.1,y+.42),xycoords='data',size=8,color='black',ha='left',va='top')
#         ax.annotate(name,xy=(x+0.1,y+.29),xycoords='data',size=8,color='black',ha='left',va='top')
            
    
#     plt.axis('equal')
#     plt.axis('off')	
#     plt.show()
    
#     #write out the figure 
#     plt.savefig(filename if filename is not None else 'periodic_heatmap.pdf',format='pdf',transparent=True,bbox_inches='tight',pad_inches=0)
    
#     #the bit below crops off extra white space.  This only works on Macs with the TexLive pdfcrop utility installed.  Comment out if not desired.
#     os.system('pdfcrop --margins -0 periodic_heatmap.pdf periodic_heatmap.pdf')
    
#     return	