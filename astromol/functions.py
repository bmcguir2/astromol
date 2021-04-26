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
matplotlib.rc("text.latex", preamble=r"\usepackage{cmbright}\usepackage[version=4]{mhchem}")

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

    trend1968 = np.polynomial.polynomial.Polynomial.fit(years[iyear(1968) :], dets[iyear(1968) :], 1).convert().coef[1]
    trend2005 = np.polynomial.polynomial.Polynomial.fit(years[iyear(2005) :], dets[iyear(2005) :], 1).convert().coef[1]

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
    ax.annotate("Nobeyama (1982)", xy=(1982, dets[iyear(1982)] - 31), xycoords="data", **args)
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
        eyear = scope.decommissioned if scope.decommissioned is not None else date.today().year

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
            label = r"\textbf" + "{" + "{}".format(int(my_list[idx][2] / my_list[idx][3] * 100)) + r"\%"

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
                    r"\textbf" + "{" + "{}".format(int(my_list[idx][2] / my_list[idx][3] * 100)) + r"\%",
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
            circle = matplotlib.patches.Circle(center, r, fill=False, edgecolor="black", linewidth=1)
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

    # If a list wasn't specified, default to all molecules
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
        np.polynomial.polynomial.Polynomial.fit(years[iyear(GBT.built) :], my_dict["GBT"][iyear(GBT.built) :], 1)
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
        np.polynomial.polynomial.Polynomial.fit(years[iyear(2006) :], my_dict["IRAM"][iyear(2006) :], 1)
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


def periodic_heatmap(mol_list=None, filename=None):

    """
    Makes a periodic table heat map.
    Defaults to using all molecules in the database, but can be passed a subset of molecules as a list of Molecule objects in mol_list.
    Default filename is 'periodic_heatmap.pdf', but that can also be overriden.
    """

    # If a list wasn't specified, default to all molecules
    if mol_list is None:
        mol_list = all_molecules

    # Close an old figure if it exists and initialize a new figure
    plt.close("Periodic Heatmap")
    plt.figure(num="Periodic Heatmap", figsize=(20, 9.5))
    plt.ion()

    # set some font defaults
    fontparams = {"size": 24, "family": "sans-serif", "sans-serif": ["Helvetica"]}
    plt.rc("font", **fontparams)
    plt.rc("mathtext", fontset="stixsans")
    # test
    # get a list of elements involved - we'll use the masses list from molecules.py, as that will be update to date
    els = list(masses.keys())
    # open a new dictionary to store detections from these
    census = {}

    # loop through, creating the dictionary entry if needed
    for mol in mol_list:
        for el in els:
            if el in mol.atoms:
                if mol.atoms[el] > 0:
                    if el in census:
                        census[el] += 1
                    else:
                        census[el] = 1

    maxdets = max([census[x] for x in census])
    map_colors = list(Color("#f1fb53").range_to(Color("#f00707"), maxdets))
    map_colors = [str(x) for x in map_colors]

    # Dictionary for the periodic table
    elements = {
        "H": [periodictable.H, 1, 6.9],
        "He": [periodictable.He, 18, 6.9],
        "Li": [periodictable.Li, 1, 5.75],
        "Be": [periodictable.Be, 2, 5.75],
        "B": [periodictable.B, 13, 5.75],
        "C": [periodictable.C, 14, 5.75],
        "N": [periodictable.N, 15, 5.75],
        "O": [periodictable.O, 16, 5.75],
        "F": [periodictable.F, 17, 5.75],
        "Ne": [periodictable.Ne, 18, 5.75],
        "Na": [periodictable.Na, 1, 4.6],
        "Mg": [periodictable.Mg, 2, 4.6],
        "Al": [periodictable.Al, 13, 4.6],
        "Si": [periodictable.Si, 14, 4.6],
        "P": [periodictable.P, 15, 4.6],
        "S": [periodictable.S, 16, 4.6],
        "Cl": [periodictable.Cl, 17, 4.6],
        "Ar": [periodictable.Ar, 18, 4.6],
        "K": [periodictable.K, 1, 3.45],
        "Ca": [periodictable.Ca, 2, 3.45],
        "Sc": [periodictable.Sc, 3, 3.45],
        "Ti": [periodictable.Ti, 4, 3.45],
        "V": [periodictable.V, 5, 3.45],
        "Cr": [periodictable.Cr, 6, 3.45],
        "Mn": [periodictable.Mn, 7, 3.45],
        "Fe": [periodictable.Fe, 8, 3.45],
        "Co": [periodictable.Co, 9, 3.45],
        "Ni": [periodictable.Ni, 10, 3.45],
        "Cu": [periodictable.Cu, 11, 3.45],
        "Zn": [periodictable.Zn, 12, 3.45],
        "Ga": [periodictable.Ga, 13, 3.45],
        "Ge": [periodictable.Ge, 14, 3.45],
        "As": [periodictable.As, 15, 3.45],
        "Se": [periodictable.Se, 16, 3.45],
        "Br": [periodictable.Br, 17, 3.45],
        "Kr": [periodictable.Kr, 18, 3.45],
        "Rb": [periodictable.Rb, 1, 2.3],
        "Sr": [periodictable.Sr, 2, 2.3],
        "Y": [periodictable.Y, 3, 2.3],
        "Zr": [periodictable.Zr, 4, 2.3],
        "Nb": [periodictable.Nb, 5, 2.3],
        "Mo": [periodictable.Mo, 6, 2.3],
        "Tc": [periodictable.Tc, 7, 2.3],
        "Ru": [periodictable.Ru, 8, 2.3],
        "Rh": [periodictable.Rh, 9, 2.3],
        "Pd": [periodictable.Pd, 10, 2.3],
        "Ag": [periodictable.Ag, 11, 2.3],
        "Cd": [periodictable.Cd, 12, 2.3],
        "In": [periodictable.In, 13, 2.3],
        "Sn": [periodictable.Sn, 14, 2.3],
        "Sb": [periodictable.Sb, 15, 2.3],
        "Te": [periodictable.Te, 16, 2.3],
        "I": [periodictable.I, 17, 2.3],
        "Xe": [periodictable.Xe, 18, 2.3],
        "Cs": [periodictable.Cs, 1, 1.15],
        "Ba": [periodictable.Ba, 2, 1.15],
        "Hf": [periodictable.Hf, 4, 1.15],
        "Ta": [periodictable.Ta, 5, 1.15],
        "W": [periodictable.W, 6, 1.15],
        "Re": [periodictable.Re, 7, 1.15],
        "Os": [periodictable.Os, 8, 1.15],
        "Ir": [periodictable.Ir, 9, 1.15],
        "Pt": [periodictable.Pt, 10, 1.15],
        "Au": [periodictable.Au, 11, 1.15],
        "Hg": [periodictable.Hg, 12, 1.15],
        "Tl": [periodictable.Tl, 13, 1.15],
        "Pb": [periodictable.Pb, 14, 1.15],
        "Bi": [periodictable.Bi, 15, 1.15],
        "Po": [periodictable.Po, 16, 1.15],
        "At": [periodictable.At, 17, 1.15],
        "Rn": [periodictable.Rn, 18, 1.15],
        "Fr": [periodictable.Fr, 1, 0.0],
        "Ra": [periodictable.Ra, 2, 0.0],
        "Rf": [periodictable.Rf, 4, 0.0],
        "Db": [periodictable.Db, 5, 0.0],
        "Sg": [periodictable.Sg, 6, 0.0],
        "Bh": [periodictable.Bh, 7, 0.0],
        "Hs": [periodictable.Hs, 8, 0.0],
        "Mt": [periodictable.Mt, 9, 0.0],
        "Ds": [periodictable.Ds, 10, 0.0],
        "Rg": [periodictable.Rg, 11, 0.0],
        "Cn": [periodictable.Cn, 12, 0.0],
        "Nh": [periodictable.Nh, 13, 0.0],
        "Fl": [periodictable.Fl, 14, 0.0],
        "Mc": [periodictable.Mc, 15, 0.0],
        "Lv": [periodictable.Lv, 16, 0.0],
        "Ts": [periodictable.Ts, 17, 0.0],
        "Og": [periodictable.Og, 18, 0.0],
    }

    # load up an axis
    ax = plt.axes([0, 0, 1, 1])
    ax.set_xlim([0, 18])
    ax.set_ylim([0, 8])

    for el in elements:
        x = elements[el][1] - 1
        y = elements[el][2]
        sym = r"\textbf{" + elements[el][0].symbol + "}"
        num = elements[el][0].number
        mass = elements[el][0].mass
        name = elements[el][0].name.capitalize()

        this_color = "white"
        if el in census:
            this_color = map_colors[census[el] - 1]
            ndets = r"\textbf{" + str(census[el]) + "}"
            ax.annotate(
                ndets,
                xy=(x + 0.8, y + 0.95),
                xycoords="data",
                size=14,
                color="black",
                ha="right",
                va="top",
            )

        rect = patches.Rectangle(
            (x, y),
            0.9,
            1.05,
            linewidth=1,
            edgecolor="black",
            facecolor=this_color,
            alpha=0.5,
        )
        ax.add_patch(rect)
        ax.annotate(
            num,
            xy=(x + 0.1, y + 0.95),
            xycoords="data",
            size=14,
            color="black",
            ha="left",
            va="top",
        )
        ax.annotate(
            sym,
            xy=(x + 0.1, y + 0.70),
            xycoords="data",
            size=20,
            color="black",
            ha="left",
            va="top",
        )
        ax.annotate(
            mass,
            xy=(x + 0.1, y + 0.42),
            xycoords="data",
            size=8,
            color="black",
            ha="left",
            va="top",
        )
        ax.annotate(
            name,
            xy=(x + 0.1, y + 0.29),
            xycoords="data",
            size=8,
            color="black",
            ha="left",
            va="top",
        )

    plt.axis("equal")
    plt.axis("off")
    plt.show()

    # write out the figure
    plt.savefig(
        filename if filename is not None else "periodic_heatmap.pdf",
        format="pdf",
        transparent=True,
        bbox_inches="tight",
        pad_inches=0,
    )

    # the bit below crops off extra white space.  This only works on Macs with the TexLive pdfcrop utility installed.
    # Comment out if not desired.
    os.system("pdfcrop --margins -0 periodic_heatmap.pdf periodic_heatmap.pdf")


def mass_by_wavelength(mol_list=None, bw=0.5, filename=None):
    """
    Makes a KDE plot of detections at each wavelength vs mass.
    Defaults to using all molecules in the database, but can be passed a subset of molecules as a list of Molecule objects in mol_list.
    Defaults to a bandwidth of 0.5, which can be overriden by any float.
    Default filename is 'mass_by_wavelengths_kde.pdf', but that can also be overriden.
    """

    # If a list wasn't specified, default to all molecules
    if mol_list is None:
        mol_list = all_molecules

    # gather the data
    my_dict = {
        "UV": [],
        "Vis": [],
        "IR": [],
        "sub-mm": [],
        "mm": [],
        "cm": [],
        "UV-Vis": [],
    }

    for x in mol_list:
        for y in my_dict:
            if y in x.wavelengths:
                my_dict[y].append(x.mass)
    for x in my_dict["UV"]:
        my_dict["UV-Vis"].append(x)
    for x in my_dict["Vis"]:
        my_dict["UV-Vis"].append(x)

    # get the plot set up
    plt.close("Detections at Wavelengths by Mass")
    fig = plt.figure(num="Detections at Wavelengths by Mass", figsize=(10, 8))
    plt.ion()

    # set some font defaults
    fontparams = {"size": 24, "family": "sans-serif", "sans-serif": ["Helvetica"]}
    plt.rc("font", **fontparams)
    plt.rc("mathtext", fontset="stixsans")

    # load up an axis
    ax = fig.add_subplot(111)
    ax.tick_params(axis="x", which="both", direction="in", length=5, width=1)
    ax.tick_params(axis="y", which="both", direction="in", length=5, width=1)

    # label the axes
    plt.xlabel("Atomic Mass (amu)")
    plt.ylabel("Probability Density Estimate")

    # Do the estimates and plot them
    xvals = np.arange(0, 160)

    density_cm = gkde(my_dict["cm"])
    density_cm.covariance_factor = lambda: bw
    density_cm._compute_covariance()

    density_mm = gkde(my_dict["mm"])
    density_mm.covariance_factor = lambda: bw
    density_mm._compute_covariance()

    density_submm = gkde(my_dict["sub-mm"])
    density_submm.covariance_factor = lambda: bw
    density_submm._compute_covariance()

    density_IR = gkde(my_dict["IR"])
    density_IR.covariance_factor = lambda: bw
    density_IR._compute_covariance()

    density_UV = gkde(my_dict["UV-Vis"])
    density_UV.covariance_factor = lambda: bw
    density_UV._compute_covariance()

    ax.plot(xvals, density_cm(xvals), color="dodgerblue")
    ax.fill_between(xvals, density_cm(xvals), 0, facecolor="dodgerblue", alpha=0.25, zorder=4)
    ax.annotate(
        "{}".format(len(my_dict["cm"])),
        xy=(75, 0.01),
        xycoords="data",
        ha="left",
        va="bottom",
        color="dodgerblue",
    )

    max_mm = max([x for x in my_dict["mm"] if x < 160])
    ax.plot(xvals[: max_mm + 1], density_mm(xvals[: max_mm + 1]), color="darkorange")
    ax.fill_between(
        xvals[: max_mm + 1],
        density_mm(xvals[: max_mm + 1]),
        0,
        facecolor="darkorange",
        alpha=0.25,
    )
    ax.annotate(
        "{}".format(len(my_dict["mm"])),
        xy=(53.5, 0.022),
        xycoords="data",
        ha="left",
        va="bottom",
        color="darkorange",
    )

    ax.plot(xvals, density_submm(xvals), color="forestgreen")
    ax.fill_between(xvals, density_submm(xvals), 0, facecolor="forestgreen", alpha=0.25)
    ax.annotate(
        "{}".format(len(my_dict["sub-mm"])),
        xy=(40, 0.038),
        xycoords="data",
        ha="left",
        va="bottom",
        color="forestgreen",
    )

    max_IR = max([x for x in my_dict["IR"] if x < 160])
    ax.plot(xvals[: max_IR + 1], density_IR(xvals[: max_IR + 1]), color="black")
    ax.fill_between(
        xvals[: max_IR + 1],
        density_IR(xvals[: max_IR + 1]),
        0,
        facecolor="black",
        alpha=0.25,
        zorder=5,
    )
    ax.annotate(
        "{}".format(len(my_dict["IR"])),
        xy=(68, 0.0025),
        xycoords="data",
        ha="left",
        va="bottom",
        color="black",
    )

    max_UV = max(my_dict["UV-Vis"])
    ax.plot(xvals[: max_UV + 1], density_UV(xvals[: max_UV + 1]), color="violet")
    ax.fill_between(
        xvals[: max_UV + 1],
        density_UV(xvals[: max_UV + 1]),
        0,
        facecolor="violet",
        alpha=0.25,
    )
    ax.annotate(
        "{}".format(len(my_dict["UV-Vis"])),
        xy=(16.5, 0.057),
        xycoords="data",
        ha="left",
        va="bottom",
        color="violet",
    )

    ax.annotate(
        r"\underline{Detection Wavelengths}",
        xy=(158, 0.06),
        xycoords="data",
        color="black",
        ha="right",
        va="top",
    )
    ax.annotate(
        "centimeter",
        xy=(158, 0.056),
        xycoords="data",
        color="dodgerblue",
        ha="right",
        va="top",
    )
    ax.annotate(
        "millimeter",
        xy=(158, 0.052),
        xycoords="data",
        color="darkorange",
        ha="right",
        va="top",
    )
    ax.annotate(
        "sub-millimeter",
        xy=(158, 0.048),
        xycoords="data",
        color="forestgreen",
        ha="right",
        va="top",
    )
    ax.annotate(
        "infrared",
        xy=(158, 0.044),
        xycoords="data",
        color="black",
        ha="right",
        va="top",
    )
    ax.annotate(
        "visible/ultraviolet",
        xy=(158, 0.04),
        xycoords="data",
        color="violet",
        ha="right",
        va="top",
    )

    plt.show()

    plt.savefig(
        filename if filename is not None else "mass_by_wavelengths_kde.pdf",
        format="pdf",
        transparent=True,
        bbox_inches="tight",
        pad_inches=0,
    )

    return


def mols_waves_by_atoms(mol_list=None, bw=0.5, filename=None):
    """
    Makes six histogram plots of molecules detected in each wavelength range by number of atoms, excepting fullerenes.
    For plots with sufficient datapoints, we'll do a KDE plot with bandwidth defaulting to 0.5 that can be overridden.
    Defaults to using all molecules in the database, but can be passed a subset of molecules as a list of Molecule objects in mol_list.
    Default filename is 'mols_waves_by_atoms.pdf', but that can also be overriden.
    """

    # If a list wasn't specified, default to all molecules
    if mol_list is None:
        mol_list = all_molecules

    # Close an old figure if it exists and initialize a new figure
    plt.close("Molecules Detected in Each Wavelength by Number of Atoms")
    plt.figure(num="Molecules Detected in Each Wavelength by Number of Atoms", figsize=(10, 8))
    plt.ion()

    # set some font defaults
    fontparams = {"size": 18, "family": "sans-serif", "sans-serif": ["Helvetica"]}
    plt.rc("font", **fontparams)
    plt.rc("mathtext", fontset="stixsans")

    # gather the data
    my_dict = {
        "UV": [],
        "Vis": [],
        "IR": [],
        "sub-mm": [],
        "mm": [],
        "cm": [],
    }

    max_n = []
    for x in mol_list:
        for y in my_dict:
            if y in x.wavelengths:
                if x.fullerene is True:
                    continue
                else:
                    my_dict[y].append(x.natoms)
                    max_n.append(x.natoms)

    max_n = max(max_n)

    # pull up some axes to plot on
    ax1 = plt.subplot(231)
    ax2 = plt.subplot(232)
    ax3 = plt.subplot(233)
    ax4 = plt.subplot(234)
    ax5 = plt.subplot(235)
    ax6 = plt.subplot(236)

    xvals = np.arange(0, max_n + 1, 0.5)

    ax1.tick_params(axis="x", which="both", direction="in", length=5, width=1)
    ax1.tick_params(axis="y", which="both", direction="in", length=5, width=1)
    ax2.tick_params(axis="x", which="both", direction="in", length=5, width=1)
    ax2.tick_params(axis="y", which="both", direction="in", length=5, width=1)
    ax3.tick_params(axis="x", which="both", direction="in", length=5, width=1)
    ax3.tick_params(axis="y", which="both", direction="in", length=5, width=1)
    ax4.tick_params(axis="x", which="both", direction="in", length=5, width=1)
    ax4.tick_params(axis="y", which="both", direction="in", length=5, width=1)
    ax5.tick_params(axis="x", which="both", direction="in", length=5, width=1)
    ax5.tick_params(axis="y", which="both", direction="in", length=5, width=1)
    ax6.tick_params(axis="x", which="both", direction="in", length=5, width=1)
    ax6.tick_params(axis="y", which="both", direction="in", length=5, width=1)

    density_cm = gkde(my_dict["cm"])
    density_cm.covariance_factor = lambda: bw
    density_cm._compute_covariance()
    ax1.plot(xvals, density_cm(xvals))
    ax1.fill_between(xvals, density_cm(xvals), 0, facecolor="dodgerblue", alpha=0.25)
    ax1.annotate("cm", xy=(0.95, 0.96), xycoords="axes fraction", ha="right", va="top", size=24)
    ax1.set_ylim([0, 0.65])
    ax1.set_xticks([0, 5, 10, 15, 20])
    ax1.set_xticklabels([])
    ax1.set_ylabel("Probability Density Estimate")

    density_mm = gkde(my_dict["mm"])
    density_mm.covariance_factor = lambda: bw
    density_mm._compute_covariance()
    ax2.plot(xvals, density_mm(xvals))
    ax2.fill_between(xvals, density_mm(xvals), 0, facecolor="dodgerblue", alpha=0.25)
    ax2.annotate("mm", xy=(0.95, 0.96), xycoords="axes fraction", ha="right", va="top", size=24)
    ax2.set_ylim([0, 0.65])
    ax2.set_xticks([0, 5, 10, 15, 20])
    ax2.set_xticklabels([])

    density_submm = gkde(my_dict["sub-mm"])
    density_submm.covariance_factor = lambda: bw
    density_submm._compute_covariance()
    ax3.plot(xvals, density_submm(xvals))
    ax3.fill_between(xvals, density_submm(xvals), 0, facecolor="dodgerblue", alpha=0.25)
    ax3.annotate("sub-mm", xy=(0.95, 0.96), xycoords="axes fraction", ha="right", va="top", size=24)
    ax3.set_ylim([0, 0.65])
    ax3.set_xticks([0, 5, 10, 15, 20])
    ax3.set_xticklabels([])

    density_IR = gkde(my_dict["IR"])
    density_IR.covariance_factor = lambda: bw
    density_IR._compute_covariance()
    ax4.plot(xvals, density_IR(xvals))
    ax4.fill_between(xvals, density_IR(xvals), 0, facecolor="dodgerblue", alpha=0.25)
    ax4.annotate("IR", xy=(0.95, 0.96), xycoords="axes fraction", ha="right", va="top", size=24)
    ax4.set_ylim([0, 0.65])
    ax4.set_xticks([0, 5, 10, 15, 20])
    ax4.set_xlabel(r"\# of Atoms")
    ax4.set_ylabel("Probability Density Estimate")

    ax5.hist(my_dict["Vis"], bins=[0.5, 1.5, 2.5, 3.5, 4.5])
    ax5.annotate("Vis", xy=(0.95, 0.96), xycoords="axes fraction", ha="right", va="top", size=24)
    ax5.set_xlim([0, 20])
    ax5.set_ylim([0, 8])
    ax5.set_xticks([0, 5, 10, 15, 20])
    ax5.set_xticklabels([])
    ax5.set_ylabel(r"\# of Detected Molecules")

    ax6.hist(my_dict["UV"], bins=[0.5, 1.5, 2.5, 3.5, 4.5])
    ax6.annotate("UV", xy=(0.95, 0.96), xycoords="axes fraction", ha="right", va="top", size=24)
    ax6.set_xlim([0, 20])
    ax6.set_ylim([0, 8])
    ax6.set_xticks([0, 5, 10, 15, 20])
    ax6.set_xticklabels([])

    plt.tight_layout()
    plt.show()
    plt.savefig(
        filename if filename is not None else "mols_waves_by_atoms.pdf",
        format="pdf",
        transparent=True,
        bbox_inches="tight",
        pad_inches=0,
    )

    return


def du_histogram(mol_list=None, filename=None):
    """
    Makes a histogram of the degree of unsaturation of molecules containing only H, O, N, C, Cl, or F.
    Defaults to using all molecules in the database, but can be passed a subset of molecules as a list of Molecule objects in mol_list.
    Default filename is 'du_histogram.pdf', but that can also be overriden.
    """

    # If a list wasn't specified, default to all molecules
    if mol_list is None:
        mol_list = all_molecules

    # Close an old figure if it exists and initialize a new figure
    plt.close("Degree of Unsaturation Histogram")
    plt.figure(num="Degree of Unsaturation Histogram", figsize=(10, 8))
    plt.ion()

    # set some font defaults
    fontparams = {"size": 24, "family": "sans-serif", "sans-serif": ["Helvetica"]}
    plt.rc("font", **fontparams)
    plt.rc("mathtext", fontset="stixsans")

    # gather the data
    dus = []
    for x in mol_list:
        # make sure we're working with molecules containing only H, D, N, C, Cl, F, O, or S.
        atoms = []
        for atom in ["H", "D", "N", "C", "Cl", "F", "S", "O"]:
            if atom in x.atoms:
                atoms.append(x.atoms[atom])
        if np.sum(atoms) == x.natoms and x.du is not None and x.fullerene is not True:
            dus.append(x.du)

    # set up a plot
    ax = plt.subplot(111)

    ax.tick_params(axis="x", which="both", direction="in", length=5, width=1)
    ax.tick_params(axis="y", which="both", direction="in", length=5, width=1)

    ax.yaxis.set_ticks_position("both")
    ax.xaxis.set_ticks_position("both")

    plt.xlabel("Degree of Unsaturation")
    plt.ylabel(r"\# of Detected Molecules")

    bins = np.arange(-0.25, 12.5, 0.5)
    (n, bins, _) = ax.hist(dus, bins=bins, facecolor="dodgerblue", alpha=0.25)
    ax.hist(dus, bins=bins, edgecolor="royalblue", linewidth=1.5, fill=False)
    ax.set_ylim(0, 35)

    ax.annotate(
        r"\ce{CH4}, \ce{CH3OH}, \ce{CH3Cl}, ...",
        xy=(0, n[0] + 1),
        xycoords="data",
        rotation=90,
        size=16,
        ha="center",
        va="bottom",
    )
    ax.annotate(r"\ce{HC11N}", xy=(12, n[12 * 2] + 1), xycoords="data", rotation=90, size=16, ha="center", va="bottom")

    plt.tight_layout()
    plt.show()

    plt.savefig(
        filename if filename is not None else "du_histogram.pdf", format="pdf", transparent=True, bbox_inches="tight"
    )

    return


def type_pie_chart(mol_list=None, filename=None):
    """
    Makes a pie chart of the fraction of interstellar molecules that are neutral, radical, cation, cyclic, pahs, fullerenes, or anions
    Defaults to using all molecules in the database, but can be passed a subset of molecules as a list of Molecule objects in mol_list.
    Default filename is 'type_pie_chart.pdf', but that can also be overriden.
    """

    # If a list wasn't specified, default to all molecules
    if mol_list is None:
        mol_list = all_molecules

    # Close an old figure if it exists and initialize a new figure
    plt.close("Type Pie Chart")
    plt.figure(num="Type Pie Chart", figsize=(10, 8))
    plt.ion()

    # set some font defaults
    fontparams = {"size": 24, "family": "sans-serif", "sans-serif": ["Helvetica"]}
    plt.rc("font", **fontparams)
    plt.rc("mathtext", fontset="stixsans")

    # gather the data
    my_dict = {
        "Neutral": [0],
        "Radical": [0],
        "Cation": [0],
        "Cyclic": [0],
        "Anion": [0],
        "Fullerene": [0],
        "PAH": [0],
    }

    for mol in mol_list:
        if mol.neutral is True:
            my_dict["Neutral"][0] += 1
        if mol.radical is True:
            my_dict["Radical"][0] += 1
        if mol.cation is True:
            my_dict["Cation"][0] += 1
        if mol.cyclic is True:
            my_dict["Cyclic"][0] += 1
        if mol.anion is True:
            my_dict["Anion"][0] += 1
        if mol.fullerene is True:
            my_dict["Fullerene"][0] += 1
        if mol.pah is True:
            my_dict["PAH"][0] += 1

    nmols = len(mol_list)
    for type in my_dict:
        my_dict[type].append(my_dict[type][0] / nmols)

    labels = ["Neutral", "Radical", "Cation", "Cyclic", "Anion", "Fullerene", "PAH"]
    fracs = [my_dict[x][1] for x in labels]

    # set up a plot
    ax = plt.subplot(111)
    size = 0.1

    def getshift(x):
        return -(90 - (360 - 360 * x) / 2)

    ax.pie(
        [fracs[0], 1.0 - fracs[0]],
        colors=["dodgerblue", "#EEEEEE"],
        radius=1,
        startangle=getshift(fracs[0]),
        wedgeprops=dict(width=size, edgecolor="w", linewidth=1),
    )
    ax.pie(
        [fracs[1], 1.0 - fracs[1]],
        colors=["darkorange", "#EEEEEE"],
        radius=1 - size - 0.02,
        startangle=getshift(fracs[1]),
        wedgeprops=dict(width=size, edgecolor="w", linewidth=1),
    )
    ax.pie(
        [fracs[2], 1.0 - fracs[2]],
        colors=["forestgreen", "#EEEEEE"],
        radius=1 - 2 * size - 0.04,
        startangle=getshift(fracs[2]),
        wedgeprops=dict(width=size, edgecolor="w", linewidth=1),
    )
    ax.pie(
        [fracs[3], 1.0 - fracs[3]],
        colors=["violet", "#EEEEEE"],
        radius=1 - 3 * size - 0.06,
        startangle=getshift(fracs[3]),
        wedgeprops=dict(width=size, edgecolor="w", linewidth=1),
    )
    ax.pie(
        [fracs[4], 1.0 - fracs[4]],
        colors=["red", "#EEEEEE"],
        radius=1 - 4 * size - 0.08,
        startangle=getshift(fracs[4]),
        wedgeprops=dict(width=size, edgecolor="w", linewidth=1),
    )
    ax.pie(
        [fracs[5], 1.0 - fracs[5]],
        colors=["goldenrod", "#EEEEEE"],
        radius=1 - 5 * size - 0.1,
        startangle=getshift(fracs[5]),
        wedgeprops=dict(width=size, edgecolor="w", linewidth=1),
    )
    ax.pie(
        [fracs[6], 1.0 - fracs[6]],
        colors=["royalblue", "#EEEEEE"],
        radius=1 - 6 * size - 0.12,
        startangle=getshift(fracs[6]),
        wedgeprops=dict(width=size, edgecolor="w", linewidth=1),
    )

    ax.annotate(
        r"\textbf{Neutrals}", xy=(0.5, 0.11), xycoords="axes fraction", color="dodgerblue", ha="center", size=14
    )
    ax.annotate(
        r"\textbf{Radicals}", xy=(0.5, 0.16), xycoords="axes fraction", color="darkorange", ha="center", size=14
    )
    ax.annotate(
        r"\textbf{Cations}", xy=(0.5, 0.205), xycoords="axes fraction", color="forestgreen", ha="center", size=14
    )
    ax.annotate(r"\textbf{Cyclics}", xy=(0.5, 0.255), xycoords="axes fraction", color="violet", ha="center", size=14)
    ax.annotate(r"\textbf{Anions}", xy=(0.5, 0.305), xycoords="axes fraction", color="red", ha="center", size=14)
    ax.annotate(
        r"\textbf{Fullerenes}", xy=(0.5, 0.3575), xycoords="axes fraction", color="goldenrod", ha="center", size=14
    )
    ax.annotate(r"\textbf{PAHs}", xy=(0.5, 0.40), xycoords="axes fraction", color="royalblue", ha="center", size=14)

    percents = [r"\textbf{" + "{:.1f}".format((x * 100)) + r"}\%" for x in fracs]

    start = 0.585
    shift = 0.0485
    ax.annotate(
        percents[0],
        xy=(start + 6 * shift, 0.5),
        xycoords="axes fraction",
        color="white",
        ha="center",
        va="center",
        size=12,
        rotation=-90,
    )
    ax.annotate(
        percents[1],
        xy=(start + 5 * shift, 0.5),
        xycoords="axes fraction",
        color="darkorange",
        ha="center",
        va="center",
        size=12,
        rotation=-90,
    )
    ax.annotate(
        percents[2],
        xy=(start + 4 * shift, 0.5),
        xycoords="axes fraction",
        color="forestgreen",
        ha="center",
        va="center",
        size=12,
        rotation=-90,
    )
    ax.annotate(
        percents[3],
        xy=(start + 3 * shift, 0.5),
        xycoords="axes fraction",
        color="violet",
        ha="center",
        va="center",
        size=12,
        rotation=-90,
    )
    ax.annotate(
        percents[4],
        xy=(start + 2 * shift, 0.5),
        xycoords="axes fraction",
        color="red",
        ha="center",
        va="center",
        size=12,
        rotation=-90,
    )
    ax.annotate(
        percents[5],
        xy=(start + 1 * shift, 0.5),
        xycoords="axes fraction",
        color="goldenrod",
        ha="center",
        va="center",
        size=12,
        rotation=-90,
    )
    ax.annotate(
        percents[6],
        xy=(start + 0 * shift, 0.5),
        xycoords="axes fraction",
        color="royalblue",
        ha="center",
        va="center",
        size=12,
        rotation=-90,
    )

    plt.tight_layout()
    plt.show()

    plt.savefig(
        filename if filename is not None else "type_pie_chart.pdf", format="pdf", transparent=True, bbox_inches="tight"
    )  # ,pad_inches=-.65)

    return


def source_pie_chart(mol_list=None, filename=None):
    """
    Makes a pie chart of the fraction of interstellar molecules detected in carbon stars, dark clouds, los clouds, star forming regions, and other types of sources
    Defaults to using all molecules in the database, but can be passed a subset of molecules as a list of Molecule objects in mol_list.
    Default filename is 'source_pie_chart.pdf', but that can also be overriden.
    """

    # If a list wasn't specified, default to all molecules
    if mol_list is None:
        mol_list = all_molecules

    # Close an old figure if it exists and initialize a new figure
    plt.close("Source Pie Chart")
    plt.figure(num="Source Pie Chart", figsize=(10, 8))
    plt.ion()

    # set some font defaults
    fontparams = {"size": 24, "family": "sans-serif", "sans-serif": ["Helvetica"]}
    plt.rc("font", **fontparams)
    plt.rc("mathtext", fontset="stixsans")

    # gather the data
    my_dict = {
        "Carbon Star": [0],
        "Dark Cloud": [0],
        "LOS Cloud": [0],
        "SFR": [0],
        "Other": [0],
    }

    # we have to be a little careful here, because for a given source, there can be two SFRs listed, and we only want to credit it once
    for mol in mol_list:

        # we'll make a dictionary here to flag if we've credited things already
        credit_dict = {
            "Carbon Star": False,
            "Dark Cloud": False,
            "LOS Cloud": False,
            "SFR": False,
            "Other": False,
        }

        for source in mol.sources:
            if source.type in my_dict:
                if credit_dict[source.type] is False:
                    my_dict[source.type][0] += 1
                    credit_dict[source.type] = True
            else:
                if credit_dict["Other"] is False:
                    my_dict["Other"][0] += 1
                    credit_dict["Other"] = True

    nmols = len(mol_list)
    for type in my_dict:
        my_dict[type].append(my_dict[type][0] / nmols)

    labels = ["SFR", "Carbon Star", "Dark Cloud", "Other", "LOS Cloud"]
    fracs = [my_dict[x][1] for x in labels]

    # set up a plot
    ax = plt.subplot(111)
    size = 0.1

    def getshift(x):
        return -(90 - (360 - 360 * x) / 2)

    ax.pie(
        [fracs[0], 1.0 - fracs[0]],
        colors=["dodgerblue", "#EEEEEE"],
        radius=1,
        startangle=getshift(fracs[0]),
        wedgeprops=dict(width=size, edgecolor="w", linewidth=1),
    )
    ax.pie(
        [fracs[1], 1.0 - fracs[1]],
        colors=["darkorange", "#EEEEEE"],
        radius=1 - size - 0.02,
        startangle=getshift(fracs[1]),
        wedgeprops=dict(width=size, edgecolor="w", linewidth=1),
    )
    ax.pie(
        [fracs[2], 1.0 - fracs[2]],
        colors=["forestgreen", "#EEEEEE"],
        radius=1 - 2 * size - 0.04,
        startangle=getshift(fracs[2]),
        wedgeprops=dict(width=size, edgecolor="w", linewidth=1),
    )
    ax.pie(
        [fracs[3], 1.0 - fracs[3]],
        colors=["violet", "#EEEEEE"],
        radius=1 - 3 * size - 0.06,
        startangle=getshift(fracs[3]),
        wedgeprops=dict(width=size, edgecolor="w", linewidth=1),
    )
    ax.pie(
        [fracs[4], 1.0 - fracs[4]],
        colors=["red", "#EEEEEE"],
        radius=1 - 4 * size - 0.08,
        startangle=getshift(fracs[4]),
        wedgeprops=dict(width=size, edgecolor="w", linewidth=1),
    )

    ax.annotate(r"\textbf{SFR}", xy=(0.5, 0.11), xycoords="axes fraction", color="dodgerblue", ha="center", size=14)
    ax.annotate(
        r"\textbf{Carbon Star}", xy=(0.5, 0.16), xycoords="axes fraction", color="darkorange", ha="center", size=14
    )
    ax.annotate(
        r"\textbf{Dark Cloud}", xy=(0.5, 0.205), xycoords="axes fraction", color="forestgreen", ha="center", size=14
    )
    ax.annotate(r"\textbf{Other}", xy=(0.5, 0.255), xycoords="axes fraction", color="violet", ha="center", size=14)
    ax.annotate(r"\textbf{LOS Cloud}", xy=(0.5, 0.305), xycoords="axes fraction", color="red", ha="center", size=14)

    percents = [r"\textbf{" + "{:.1f}".format((x * 100)) + r"}\%" for x in fracs]

    ax.annotate(
        percents[4],
        xy=(0.51, 0.68),
        xycoords="axes fraction",
        color="white",
        ha="center",
        size=12,
    )
    ax.annotate(
        percents[3],
        xy=(0.51, 0.725),
        xycoords="axes fraction",
        color="white",
        ha="center",
        size=12,
    )
    ax.annotate(
        percents[2],
        xy=(0.51, 0.775),
        xycoords="axes fraction",
        color="white",
        ha="center",
        size=12,
    )
    ax.annotate(
        percents[1],
        xy=(0.51, 0.825),
        xycoords="axes fraction",
        color="white",
        ha="center",
        size=12,
    )
    ax.annotate(
        percents[0],
        xy=(0.51, 0.87),
        xycoords="axes fraction",
        color="white",
        ha="center",
        size=12,
    )

    plt.tight_layout()
    plt.show()

    plt.savefig(
        filename if filename is not None else "source_pie_chart.pdf",
        format="pdf",
        transparent=True,
        bbox_inches="tight",
        pad_inches=-0.65,
    )

    return


def indiv_source_pie_chart(mol_list=None, filename=None):

    """
    Makes a pie chart of the fraction of interstellar molecules detected in IRC+10216, TMC-1, Orion, and Sgr
    Defaults to using all molecules in the database, but can be passed a subset of molecules as a list of Molecule objects in mol_list.
    Default filename is 'indiv_source_pie_chart.pdf', but that can also be overriden.
    """

    # If a list wasn't specified, default to all molecules
    if mol_list is None:
        mol_list = all_molecules

    # Close an old figure if it exists and initialize a new figure
    plt.close("Individual Source Pie Chart")
    plt.figure(num="Individual Source Pie Chart", figsize=(10, 8))
    plt.ion()

    # set some font defaults
    fontparams = {"size": 24, "family": "sans-serif", "sans-serif": ["Helvetica"]}
    plt.rc("font", **fontparams)
    plt.rc("mathtext", fontset="stixsans")

    # gather the data
    my_dict = {
        "IRC+10216": [0],
        "TMC-1": [0],
        "Orion": [0],
        "Sgr B2": [0],
        "Other": [0],
    }

    for mol in mol_list:
        for source in mol.sources:
            other = False
            if source == IRC10216:
                my_dict["IRC+10216"][0] += 1
                other = True
            if source == TMC1:
                my_dict["TMC-1"][0] += 1
                other = True
            if source == Orion:
                my_dict["Orion"][0] += 1
                other = True
            if source == SgrB2:
                my_dict["Sgr B2"][0] += 1
                other = True
            if other is False:
                my_dict["Other"][0] += 1

    nmols = len(mol_list)
    for type in my_dict:
        my_dict[type].append(my_dict[type][0] / nmols)
    labels = ["Other", "Sgr B2", "IRC+10216", "TMC-1", "Orion"]
    fracs = [my_dict[x][1] for x in labels]

    # set up a plot
    ax = plt.subplot(111)
    size = 0.1

    def getshift(x):
        return -(90 - (360 - 360 * x) / 2)

    ax.pie(
        [fracs[0], 1.0 - fracs[0]],
        colors=["dodgerblue", "#EEEEEE"],
        radius=1,
        startangle=getshift(fracs[0]),
        wedgeprops=dict(width=size, edgecolor="w", linewidth=1),
    )
    ax.pie(
        [fracs[1], 1.0 - fracs[1]],
        colors=["darkorange", "#EEEEEE"],
        radius=1 - size - 0.02,
        startangle=getshift(fracs[1]),
        wedgeprops=dict(width=size, edgecolor="w", linewidth=1),
    )
    ax.pie(
        [fracs[2], 1.0 - fracs[2]],
        colors=["forestgreen", "#EEEEEE"],
        radius=1 - 2 * size - 0.04,
        startangle=getshift(fracs[2]),
        wedgeprops=dict(width=size, edgecolor="w", linewidth=1),
    )
    ax.pie(
        [fracs[3], 1.0 - fracs[3]],
        colors=["violet", "#EEEEEE"],
        radius=1 - 3 * size - 0.06,
        startangle=getshift(fracs[3]),
        wedgeprops=dict(width=size, edgecolor="w", linewidth=1),
    )
    ax.pie(
        [fracs[4], 1.0 - fracs[4]],
        colors=["red", "#EEEEEE"],
        radius=1 - 4 * size - 0.08,
        startangle=getshift(fracs[4]),
        wedgeprops=dict(width=size, edgecolor="w", linewidth=1),
    )

    ax.annotate(r"\textbf{Other}", xy=(0.5, 0.11), xycoords="axes fraction", color="dodgerblue", ha="center", size=14)
    ax.annotate(r"\textbf{Sgr B2}", xy=(0.5, 0.16), xycoords="axes fraction", color="darkorange", ha="center", size=14)
    ax.annotate(
        r"\textbf{IRC+10216}", xy=(0.5, 0.205), xycoords="axes fraction", color="forestgreen", ha="center", size=14
    )
    ax.annotate(r"\textbf{TMC-1}", xy=(0.5, 0.255), xycoords="axes fraction", color="violet", ha="center", size=14)
    ax.annotate(r"\textbf{Orion}", xy=(0.5, 0.305), xycoords="axes fraction", color="red", ha="center", size=14)

    percents = [r"\textbf{" + "{:.1f}".format((x * 100)) + r"}\%" for x in fracs]

    ax.annotate(
        percents[4],
        xy=(0.51, 0.68),
        xycoords="axes fraction",
        color="white",
        ha="center",
        size=12,
    )
    ax.annotate(
        percents[3],
        xy=(0.51, 0.725),
        xycoords="axes fraction",
        color="white",
        ha="center",
        size=12,
    )
    ax.annotate(
        percents[2],
        xy=(0.51, 0.775),
        xycoords="axes fraction",
        color="white",
        ha="center",
        size=12,
    )
    ax.annotate(
        percents[1],
        xy=(0.51, 0.825),
        xycoords="axes fraction",
        color="white",
        ha="center",
        size=12,
    )
    ax.annotate(
        percents[0],
        xy=(0.51, 0.87),
        xycoords="axes fraction",
        color="white",
        ha="center",
        size=12,
    )

    plt.tight_layout()
    plt.show()

    plt.savefig(
        filename if filename is not None else "indiv_source_pie_chart.pdf",
        format="pdf",
        transparent=True,
        bbox_inches="tight",
        pad_inches=-0.65,
    )

    return


def mol_type_by_source_type(mol_list=None, filename=None):

    """
    Generates four pie charts, one for each generalized source type, with the wedges for the types of molecules detected first in each type
    Defaults to using all molecules in the database, but can be passed a subset of molecules as a list of Molecule objects in mol_list.
    Default filename is 'mol_type_by_source_type.pdf', but that can also be overriden.
    """

    # If a list wasn't specified, default to all molecules
    if mol_list is None:
        mol_list = all_molecules

    # Close an old figure if it exists and initialize a new figure
    plt.close("Molecule Type by Source Type")
    fig, axs = plt.subplots(2, 2, num="Molecule Type by Source Type", figsize=(15, 12))
    plt.ion()

    # set some font defaults
    fontparams = {"size": 26, "family": "sans-serif", "sans-serif": ["Helvetica"]}
    plt.rc("font", **fontparams)
    plt.rc("mathtext", fontset="stixsans")

    # collect the data
    type_dict = {
        "Carbon Star": {"Anion": 0, "Cation": 0, "Cyclic": 0, "Neutral": 0, "Radical": 0},
        "Dark Cloud": {"Anion": 0, "Cation": 0, "Cyclic": 0, "Neutral": 0, "Radical": 0},
        "LOS Cloud": {"Anion": 0, "Cation": 0, "Cyclic": 0, "Neutral": 0, "Radical": 0},
        "SFR": {"Anion": 0, "Cation": 0, "Cyclic": 0, "Neutral": 0, "Radical": 0},
    }

    for mol in mol_list:
        # we'll make a dictionary here to flag if we've credited things already
        credit_dict = {
            "Carbon Star": False,
            "Dark Cloud": False,
            "LOS Cloud": False,
            "SFR": False,
        }

        for source in mol.sources:
            if source.type in type_dict:
                if credit_dict[source.type] is False:
                    if mol.anion is True:
                        type_dict[source.type]["Anion"] += 1
                    if mol.cation is True:
                        type_dict[source.type]["Cation"] += 1
                    if mol.cyclic is True:
                        type_dict[source.type]["Cyclic"] += 1
                    if mol.neutral is True:
                        type_dict[source.type]["Neutral"] += 1
                    if mol.radical is True:
                        type_dict[source.type]["Radical"] += 1
                    credit_dict[source.type] = True

    # make the pie charts
    # Carbon Stars
    carbon_data = [
        type_dict["Carbon Star"]["Anion"],
        # type_dict['Carbon Star']['Cation'],
        type_dict["Carbon Star"]["Cyclic"],
        type_dict["Carbon Star"]["Neutral"],
        type_dict["Carbon Star"]["Radical"],
    ]

    carbon_colors = [
        "darkorange",
        #'forestgreen',
        "violet",
        "dodgerblue",
        "red",
    ]

    # Dark Clouds
    dark_data = [
        type_dict["Dark Cloud"]["Anion"],
        type_dict["Dark Cloud"]["Cation"],
        type_dict["Dark Cloud"]["Cyclic"],
        type_dict["Dark Cloud"]["Neutral"],
        type_dict["Dark Cloud"]["Radical"],
    ]

    dark_colors = [
        "darkorange",
        "forestgreen",
        "violet",
        "dodgerblue",
        "red",
    ]

    # LOS Clouds
    los_data = [
        # type_dict['LOS Cloud']['Anion'],
        type_dict["LOS Cloud"]["Cation"],
        # type_dict['LOS Cloud']['Cyclic'],
        type_dict["LOS Cloud"]["Neutral"],
        type_dict["LOS Cloud"]["Radical"],
    ]

    los_colors = [
        #'darkorange',
        "forestgreen",
        #'violet',
        "dodgerblue",
        "red",
    ]

    # SFRs
    sfr_data = [
        # type_dict['SFR']['Anion'],
        type_dict["SFR"]["Cation"],
        type_dict["SFR"]["Cyclic"],
        type_dict["SFR"]["Neutral"],
        type_dict["SFR"]["Radical"],
    ]

    sfr_colors = [
        #'darkorange',
        "forestgreen",
        "violet",
        "dodgerblue",
        "red",
    ]

    types = [
        "Anion",
        "Cation",
        "Cyclic",
        "Neutral",
        "Radical",
    ]

    axs[0, 0].pie(
        carbon_data,
        labels=carbon_data,
        colors=carbon_colors,
        labeldistance=0.8,
        wedgeprops={"linewidth": 1.0, "edgecolor": "black", "alpha": 0.5},
    )
    axs[0, 0].annotate(
        r"\textbf{Carbon Stars}", xy=(0.5, 1.0), xycoords="axes fraction", color="black", ha="center", va="top", size=30
    )

    wedges, _ = axs[0, 1].pie(
        dark_data,
        labels=dark_data,
        colors=dark_colors,
        labeldistance=0.8,
        wedgeprops={"linewidth": 1.0, "edgecolor": "black", "alpha": 0.5},
    )
    axs[0, 1].annotate(
        r"\textbf{Dark Clouds}", xy=(0.5, 1.0), xycoords="axes fraction", color="black", ha="center", va="top", size=30
    )

    axs[1, 0].pie(
        los_data,
        labels=los_data,
        colors=los_colors,
        labeldistance=0.8,
        wedgeprops={"linewidth": 1.0, "edgecolor": "black", "alpha": 0.5},
    )
    axs[1, 0].annotate(
        r"\textbf{LOS Clouds}", xy=(0.5, 1.0), xycoords="axes fraction", color="black", ha="center", va="top", size=30
    )

    axs[1, 1].pie(
        sfr_data,
        labels=sfr_data,
        colors=sfr_colors,
        labeldistance=0.8,
        wedgeprops={"linewidth": 1.0, "edgecolor": "black", "alpha": 0.5},
    )
    axs[1, 1].annotate(
        r"\textbf{SFRs}", xy=(0.5, 1.0), xycoords="axes fraction", color="black", ha="center", va="top", size=30
    )

    axs[0, 1].legend(wedges, types, title="Molecule Types", loc="center left", bbox_to_anchor=(1, 0, 0.5, 1))

    fig.tight_layout()
    fig.subplots_adjust(wspace=-0.3, hspace=0)
    plt.show()

    plt.savefig(
        filename if filename is not None else "mol_type_by_source_type.pdf",
        format="pdf",
        transparent=True,
        bbox_inches="tight",
    )

    return


def du_by_source_type(mol_list=None, bw=0.5, filename=None):

    """
    Makes a KDE plot of the dus in each source type
    Defaults to using all molecules in the database, but can be passed a subset of molecules as a list of Molecule objects in mol_list.
    Defaults to a bandwidth of 0.5, which can be overriden by any float.
    Default filename is 'du_by_source_type_kde.pdf', but that can also be overriden.
    """

    # If a list wasn't specified, default to all molecules
    if mol_list is None:
        mol_list = all_molecules

    # gather the data
    my_dict = {
        "Carbon Star": [],
        "Dark Cloud": [],
        "LOS Cloud": [],
        "SFR": [],
    }

    # we have to be a little careful here, because for a given source, there can be two SFRs listed, and we only want to credit it once
    for mol in mol_list:
        # we'll make a dictionary here to flag if we've credited things already
        credit_dict = {
            "Carbon Star": False,
            "Dark Cloud": False,
            "LOS Cloud": False,
            "SFR": False,
        }

        if mol.du is None or mol.fullerene is True:
            continue
        for source in mol.sources:
            if source.type in my_dict:
                if credit_dict[source.type] is False:
                    my_dict[source.type].append(mol.du)
                    credit_dict[source.type] = True

    plt.close("DU by Source Type")
    fig = plt.figure(num="DU by Source Type", figsize=(10, 8))
    plt.ion()

    # set some font defaults
    fontparams = {"size": 24, "family": "sans-serif", "sans-serif": ["Helvetica"]}
    plt.rc("font", **fontparams)
    plt.rc("mathtext", fontset="stixsans")

    # load up an axis
    ax = fig.add_subplot(111)
    ax.tick_params(axis="x", which="both", direction="in", length=5, width=1)
    ax.tick_params(axis="y", which="both", direction="in", length=5, width=1)

    # label the axes
    plt.xlabel("Degree of Unsaturation")
    plt.ylabel("Probably Density Estimate")

    # Do the estimates and plot them
    xvals = np.arange(0, 15, 0.1)

    density_carbon = gkde(my_dict["Carbon Star"])
    density_carbon.covariance_factor = lambda: bw
    density_carbon._compute_covariance()

    density_dark = gkde(my_dict["Dark Cloud"])
    density_dark.covariance_factor = lambda: bw
    density_dark._compute_covariance()

    density_los = gkde(my_dict["LOS Cloud"])
    density_los.covariance_factor = lambda: bw
    density_los._compute_covariance()

    density_sfr = gkde(my_dict["SFR"])
    density_sfr.covariance_factor = lambda: bw
    density_sfr._compute_covariance()

    x_ann = 0.97
    y_ann = 0.96
    y_sep = 0.06

    ax.plot(xvals, density_carbon(xvals), color="darkorange")
    ax.fill_between(xvals, density_carbon(xvals), 0, facecolor="darkorange", alpha=0.25, zorder=4)

    ax.plot(xvals, density_dark(xvals), color="forestgreen")
    ax.fill_between(xvals, density_dark(xvals), 0, facecolor="forestgreen", alpha=0.25, zorder=4)

    ax.plot(xvals, density_los(xvals), color="red")
    ax.fill_between(xvals, density_los(xvals), 0, facecolor="red", alpha=0.25, zorder=4)

    ax.plot(xvals, density_sfr(xvals), color="dodgerblue")
    ax.fill_between(xvals, density_sfr(xvals), 0, facecolor="dodgerblue", alpha=0.25, zorder=4)

    ax.annotate(
        r"\underline{Source Types}",
        xy=(x_ann, y_ann - 0 * y_sep),
        xycoords="axes fraction",
        color="black",
        ha="right",
        va="top",
    )

    ax.annotate("SFR", xy=(x_ann, y_ann - y_sep), xycoords="axes fraction", color="dodgerblue", ha="right", va="top")
    ax.annotate(
        "{}".format(len(my_dict["SFR"])), xy=(2.3, 0.3), xycoords="data", ha="left", va="bottom", color="dodgerblue"
    )

    ax.annotate(
        "Carbon Star", xy=(x_ann, y_ann - 2 * y_sep), xycoords="axes fraction", color="darkorange", ha="right", va="top"
    )
    ax.annotate(
        "{}".format(len(my_dict["Carbon Star"])),
        xy=(6, 0.14),
        xycoords="data",
        ha="left",
        va="bottom",
        color="darkorange",
    )

    ax.annotate(
        "Dark Cloud", xy=(x_ann, y_ann - 3 * y_sep), xycoords="axes fraction", color="forestgreen", ha="right", va="top"
    )
    ax.annotate(
        "{}".format(len(my_dict["Dark Cloud"])),
        xy=(10.9, 0.02),
        xycoords="data",
        ha="left",
        va="bottom",
        color="forestgreen",
    )

    ax.annotate("LOS Cloud", xy=(x_ann, y_ann - 4 * y_sep), xycoords="axes fraction", color="red", ha="right", va="top")
    ax.annotate(
        "{}".format(len(my_dict["LOS Cloud"])), xy=(0.65, 0.36), xycoords="data", ha="left", va="bottom", color="red"
    )

    plt.show()

    plt.savefig(
        filename if filename is not None else "du_by_source_type_kde.pdf",
        format="pdf",
        transparent=True,
        bbox_inches="tight",
        pad_inches=0,
    )

    return


def rel_du_by_source_type(mol_list=None, bw=0.5, filename=None):

    """
    Makes a KDE plot of the relative dus in each source type
    Defaults to using all molecules in the database, but can be passed a subset of molecules as a list of Molecule objects in mol_list.
    Defaults to a bandwidth of 0.5, which can be overriden by any float.
    Default filename is 'relative_du_by_source_type_kde.pdf', but that can also be overriden.
    """

    # If a list wasn't specified, default to all molecules
    if mol_list is None:
        mol_list = all_molecules

    # gather the data
    my_dict = {
        "Carbon Star": [],
        "Dark Cloud": [],
        "LOS Cloud": [],
        "SFR": [],
    }

    # we have to be a little careful here, because for a given source, there can be two SFRs listed, and we only want to credit it once
    for mol in mol_list:
        # we'll make a dictionary here to flag if we've credited things already
        credit_dict = {
            "Carbon Star": False,
            "Dark Cloud": False,
            "LOS Cloud": False,
            "SFR": False,
        }

        if mol.du is None or mol.fullerene is True:
            continue
        for source in mol.sources:
            if source.type in my_dict:
                if credit_dict[source.type] is False:
                    my_dict[source.type].append(mol.du / mol.maxdu)
                    credit_dict[source.type] = True

    plt.close("Relative DU by Source Type")
    _, axs = plt.subplots(2, 2, num="Relative DU by Source Type", figsize=(10, 8))
    plt.ion()

    # set some font defaults
    fontparams = {"size": 18, "family": "sans-serif", "sans-serif": ["Helvetica"]}
    plt.rc("font", **fontparams)
    plt.rc("mathtext", fontset="stixsans")

    # axes ticks and limits
    axs[0, 0].xaxis.set_visible(False)
    axs[0, 0].yaxis.set_visible(False)
    axs[0, 1].xaxis.set_visible(False)
    axs[0, 1].yaxis.set_visible(False)
    axs[1, 1].xaxis.set_visible(False)
    axs[1, 1].yaxis.set_visible(False)

    # label the axes
    axs[1, 0].tick_params(axis="y", which="both", direction="in", length=5, width=1, labelsize=18)
    axs[1, 0].tick_params(axis="x", which="both", direction="in", length=5, width=1, labelsize=18)
    axs[1, 0].set_xlabel("Relative Degree of Unsaturation", fontsize=18)
    axs[1, 0].set_ylabel("Probability Density Estimate", fontsize=18)

    # Do the estimates and plot them
    xvals = np.arange(0, 1, 0.01)

    density_carbon = gkde(my_dict["Carbon Star"])
    density_carbon.covariance_factor = lambda: bw
    density_carbon._compute_covariance()

    density_dark = gkde(my_dict["Dark Cloud"])
    density_dark.covariance_factor = lambda: bw
    density_dark._compute_covariance()

    density_los = gkde(my_dict["LOS Cloud"])
    density_los.covariance_factor = lambda: bw
    density_los._compute_covariance()

    density_sfr = gkde(my_dict["SFR"])
    density_sfr.covariance_factor = lambda: bw
    density_sfr._compute_covariance()

    axs[0, 0].plot(xvals, density_carbon(xvals), color="darkorange")
    axs[0, 0].fill_between(xvals, density_carbon(xvals), 0, facecolor="darkorange", alpha=0.25, zorder=4)
    axs[0, 0].annotate(
        "Carbon Star", xy=[0.04, 0.96], xycoords="axes fraction", ha="left", va="top", size=24, color="darkorange"
    )

    axs[1, 0].plot(xvals, density_dark(xvals), color="forestgreen")
    axs[1, 0].fill_between(xvals, density_dark(xvals), 0, facecolor="forestgreen", alpha=0.25, zorder=4)
    axs[1, 0].annotate(
        "Dark Cloud", xy=[0.04, 0.96], xycoords="axes fraction", ha="left", va="top", size=24, color="forestgreen"
    )

    axs[0, 1].plot(xvals, density_los(xvals), color="red")
    axs[0, 1].fill_between(xvals, density_los(xvals), 0, facecolor="red", alpha=0.25, zorder=4)
    axs[0, 1].annotate(
        "LOS Cloud", xy=[0.04, 0.96], xycoords="axes fraction", ha="left", va="top", size=24, color="red"
    )

    axs[1, 1].plot(xvals, density_sfr(xvals), color="dodgerblue")
    axs[1, 1].fill_between(xvals, density_sfr(xvals), 0, facecolor="dodgerblue", alpha=0.25, zorder=4)
    axs[1, 1].annotate(
        "SFR", xy=[0.04, 0.96], xycoords="axes fraction", ha="left", va="top", size=24, color="dodgerblue"
    )

    plt.subplots_adjust(wspace=0, hspace=0)
    plt.show()

    plt.savefig(
        filename if filename is not None else "relative_du_by_source_type_kde.pdf",
        format="pdf",
        transparent=True,
        bbox_inches="tight",
        pad_inches=0,
    )

    return


def mass_by_source_type(mol_list=None, bw=0.5, filename=None):

    """
    Makes a KDE plot of the masses in each source type
    Defaults to using all molecules in the database, but can be passed a subset of molecules as a list of Molecule objects in mol_list.
    Defaults to a bandwidth of 0.5, which can be overriden by any float.
    Default filename is 'mass_by_source_type_kde.pdf', but that can also be overriden.
    """

    # If a list wasn't specified, default to all molecules
    if mol_list is None:
        mol_list = all_molecules

    # gather the data
    my_dict = {
        "Carbon Star": [],
        "Dark Cloud": [],
        "LOS Cloud": [],
        "SFR": [],
    }

    # we have to be a little careful here, because for a given source, there can be two SFRs listed, and we only want to credit it once
    masses = []
    for mol in mol_list:
        # we'll make a dictionary here to flag if we've credited things already
        credit_dict = {
            "Carbon Star": False,
            "Dark Cloud": False,
            "LOS Cloud": False,
            "SFR": False,
        }

        # drop the fullerenes
        if mol.fullerene is True:
            continue
        # add the mass to the list for axis limit purposes
        masses.append(mol.mass)
        for source in mol.sources:
            if source.type in my_dict:
                if credit_dict[source.type] is False:
                    my_dict[source.type].append(mol.mass)
                    credit_dict[source.type] = True

    plt.close("Mass by Source Type")
    plt.figure(num="Mass by Source Type", figsize=(10, 8))
    plt.ion()

    # set some font defaults
    fontparams = {"size": 24, "family": "sans-serif", "sans-serif": ["Helvetica"]}
    plt.rc("font", **fontparams)
    plt.rc("mathtext", fontset="stixsans")

    # axes ticks and limits
    ax = plt.subplot(111)
    ax.set_xlim([min(masses), max(masses)])
    ax.tick_params(axis="y", which="both", direction="in", length=5, width=1, labelleft="off")
    ax.tick_params(axis="x", which="both", direction="in", length=5, width=1, labelbottom="off")

    # label the axes
    ax.set(xlabel="Molecular Mass (amu)")
    ax.set(ylabel="Probability Density Estimate")

    # Do the estimates and plot them
    xvals = np.arange(0, max(masses), 1)

    density_carbon = gkde(my_dict["Carbon Star"])
    density_carbon.covariance_factor = lambda: 0.5
    density_carbon._compute_covariance()

    density_dark = gkde(my_dict["Dark Cloud"])
    density_dark.covariance_factor = lambda: 0.5
    density_dark._compute_covariance()

    density_los = gkde(my_dict["LOS Cloud"])
    density_los.covariance_factor = lambda: 0.5
    density_los._compute_covariance()

    density_sfr = gkde(my_dict["SFR"])
    density_sfr.covariance_factor = lambda: 0.5
    density_sfr._compute_covariance()

    ax.plot(xvals[min(masses) : max(masses)], density_carbon(xvals[min(masses) : max(masses)]), color="darkorange")
    ax.fill_between(
        xvals[min(masses) : max(masses)],
        density_carbon(xvals[min(masses) : max(masses)]),
        0,
        facecolor="darkorange",
        alpha=0.25,
        zorder=4,
    )

    ax.plot(xvals[min(masses) : max(masses)], density_dark(xvals[min(masses) : max(masses)]), color="forestgreen")
    ax.fill_between(
        xvals[min(masses) : max(masses)],
        density_dark(xvals[min(masses) : max(masses)]),
        0,
        facecolor="forestgreen",
        alpha=0.25,
        zorder=4,
    )

    ax.plot(xvals[min(masses) : max(masses)], density_los(xvals[min(masses) : max(masses)]), color="red")
    ax.fill_between(
        xvals[min(masses) : max(masses)],
        density_los(xvals[min(masses) : max(masses)]),
        0,
        facecolor="red",
        alpha=0.25,
        zorder=4,
    )

    ax.plot(xvals[min(masses) : max(masses)], density_sfr(xvals[min(masses) : max(masses)]), color="dodgerblue")
    ax.fill_between(
        xvals[min(masses) : max(masses)],
        density_sfr(xvals[min(masses) : max(masses)]),
        0,
        facecolor="dodgerblue",
        alpha=0.25,
        zorder=4,
    )

    x_ann = 0.97
    y_ann = 0.96
    y_sep = 0.06

    ax.annotate(
        r"\underline{Source Types}",
        xy=(x_ann, y_ann - 0 * y_sep),
        xycoords="axes fraction",
        color="black",
        ha="right",
        va="top",
    )
    ax.annotate("SFR", xy=(x_ann, y_ann - y_sep), xycoords="axes fraction", color="dodgerblue", ha="right", va="top")
    ax.annotate(
        "Carbon Star", xy=(x_ann, y_ann - 2 * y_sep), xycoords="axes fraction", color="darkorange", ha="right", va="top"
    )
    ax.annotate(
        "Dark Cloud", xy=(x_ann, y_ann - 3 * y_sep), xycoords="axes fraction", color="forestgreen", ha="right", va="top"
    )
    ax.annotate("LOS Cloud", xy=(x_ann, y_ann - 4 * y_sep), xycoords="axes fraction", color="red", ha="right", va="top")

    plt.show()

    plt.savefig(
        filename if filename is not None else "mass_by_source_type_kde.pdf",
        format="pdf",
        transparent=True,
        bbox_inches="tight",
        pad_inches=0,
    )

    return


def waves_by_source_type(mol_list=None, filename=None):

    """
    Generates four pie charts, one for each generalized source type, with the wedges for the wavelengths used for first detections in those sources
    Defaults to using all molecules in the database, but can be passed a subset of molecules as a list of Molecule objects in mol_list.
    Default filename is 'waves_by_source_type.pdf', but that can also be overriden.
    """

    # If a list wasn't specified, default to all molecules
    if mol_list is None:
        mol_list = all_molecules

    # Close an old figure if it exists and initialize a new figure
    plt.close("Wavelength by Source Type")
    fig, axs = plt.subplots(2, 2, num="Wavelength by Source Type", figsize=(15, 12))
    plt.ion()

    # set some font defaults
    fontparams = {"size": 26, "family": "sans-serif", "sans-serif": ["Helvetica"]}
    plt.rc("font", **fontparams)
    plt.rc("mathtext", fontset="stixsans")

    # collect the data
    type_dict = {
        "Carbon Star": {"cm": 0, "mm": 0, "sub-mm": 0, "IR": 0, "UV": 0, "Vis": 0},
        "Dark Cloud": {"cm": 0, "mm": 0, "sub-mm": 0, "IR": 0, "UV": 0, "Vis": 0},
        "LOS Cloud": {"cm": 0, "mm": 0, "sub-mm": 0, "IR": 0, "UV": 0, "Vis": 0},
        "SFR": {"cm": 0, "mm": 0, "sub-mm": 0, "IR": 0, "UV": 0, "Vis": 0},
    }

    for mol in mol_list:
        # we'll make a dictionary here to flag if we've credited things already
        credit_dict = {
            "Carbon Star": False,
            "Dark Cloud": False,
            "LOS Cloud": False,
            "SFR": False,
        }

        for source in mol.sources:
            if source.type in type_dict:
                if credit_dict[source.type] is False:
                    for wave in mol.wavelengths:
                        type_dict[source.type][wave] += 1
                    credit_dict[source.type] = True

    # make the pie charts
    # Carbon Stars
    carbon_data = [
        type_dict["Carbon Star"]["cm"],
        type_dict["Carbon Star"]["mm"],
        type_dict["Carbon Star"]["sub-mm"],
        type_dict["Carbon Star"]["IR"],
        # type_dict['Carbon Star']['UV'] + type_dict['Carbon Star']['Vis'],
    ]

    carbon_colors = [
        "dodgerblue",
        "darkorange",
        "forestgreen",
        "black",
        #'violet',
    ]

    # Dark Clouds
    dark_data = [
        type_dict["Dark Cloud"]["cm"],
        type_dict["Dark Cloud"]["mm"],
        # type_dict['Dark Cloud']['sub-mm'],
        # type_dict['Dark Cloud']['IR'],
        # type_dict['Dark Cloud']['UV'] + type_dict['Dark Cloud']['Vis'],
    ]

    dark_colors = [
        "dodgerblue",
        "darkorange",
        #'forestgreen',
        #'black',
        #'violet',
    ]

    # LOS Clouds

    los_data = [
        type_dict["LOS Cloud"]["cm"],
        type_dict["LOS Cloud"]["mm"],
        type_dict["LOS Cloud"]["sub-mm"],
        type_dict["LOS Cloud"]["IR"],
        type_dict["LOS Cloud"]["UV"] + type_dict["LOS Cloud"]["Vis"],
    ]

    los_colors = [
        "dodgerblue",
        "darkorange",
        "forestgreen",
        "black",
        "violet",
    ]

    # SFRs

    sfr_data = [
        type_dict["SFR"]["cm"],
        type_dict["SFR"]["mm"],
        type_dict["SFR"]["sub-mm"],
        # type_dict['SFR']['IR'],
        # type_dict['SFR']['UV'] + type_dict['SFR']['Vis'],
    ]

    sfr_colors = [
        "dodgerblue",
        "darkorange",
        "forestgreen",
        #'black',
        #'violet',
    ]

    types = [
        "cm",
        "mm",
        "sub-mm",
        "IR",
        "UV/Vis",
    ]

    def make_labels(list):
        new_labels = []
        total = sum(list)
        for i in list:
            percent = 100 * i / total
            new_labels.append(r"{:.1f}\%".format(percent))
        return new_labels

    axs[0, 0].pie(
        carbon_data,
        labels=make_labels(carbon_data),
        colors=carbon_colors,
        labeldistance=1.1,
        wedgeprops={"linewidth": 1.0, "edgecolor": "black", "alpha": 0.5},
    )
    axs[0, 0].annotate(
        r"\textbf{Carbon Stars}",
        xy=(0.5, 1.05),
        xycoords="axes fraction",
        color="black",
        ha="center",
        va="top",
        size=30,
    )

    axs[0, 1].pie(
        dark_data,
        labels=make_labels(dark_data),
        colors=dark_colors,
        labeldistance=1.1,
        wedgeprops={"linewidth": 1.0, "edgecolor": "black", "alpha": 0.5},
    )
    axs[0, 1].annotate(
        r"\textbf{Dark Clouds}", xy=(0.5, 1.05), xycoords="axes fraction", color="black", ha="center", va="top", size=30
    )

    wedges, _ = axs[1, 0].pie(
        los_data,
        labels=make_labels(los_data),
        colors=los_colors,
        labeldistance=1.1,
        wedgeprops={"linewidth": 1.0, "edgecolor": "black", "alpha": 0.5},
    )
    axs[1, 0].annotate(
        r"\textbf{LOS Clouds}", xy=(0.5, 1.05), xycoords="axes fraction", color="black", ha="center", va="top", size=30
    )

    axs[1, 1].pie(
        sfr_data,
        labels=make_labels(sfr_data),
        colors=sfr_colors,
        labeldistance=1.1,
        wedgeprops={"linewidth": 1.0, "edgecolor": "black", "alpha": 0.5},
    )
    axs[1, 1].annotate(
        r"\textbf{SFRs}", xy=(0.5, 1.05), xycoords="axes fraction", color="black", ha="center", va="top", size=30
    )

    axs[0, 1].legend(wedges, types, title="Wavelengths", loc="center left", bbox_to_anchor=(1, 0, 0.5, 1))

    fig.tight_layout()
    fig.subplots_adjust(wspace=0.1, hspace=0.1)
    plt.show()

    plt.savefig(
        filename if filename is not None else "waves_by_source_type.pdf",
        format="pdf",
        transparent=True,
        bbox_inches="tight",
    )

    return
