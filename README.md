# astromol

## A Database of Molecules Detected in Space

`astromol` is a Python 3 package that provides a database of molecules detected in space and an object-oriented interface for interacting with the database and generating the figures and tables from the [_Census of Interstellar, Circumstellar, Extragalactic, Protoplanetary Disk, and Exoplanetary Molecules_ by B. McGuire](https://iopscience.iop.org/article/10.3847/1538-4365/aae5d2).

If you use `astromol` for your own work, please cite the Zenodo entry: [![DOI]()]().

## Setup Instructions

Clone the `astromol` repository to your computer.  It is **not** recommended to clone the `development` branch.

In your terminal, navigate to the `astromol` directory, then install using

```
pip install -e .
```

The use of the `-e` flag creates a symlink to the package which will enable you to easily check for and download updates to the code by simply using

```
git pull
```

in the `astromol` directory without needing to re-install.

## Version Control

The version number of `astromol` is given in the format XXXX.Y.Z, with each number corresponding to a different level of update:
* XXXX is a year (i.e. 2018 or 2021) corresponding to the most recent major update of the McGuire 2018 living census paper.
* Y is reset to 0 with each update of XXXX, and is incremented by 1 anytime a new molecule (or molecules) is added to the database between census releases.
* Z is reset to 0 with each update of XXXX or Y, and is incremented by 1 anytime an update is made to the code base that is not related to the addition of new molecules or the release of a new census.

The version number for `astromol` is accessible by:

```python
from astromol import version
version()
```

It is also stored in the `__version__` variable.

The date of the last update to `astromol` is available as well by:

```python
from astromol import updated
updated()
```

It is also stored in the `__updated__` variable as a Python `datetime` object.

Updates to documentation, such as this readme file, that do not accompany an actual change in the codebase are _not_ incremented in the version numbers, but are of course recorded as commits in the GitHub repository.

## Overall Structure and Usage

This package was written primarily to faciliate the McGuire 2018 living census paper.  As such, the figure- and table-making functions contained in it are tailor-made to that purpose, with limited flexibility or functionality for alteration.  However, to make the data as widely useful and accessible as possible, the package has been written to be object-oriented such that it is relatively straigthforward for anyone else to manipulate the dataset to suit individual needs.

The most useful way to use the package is to do the (dreaded) import \*:

```Python
from astromol import *
```

This will pre-load all of the information on detected species, facilities, and astronomical sources as variables into a Python session to work with.

### Object Classes

Data in `astromol` is stored in one of three classes: `Molecule`, `Telescope`, or `Source`.  Each molecule, telescope, or astronomical source entry in the database is a variable that comes preloaded with the `astromol` package.  Accessing these information therefore requires knowing these variable names.  These can be found either by inspecting the corresponding files (_molecules.py_, _telescopes.py_, and _sources.py_) or by using the helper function:

```Python
print_variables()
```

This function takes two optional arguments: `type = ` and `natoms = `.  The former can be set to `molecules`, `telescopes`, or `sources`, and the later can be used by specifying an integer and only molecules with that number of atoms will be printed.  If nothing is set, all variables will be printed.

Some helper lists have been pre-loaded, these are:

```Python
all_molecules
all_telescopes
all_sources
```

These are just lists containing all variables in the database of the corresponding type.  This isn't useful for direct inspection (as it will just return Object IDs), but is very useful for looping over.

#### The Molecule Class

This is the primary data container for the package.  Each molecular species is represented by a `Molecule` object.  These objects have several dozen attributes, a full list of which can be found by calling:

```Python
help(Molecule)
```

Many of these attributes are placeholders and are not yet filled with data.  To see a raw display of all the possible attributes for a molecule, and see which portions have data for that particular entry, use the Molecule Class method `inspect`.  For example, to see the data for methanol (CH3OH), use:

```Python
CH3OH.inspect()
```

Alternatively, the function ```inspect``` can be called, although this just calls the underlying class method anyway:

```Python
inspect(CH3OH)
```

A more nicely formatted summary of much of the pertinent data can be achieved using the Molecule Class method `summary`:

```Python
CH3OH.summary()
```

or, again, the function ```summary``` is available to call the underlying class method:

```Python
summary(CH3OH)
```

#### The Telescope Class

This data container holds information on the telescopes used to detect molecules.  These include currently:

* Name
* Shorter name / abbreviation
* Type of facility
* Generalized operational wavelength ranges
* Latitude and longitude in degrees
* Diameter of the dish (when appropriate)
* Dates of commissioning and decommissioning

This class has the same inspection command to view the entire contents:

```Python
ALMA.inspect()
```

As well, the overall ```inspect``` function can accept Telescope objects as an argument.

#### The Source Class

This data container holds information on the sources in which molecules are detected.  Right now these are only for ISM/CSM species.  The data gathered includes:

* Name
* Generalized type
* RA and Dec (hh:mm:ss and deg:min:sec)
* Direct link to the simbad entry for this source [may actually show up as a clickable hyperlink if view in Jupyter Notebooks]

This class has the same inspection command to view the entire contents:

```Python
SgrB2.inspect()
```

As well, the overall ```inspect``` function can accept Source objects as an argument.

### Functions

As discussed above, nearly all functions provided in the `functions.py` file are used for generating the figures, tables, and other minutae for the census document.  Also as mentioned earlier, some limited customizability has been built into these functions.  For all plotting functions, it is possible, for example, to specify a custom list of molecules on which to operate.

For example:

```Python
cumu_det_plot()
```

will generate the plot of the cumulative number of detections over time, using all molecules in the database.  However, one could choose to instead plot the cumulative number of ionic species detected in the ISM by:

```Python
ions = [x for x in all_molecules if (x.cation or x.anion)]
cumu_det_plot(mol_list = ions)
```

It's also possible to modify the filename that is used for output:

```Python
ions = [x for x in all_molecules if (x.cation or x.anion)]
cumu_det_plot(mol_list = ions, filename = 'cumulative_ions.pdf')
```

Custom plotting functions can of course be written to generate whatever plots are desired using the information in the database.  As well, one can always modify the built-in functions to change labels, colors, etc. to suit a need or preference.

### PowerPoint Slide of Detected ISM Molecules

The ```make_mols_slide()``` function will generate a slide containing a formatted dispaly of all detected ISM/CSM molecules to date, in widescreen PowerPoint format, sorted by the number of atoms.  It will display as well the total number of species, the date of the most recent update, the version of `astromol` used to generate the slide, and the appropriate literature reference.  This, too, can take a modified list of species with the ```mol_list = ``` optional argument.

Note that in the current implementation, the placement of the columns is done manually.  While it can adapt a bit to changes in the list, the adaptation is not perfect so some adjustment afterward might be needed.  Additionally, if the custom list doesn't contain molecules of a given number of atoms, PowerPoint will likely say the file is "broken" and ask to "repair" it.  Choosing the repair option is fine, and will produce the slide as normal.

Updates are planned for this feature down the road to make it more versatile.

## BibTeX Databases

Alongside the codebase itself, the installation comes with (currently) three \*.bib files in the `bibtex` folder:

```
exgal_refs.bib
exo_refs.bib
ppd_refs.bib
```

Again, the primary purpose of these is to interface with the functions to generate LaTeX tables for the census, but they may be useful.

These are also maintained as libraries in the NASA ADS system, accessible here:

* [Extragalactic Molecule Detections](https://ui.adsabs.harvard.edu/public-libraries/uEvz0RU1TvyPrxPqSgVzyw)
* [Exoplanet Atmosphere Molecule Detections](https://ui.adsabs.harvard.edu/public-libraries/m_TQ7PZPQMGNOKhAR1Zb9w)
* [Protoplanetary Disk Molecule Detections](https://ui.adsabs.harvard.edu/public-libraries/D8GfNtC7RLO3LoSICiQLbw)

## Planned Development

A number of upgrades are planned to the code, to the structure of the database, and to the content of the database.  If an item is on this list, please do not open a new issue on GitHub to request it; emails to [brettmc@mit.edu](mailto:brettmc@mit.edu) are welcome to express a prioritization desire. Current planned upgrades include the following:

* Inclusion of [SMILES](https://en.wikipedia.org/wiki/Simplified_molecular-input_line-entry_system) strings for each molecule.  There is a `smiles` attribute for the `Molecule` class for this purpose, but the data haven't been aggregated and input yet.
* Inclusion of links to relevant spectroscopic data -- largely CDMS and JPL rotational catalogs -- for each molecule.  
* Inclusion of outputs of quantum chemical structure calculations for each molecule, using the methods recommended in [Lee & McCarthy 2020](https://pubs.acs.org/doi/abs/10.1021/acs.jpca.9b09982).
* Migration of references from strings to entirely bibtex citekeys.  Implementation has particually begun for protoplanetary disk, extragalactic, and exoplanetary atmosphere molecules.
* Move to dealing with isotopologues as `Molecule` objects themselves.  This has already been implemented for protoplanetary disk species partially; references are still contained as a string in the main molecule.  A future update will deprecate strings entirely in favor of bibtex entries.
* Incorporation of isotopologue data for molecules other than those found in protoplanetary disks.  Likely to be a gradual effort.
* Addition of tracking for cometary species as well.
* Addition of data on telescopes and detection sources for non-ISM/CSM species.
* Coordinates for sources will eventually be updated to be `astropy` objects
* Links to telescope facility websites
* Increase versatility of the generation of the PowerPoint slide by allowing selective coloring or bolding of molecules given certain properties.
* A database of publication objects and/or authors is being considered, but is a long way off.

## Bugs, Errors, and Feature Requests
 
Should you encounter a bug, discover a factual error, or have a request for a new feature, please submit an issue using the GitHub issues tracker and providing as much description as possible.  

## Pull Requests

Direct contributions to the codebase are not being accepted at this time.
