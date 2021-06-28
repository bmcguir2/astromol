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
 
## Bugs, Errors, and Feature Requests
 
Should you encounter a bug, discover a factual error, or have a request for a new feature, please submit an issue using the GitHub issues tracker and providing as much description as possible.  

## Pull Requests

Direct contributions to the code are not being accepted at this time.

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

## Overall Structure

This package was written primarily to faciliate the McGuire 2018 living census paper.  As such, the figure- and table-making functions contained in it are tailor-made to that purpose, with limited flexibility or functionality for alteration.  However, to make the data as widely useful and accessible as possible, the package has been written to be object-oriented such that it is relatively straigthforward for anyone else to manipulate the dataset to suit individual needs.

### Object Classes

Data in `astromol` is stored in one of three classes: `Molecule`, `Telescope`, or `Source`.  Each molecule, telescope, or astronomical source entry in the database is a variable that comes preloaded with the `astromol` package.  Accessing these information therefore requires knowing these variable names.  These can be found either by inspecting the corresponding files (_molecules.py_, _telescopes.py_, and _sources.py_) or by using the helper function:

```Python
print_variables()
```

This function takes two optional arguments: `type = ` and `natoms = `.  The former can be set to `molecules`, `telescopes`, or `sources`, and the later can be used by specifying an integer and only molecules with that number of atoms will be printed.  If nothing is set, all variables will be printed.

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




## Planned Development

A number of upgrades are planned to the code, to the structure of the database, and to the content of the database.  These include the following:

* Inclusion of [SMILES](https://en.wikipedia.org/wiki/Simplified_molecular-input_line-entry_system) strings for each molecule.  There is a `smiles` attribute for the `Molecule` class for this purpose, but the data haven't been aggregated and input yet.
* Inclusion of links to relevant spectroscopic data -- largely CDMS and JPL rotational catalogs -- for each molecule.  
* Inclusion of outputs of quantum chemical structure calculations for each molecule, using the methods recommended in [Lee & McCarthy 2020](https://pubs.acs.org/doi/abs/10.1021/acs.jpca.9b09982).
* Migration of references from strings to entirely bibtex citekeys.  Implementation has particually begun for protoplanetary disk, extragalactic, and exoplanetary atmosphere molecules.
* Move to dealing with isotopologues as `Molecule` objects themselves.  This has already been implemented for protoplanetary disk species partially; references are still contained as a string in the main molecule.  A future update will deprecate strings entirely in favor of bibtex entries.
* Incorporation of isotopologue data for molecules other than those found in protoplanetary disks.  Likely to be a gradual effort.
* Addition of tracking for cometary species as well.
* Addition of data on telescopes and detection sources for non-ISM/CSM species.
