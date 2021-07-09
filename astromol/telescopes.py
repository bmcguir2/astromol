"""
An exceptional resource for these is 'Observatories and Telescopes of Modern Times' 
by David Leverington.
"""


class Telescope(object):
    """
    A class used to represent a telescope facility

    ...

    Attributes
    ----------
    name : str
        The full name of the facility.  All dimensions (i.e. 100-m) should be
        separated from units by a dash, not a space.
    shortname : str
        A shorter name for the facility, e.g.: GBT 100-m
    astromol_name : str
        A string representation of the variable representing this telescope
        in astromol
    type : str
        The type of facility, chosen from:
            "Single Dish"
            "Space"
            "Interferometer"
            "Optical"
            "Airborne"
        Note that any space-based telescope overrides other descriptors, any 
        airborne facility overrides other descriptors, and optical includes
        anything not radio (and neither airborne nor space-based)
    wavelength : list
        The primary wavelengths of operation of the facility, given as a
        a list of the following strings:
            "cm"
            "mm"
            "sub-mm"
            "IR"
            "Vis"
            "UV"
    latitude : float
        The latitude of the facility in degrees
    longitude : float
        The signed longitude of the facility in degrees
    diameter : float, int
        The diamater of the facility in meters
    built : int
        The year in which the facility came online
    decommissioned : int or None
        The year in which the facility ceased operation.
        (Default value of None is for facilities still in operation)
    notes : str
        General notes

    Properties
    ----------
    ndetects(mol_list=None)
        Returns the number of molecules from 'mol_list' this telescope was 
        used to detect (default is 'all_molecules')

    mols(mol_list=None)
        Returns a list of molecules from 'mol_list' this telescope was 
        used to detect (default is 'all_molecules')

    Methods
    -------
    inspect()
        Prints out a list to the terminal all of the attributes for this telescope and their values.

    summary()
        Prints out a nicely formatted list of properties and detected molecules to the terminal.          
    """

    def __init__(
        self,
        name=None,
        shortname=None,
        astromol_name = None,
        type=None,
        wavelength=None,
        latitude=None,
        longitude=None,
        diameter=None,
        built=None,
        decommissioned=None,
        notes=None,
    ):

        """
        Parameters
        ----------
        name : str
            The full name of the facility.  All dimensions (i.e. 100-m) should be
            separated from units by a dash, not a space.
        shortname : str
            A shorter name for the facility, e.g.: GBT 100-m
        astromol_name : str
            A string representation of the variable representing this telescope
            in astromol
        type : str
            The type of facility, chosen from:
                "Single Dish"
                "Space"
                "Interferometer"
                "Optical"
                "Airborne"
            Note that any space-based telescope overrides other descriptors, any 
            airborne facility overrides other descriptors, and optical includes
            anything not radio (and neither airborne nor space-based)
        wavelength : list
            The primary wavelengths of operation of the facility, given as a
            a list of the following strings:
                "cm"
                "mm"
                "sub-mm"
                "IR"
                "Vis"
                "UV"
        latitude : float
            The latitude of the facility in degrees
        longitude : float
            The signed longitude of the facility in degrees
        diameter : float, int
            The diamater of the facility in meters
        built : int
            The year in which the facility came online
        decommissioned : int or None
            The year in which the facility ceased operation.
            (Default value of None is for facilities still in operation)
        notes : str
            General notes   
        """ 

        self.name = name
        self.shortname = shortname
        self.astromol_name = astromol_name
        self.type = type
        self.wavelength = wavelength
        self.latitude = latitude
        self.longitude = longitude
        self.diameter = diameter
        self.built = built
        self.decommissioned = decommissioned
        self.notes = notes

    @property
    def ndetects(self, mol_list=None):
        """
        Returns the number of molecules from 'mol_list' that this telescope 
        was used to detect

        Parameters
        ----------
        mol_list : list, optional
            A list of molecule objects to sort through (default is 'all_molecules')

        Returns
        -------
        len(my_mols) : int
            The number of molecules this telescope was used to detect
        """

        if mol_list is None:
            from astromol.molecules import all_molecules
            mol_list = all_molecules

        my_mols = []
        for mol in mol_list:
            if self in mol.telescopes:
                my_mols.append(mol)
        return len(my_mols)

    @property
    def mols(self, mol_list=None, formulas=False):
        """
        Returns a list of molecules from 'mol_list' that this telescope was
        used to detect

        Parameters
        ----------
        mol_list : list
            A list of molecule objects to sort through (default is 'all_molecules')

        Returns
        -------
        my_mols : list
            A list of molecule objects that this telescope was used to detect
        """

        if mol_list is None:
            from astromol.molecules import all_molecules
            mol_list = all_molecules

        my_mols = []
        for mol in mol_list:
            if self in mol.telescopes:
                my_mols.append(mol)
        if formulas is True:
            return [mol.formula for mol in my_mols]
        else:
            return my_mols

    def inspect(self):
        """
        Prints out a list to the terminal all of the attributes for this telescope and their values.

        Note that this will be a largely unformatted output, and is really only intended to be used
        when a detailed look at what information is in the database for a given species is desired.
        """        

        # Get the attributes and values as a dictionary
        attr_dict = self.__dict__
        # Print them out
        for attr in attr_dict:
            print(f"{attr:20}: {attr_dict[attr] if attr_dict[attr] is not None else ''}")

    def summary(self):
        """
        Prints out a nicely formatted list of properties and detected molecules to the terminal.
        """

        from astromol.molecules import all_molecules

        print('-'*len(self.name))
        print(f"{self.name}")
        print('-'*len(self.name) + "\n")

        print(f"{'Type':11}\t{self.type}")
        if self.diameter:
            print(f"{'Diameter':11}\t{self.diameter}")
        print(f"{'Wavelengths':11}\t{', '.join(self.wavelength)}")
        if self.latitude:
            print(f"{'Latitude':11}\t{self.latitude}")
        if self.longitude:
            print(f"{'Longitude':11}\t{self.longitude}")
        print(f"{'Built':11}\t{self.built}")
        if self.decommissioned:
            print(f"{'Decommissioned':11}\t{self.decommissioned}")

        mol_str = f"Molecules Detected With {self.name}"
        print("\n" + "-"*len(mol_str))
        print(mol_str)
        print("-"*len(mol_str))

        detects = ', '.join([x.formula for x in all_molecules if self in x.telescopes])
        print(detects)        

GBT = Telescope(
    name="Green Bank Telescope",
    shortname="GBT 100-m",
    astromol_name="GBT",
    type="Single Dish",
    wavelength=["cm", "mm"],
    latitude=38.433056,
    longitude=-79.839722,
    diameter=100,
    built=2004,
)

IRAM30 = Telescope(
    name="IRAM 30-m",
    shortname="IRAM 30-m",
    astromol_name="IRAM30",
    type="Single Dish",
    wavelength=["mm", "sub-mm"],
    latitude=37.066161,
    longitude=-3.392719,
    diameter=30,
    built=1984,
)

Spitzer = Telescope(
    name="Spizter",
    shortname="Spitzer",
    astromol_name="Spitzer",
    type="Space",
    wavelength=["IR"],
    diameter=0.85,
    built=2003,
    decommissioned=2020,
)

Hubble = Telescope(
    name="Hubble Space Telescope",
    shortname="Hubble",
    astromol_name="Hubble",
    type="Space",
    wavelength=["IR", "Vis", "UV"],
    diameter=2.4,
    built=1990,
)

ALMA = Telescope(
    name="Atacama Large Millimeter/sub-millimeter Array",
    shortname="ALMA",
    astromol_name="ALMA",
    type="Interferometer",
    wavelength=["mm", "sub-mm"],
    latitude=-23.0193,
    longitude=-67.7532,
    built=2011,
)

ISO = Telescope(
    name="Infrared Space Observatory",
    shortname="ISO",
    astromol_name="ISO",
    type="Space",
    wavelength=["IR"],
    diameter=0.6,
    built=1995,
    decommissioned=1998,
)

NRAO140 = Telescope(
    name="NRAO 140-ft",
    shortname="NRAO 140-ft",
    astromol_name="NRAO140",
    type="Single Dish",
    wavelength=["cm"],
    latitude=38.433056,
    longitude=-79.839722,
    diameter=43,
    built=1965,
    decommissioned=2008,
    notes="Technically started operations again in 2014, but not for PI science.",
)

Algonquin46 = Telescope(
    name="Algonquin 46-m Telescope",
    shortname="Algonquin 46-m",
    astromol_name="Algonquin46",
    type="Single Dish",
    wavelength=["cm", "mm"],
    latitude=45.955503,
    longitude=-78.073042,
    diameter=46,
    built=1966,
    decommissioned=1987,
    notes="Technically in operation much longer, but seems to have ceased PI science in 1987.",
)

NRAOARO12 = Telescope(
    name="NRAO/ARO 12-m Telescope",
    shortname="NRAO/ARO 12-m",
    astromol_name="NRAOARO12",
    type="Single Dish",
    wavelength=["mm"],
    latitude=31.9533,
    longitude=-111.615,
    diameter=12,
    built=1984,
    notes="Originally the NRAO 36-ft (11 m) telescope until 1984, when the dish was replaced (these are counted as separate facilities).  The observatory was handed over to the ARO in 2000 and renamed.  In 2013, the antenna was replaced with a 12-m ALMA prototype antenna.",
)

NRAO36 = Telescope(
    name="NRAO 36-ft Telescope",
    shortname="NRAO 36-ft",
    astromol_name="NRAO36",
    type="Single Dish",
    wavelength=["mm"],
    latitude=31.9533,
    longitude=-111.615,
    diameter=11,
    built=1967,
    decommissioned=1984,
    notes="Became the NRAO/ARO 12-m telescope in 1984.",
)

Nobeyama45 = Telescope(
    name="Nobeyama 45-m Telescope",
    shortname="Nobeyama 45-m",
    astromol_name="Nobeyama45",
    type="Single Dish",
    wavelength=["cm", "mm"],
    latitude=35.9417,
    longitude=138.4758,
    diameter=45,
    built=1982,
)

Effelsberg100 = Telescope(
    name="Effelsberg 100-m Telescope",
    shortname="Effelsberg 100-m",
    astromol_name="Effelsberg100",
    type="Single Dish",
    wavelength=["cm"],
    latitude=50.5247,
    longitude=-6.8828,
    diameter=100,
    built=1972,
)

Haystack37 = Telescope(
    name="Haystack 37-m Telescope",
    shortname="Haystack",
    astromol_name="Haystack37",
    type="Single Dish",
    wavelength=["cm", "mm"],
    latitude=42.6233,
    longitude=-71.4882,
    diameter=37,
    built=1964,
)

PdBI = Telescope(
    name="Plateu de Bure Interferometer",
    shortname="PdBI",
    astromol_name="PdBI",
    type="Interferometer",
    wavelength=["mm"],
    latitude=44.63389,
    longitude=5.90792,
    built=1988,
    decommissioned=2016,
)

# NOEMA = Telescope(name='Northern Extended Millimeter Array',
#                     shortname='NOEMA',
#                     type='Interferometer',
#                     wavelength=['mm'],
#                     latitude=44.63389,
#                     longitude=5.90792,
#                     built=2016)

BIMA = Telescope(
    name="Berkeley-Illinois-Maryland Array",
    shortname="BIMA",
    astromol_name="BIMA",
    type="Interferometer",
    wavelength=["mm"],
    latitude=40.8178,
    longitude=-121.473,
    built=1986,
    decommissioned=2005,
    notes="Became part of CARMA.",
)

OVRO = Telescope(
    name="Caltech Owens Valley Radio Observatory Millimeter Array",
    shortname="OVRO",
    astromol_name="OVRO",
    type="Interferometer",
    wavelength=["mm"],
    latitude=37.2339,
    longitude=-118.282,
    built=1984,
    decommissioned=2005,
    notes="Became part of CARMA.",
)

Yebes40 = Telescope(
    name="Yebes RT40-m Telescope",
    shortname="Yebes 40-m",
    astromol_name="Yebes40",
    type="Single Dish",
    wavelength=["cm", "mm"],
    latitude=40525208,
    longitude=-3.088725,
    built=2007,
)

NRL85 = Telescope(
    name="Maryland Point Observatory Naval Research Lab 85-foot Telescope",
    shortname="NRL 85-ft",
    astromol_name="NRL85",
    type="Single Dish",
    wavelength=["cm"],
    diameter=26,
    latitude=38.3741667,
    longitude=-77.230833,
    built=1965,
    decommissioned=1994,
    notes="Primarily used for VLBI for much of its later years.",
)

ATCA = Telescope(
    name="Australia Telescope Compact Array",
    shortname="ATCA",
    astromol_name="ATCA",
    type="Interferometer",
    wavelength=["cm"],
    latitude=-30.312778,
    longitude=149.550278,
    built=1988,
)

Parkes64 = Telescope(
    name="Parkes 64-m Telescope",
    shortname="Parkes",
    astromol_name="Parkes64",
    type="Single Dish",
    wavelength=["cm"],
    diameter=64,
    latitude=-32.99778,
    longitude=148.26292,
    built=1961,
)

SMT10 = Telescope(
    name="ARO 10-m Submillimeter Telescope",
    shortname="SMT",
    astromol_name="SMT10",
    type="Single Dish",
    wavelength=["mm", "sub-mm"],
    diameter=10,
    latitude=32.701658,
    longitude=-109.871391,
    built=1993,
)

SEST15 = Telescope(
    name="Swedish-ESO 15-m Submillimetre Telescope",
    shortname="SEST",
    astromol_name="SEST15",
    type="Single Dish",
    wavelength=["mm", "sub-mm"],
    diameter=15,
    latitude=-29.26,
    longitude=-70.73,
    built=1987,
    decommissioned=2003,
)

Goldstone70 = Telescope(
    name='Goldstone 72-m (DSS-14; "Mars")',
    shortname="Goldstone",
    astromol_name="Goldstone70",
    type="Single Dish",
    wavelength=["cm"],
    diameter=70,
    latitude=35.426667,
    longitude=-116.89,
    built=1966,
    notes="Originally a 64-m dish; become 70-m in 1988. Conceivably still PI Science Capable?",
)

Mitaka6 = Telescope(
    name="Tokyo Astronomical Observatory Mitaka 6-m",
    shortname="Mitaka 6-m",
    astromol_name="Mitaka6",
    type="Single Dish",
    wavelength=["mm"],
    diameter=6,
    latitude=35.675217,
    longitude=139.538083,
    built=1970,
    decommissioned=2018,
    notes="Moved around quite a bit within Japan until returning (and retiring) in 2018.",
)

McMath = Telescope(
    name="McMath-Pierce Solar Telescope",
    shortname="McMath Solar Telescope",
    astromol_name="McMath",
    type="Optical",
    wavelength=["IR", "Vis", "UV"],
    diameter=1.6,
    latitude=31.9584,
    longitude=-111.595,
    built=1962,
)

Bell7m = Telescope(
    name="AT&T Bell Laboratories 7-m Telescope",
    shortname="Bell 7-m",
    astromol_name="Bell7m",
    type="Single Dish",
    wavelength=["cm"],
    diameter=7,
    built=1976,
    decommissioned=1992,
)

IRTF = Telescope(
    name="NASA Infrared Telescope Facility",
    shortname="IRTF",
    astromol_name="IRTF",
    type="Optical",
    wavelength=["IR"],
    diameter=3,
    latitude=19.8263,
    longitude=-155.473,
    built=1974,
)

KPNO4m = Telescope(
    name="Mayall 4-m Telescope",
    shortname="KPNO 4-m",
    astromol_name="KPNO4m",
    type="Optical",
    wavelength=["IR"],
    diameter=4,
    latitude=31.9583,
    longitude=-111.5967,
    built=1973,
)

Onsala20m = Telescope(
    name="Onsala 20-m Telescope",
    shortname="Onsala 20-m",
    astromol_name="Onsala20m",
    type="Single Dish",
    wavelength=["cm", "mm"],
    diameter=20,
    latitude=57.393056,
    longitude=11.917778,
    built=1976,
)

FCRAO14m = Telescope(
    name="Five College Radio Observatory 14-m Telescope",
    shortname="FCRAO 14-m",
    astromol_name="FCRAO14m",
    type="Single Dish",
    wavelength=["cm", "mm"],
    latitude=42.391925,
    longitude=-72.344097,
    built=1976,
    decommissioned=2005,
)

APEX = Telescope(
    name="Atacama Pathfinder Experiment",
    shortname="APEX",
    astromol_name="APEX",
    type="Single Dish",
    wavelength=["mm", "sub-mm"],
    diameter=12,
    latitude=-23.0058,
    longitude=-67.7592,
    built=2005,
)

CSO = Telescope(
    name="Caltech Submillimeter Observatory",
    shortname="CSO",
    astromol_name="CSO",
    type="Single Dish",
    wavelength=["mm", "sub-mm"],
    diameter=10.4,
    latitude=19.8225,
    longitude=-155.70694,
    built=1986,
    decommissioned=2015,
)

MWO4m = Telescope(
    name="University of Texas Millimeter Wave Observatory 4.9-m Telescope",
    shortname="MWO 4.9-m",
    astromol_name="MWO4m",
    type="Single Dish",
    wavelength=["mm"],
    latitude=30.3866,
    longitude=-97.7269,
    built=1971,
    decommissioned=1988,
)

HatCreek = Telescope(
    name="Hat Creek Station 20-ft Telescope",
    shortname="Hat Creek 20-ft",
    astromol_name="HatCreek",
    type="Single Dish",
    wavelength=["cm", "mm"],
    latitude=40.8178,
    longitude=-121.473,
    built=1965,
    decommissioned=1983,
    notes='Best build date found was "mid 1960s", so 1965 is an estimate.  It appears to have either been subsummed into BIMA or decommissioned when BIMA came online.  The decommissioning date is an estimate.',
)

SMA = Telescope(
    name="Submillimeter Array",
    shortname="SMA",
    astromol_name="SMA",
    type="Interferometer",
    wavelength=["mm"],
    latitude=19.8225,
    longitude=-155.70694,
    built=2003,
)

Herschel = Telescope(
    name="Herschel Space Telescope",
    shortname="Herschel",
    astromol_name="Herschel",
    type="Space",
    wavelength=["sub-mm", "IR"],
    built=2009,
    decommissioned=2013,
)

UKIRT = Telescope(
    name="United Kingdom Infrared Telescope",
    shortname="UKIRT",
    astromol_name="UKIRT",
    type="Optical",
    wavelength=["IR"],
    diameter=3.8,
    latitude=19.8225,
    longitude=-155.70694,
    built=1979,
)

SOFIA = Telescope(
    name="Stratospheric Observatory for Infrared Astronomy",
    shortname="SOFIA",
    astromol_name="SOFIA",
    type="Airborne",
    wavelength=["sub-mm", "IR"],
    diameter=2.5,
    built=2010,
)

Odin = Telescope(
    name="Odin",
    shortname="Odin",
    astromol_name="Odin",
    type="Space",
    wavelength=["sub-mm"],
    diameter=1.1,
    built=2001,
)

FUSE = Telescope(
    name="Far Ultraviolet Spectroscopic Explorer",
    shortname="FUSE",
    astromol_name="FUSE",
    type="Space",
    wavelength=["UV"],
    built=1999,
    decommissioned=2007,
)

Kuiper = Telescope(
    name="Kuiper Airborne Observatory",
    shortname="KAO",
    astromol_name="Kuiper",
    type="Airborne",
    wavelength=["sub-mm", "IR"],
    built=1974,
    decommissioned=1995,
)

MtHopkins = Telescope(
    name="Tillinghast 60 inch",
    shortname="Mt. Hopkins 60-in",
    astromol_name="MtHopkins",
    type="Optical",
    diameter=1.5,
    wavelength=["IR"],
    built=1969,
    latitude=31.6811,
    longitude=-110.878,
)

Aerobee = Telescope(
    name="Aerobee-150 Rocket",
    shortname="Aerobee-150 Rocket",
    astromol_name="Aerobee",
    type="Airborne",
    wavelength=["UV"],
    built=1970,
    decommissioned=1970,
    notes="They literally put a spectrometer on a rocket and shot it into the sky, then used the same spectrometer to measure H2 in the laboratory.",
)

Millstone = Telescope(
    name="Lincoln Laboratory Millstone Hill Observatory 84-ft",
    shortname="Millstone Hill 84-ft",
    astromol_name="Millstone",
    type="Single Dish",
    wavelength=["cm"],
    built=1956,
    decommissioned=1978,
    diameter=26,
    latitude=42.6233,
    longitude=-71.4882,
    notes="Originally the Ballistic Missile Early Warning System radar antenna.  Decommissioning date is a best guess based on the installation of larger telescopes to the site at that time.",
)

MtWilson = Telescope(
    name="Mount Wilson 100-in",
    shortname="Mt. Wilson",
    astromol_name="MtWilson",
    type="Optical",
    wavelength=["UV", "Vis"],
    diameter=2.54,
    built=1917,
    decommissioned=1989,
    notes="Now known as the Hooker Telescope.",
)

IRAS = Telescope(
    name="Infrared Astronomical Satellite",
    shortname="IRAS",
    astromol_name="IRAS",
    type="Space",
    wavelength=["IR"],
    diameter=0.60,
    built=1983,
    decommissioned=1983,
)


all_telescopes = [
    GBT,
    IRAM30,
    Spitzer,
    Hubble,
    ALMA,
    ISO,
    NRAO140,
    Algonquin46,
    NRAOARO12,
    NRAO36,
    Nobeyama45,
    Effelsberg100,
    Haystack37,
    PdBI,
    # NOEMA,
    BIMA,
    OVRO,
    Yebes40,
    NRL85,
    ATCA,
    Parkes64,
    SMT10,
    SEST15,
    Goldstone70,
    Mitaka6,
    McMath,
    Bell7m,
    IRTF,
    KPNO4m,
    Onsala20m,
    FCRAO14m,
    APEX,
    CSO,
    MWO4m,
    HatCreek,
    SMA,
    Herschel,
    UKIRT,
    SOFIA,
    Odin,
    FUSE,
    Kuiper,
    MtHopkins,
    Aerobee,
    Millstone,
    MtWilson,
    IRAS,
]
