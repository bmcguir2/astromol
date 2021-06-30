class Source(object):
    """
    A class used to represent an astronomical source

    ...

    Attributes
    ----------
    name : str
        The common name of the sources
    astromol_name : str
        A string representation of the variable representing this source
        in astromol        
    type : str
        The generalized source type, chosen from the following 
        (descriptions in [], do not include in string):
            LOS Cloud
            Dark Cloud
            Carbon Star
            SNR [Supernova Remnant]
            SFR [Star-forming Region]
            Shock
            PDR [Photodissociation Region]
            Protostar
            HII
            YSO [Young Stellar Object]
            PN [Planetary Nebula]
            Sgr A
            Oxygen Star
    ra : str
        The right ascension coordinates of the source in the format "hh:mm:ss.sss".
        hh and mm should be two digitse each, ss should be at least two, and may 
        include as many decimals as desired
    dec : str
        The declination coordinates of the source in the format "deg:min:sec" deg 
        and min should be two digits each deg can be signed, sec should be at least 
        two digits and may include as many decimals as desired
    simbad_url : str
        Url of the simbad entry for this source.  Should be truncated after 
        the Ident= string.

    Properties
    ----------
    ndetects(mol_list=None)
        Returns the number of molecules from 'mol_list' detected in the source
        (default is 'all_molecules')

    mols(mol_list=None)
        Returns a list of molecules from 'mol_list' detected in the source
        (default is 'all_molecules')

    ism_csm
        Is there a detection of an ISM/CSM molecule in this source? 

    Methods
    -------
    inspect()
        Prints out a list to the terminal all of the attributes for this source and their values.

    summary()
        Prints out a nicely formatted list of properties and detected molecules to the terminal.    
    """

    def __init__(
        self,
        name=None,
        astromol_name=None,
        type=None,
        ra=None,
        dec=None,
        simbad_url=None,
    ):
        """
        Parameters
        ----------
        name : str
            The common name of the sources
        astromol_name : str
            A string representation of the variable representing this source
            in astromol             
        type : str
            The generalized source type, chosen from the following 
            (descriptions in [], do not include in string):
                "LOS Cloud"
                "Dark Cloud"
                "Carbon Star"
                "SNR" [Supernova Remnant]
                "SFR" [Star-forming Region]
                "Shock"
                "PDR" [Photodissociation Region]
                "Protostar"
                "HII"
                "YSO" [Young Stellar Object]
                "PN" [Planetary Nebula]
                "Sgr A"
                "Oxygen Star"
        ra : str
            The right ascension coordinates of the source in the format "hh:mm:ss.sss".
            hh and mm should be two digitse each, ss should be at least two, and may 
            include as many decimals as desired
        dec : str
            The declination coordinates of the source in the format "deg:min:sec" deg 
            and min should be two digits each deg can be signed, sec should be at least 
            two digits and may include as many decimals as desired
        simbad_url : str
            Url of the simbad entry for this source.  Should be truncated after 
            the Ident= string.
        """

        self.name = name
        self.astromol_name = astromol_name
        self.type = type
        self.ra = ra
        self.dec = dec
        self.simbad_url = simbad_url

    @property
    def ndetects(self, mol_list=None):

        """
        Returns the number of molecules from 'mol_list' detected in the source

        Parameters
        ----------
        mol_list : list, optional
            A list of molecule objects to sort through (default is 'all_molecules')

        Returns
        -------
        len(my_mols) : int
            The number of molecules which have this source in their list of 
            detected locations
        """

        if mol_list is None:
            from astromol.molecules import all_molecules
            mol_list = all_molecules

        my_mols = []
        for mol in mol_list:
            if self in mol.sources:
                my_mols.append(mol)
        return len(my_mols)

    @property
    def mols(self, mol_list=None):
        """
        Returns a list of molecules from 'mol_list' detected in the source

        Parameters
        ----------
        mol_list : list
            A list of molecule objects to sort through (default is 'all_molecules')

        Returns
        -------
        my_mols : list
            A list of molecule objects which have this source in their list of 
            detected locations
        """
 
        if mol_list is None:
            from astromol.molecules import all_molecules
            mol_list = all_molecules

        my_mols = []
        for mol in mol_list:
            if self in mol.sources:
                my_mols.append(mol)
        return my_mols

    @property
    def ism_csm(self):
        """
        Is there a detection of an ISM/CSM molecule in this source?

        Returns
        -------
        bool
        """ 

        from astromol.molecules import all_molecules
        return any([self in x.sources for x in all_molecules])          

    def inspect(self):
        """
        Prints out a list to the terminal all of the attributes for this source and their values.

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
        print(f"{'RA (J2000)':11}\t{self.ra}")
        print(f"{'DEC (J2000)':11}\t{self.dec}")
        print(f"{'Simbad URL':11}\t{self.simbad_url}")

        mol_str = f"Molecules Detected in {self.name}"
        print("\n" + "-"*len(mol_str))
        print(mol_str)
        print("-"*len(mol_str))

        detects = ', '.join([x.formula for x in all_molecules if self in x.sources])
        print(detects)


"""
Source coordinates are generalized for simplicity, and individual detections may have 
been made (read: absolutely have been made) toward various pointing positions within 
these sources.  

Similarly, source types are highly generalized, to allow for some aggregate analysis.  

Finally, some detections have been made along the line of sight to these sources.  In 
these cases,  the source is catagorized as a 'LOS Cloud', regardless of the actual source
type, as it was only used as an absorbing background.

Diffuse, translucent, and dense clouds are all classified this way for simplicity.

Please check the individual papers for details.
"""

AFGL890LOS = Source(
    name="AFGL 890 LOS",
    astromol_name="AFGL890LOS",
    type="LOS Cloud",
    ra="06:10:48.0",
    dec="-06:12:00",
    simbad_url="http://simbad.u-strasbg.fr/simbad/sim-id?Ident=AFGL+890",
)

AFGL961LOS = Source(
    name="AFGL 961 LOS",
    astromol_name="AFGL961LOS",
    type="LOS Cloud",
    ra="06:34:37.741",
    dec="+04:12:44.20",
    simbad_url="http://simbad.u-strasbg.fr/simbad/sim-id?Ident=AFGL961",
)

AFGL989LOS = Source(
    name="AFGL 989 LOS",
    astromol_name="AFGL989LOS",
    type="LOS Cloud",
    ra="06:41:10.06",
    dec="+09:29:35.8",
    simbad_url="http://simbad.u-strasbg.fr/simbad/sim-id?Ident=AFGL989",
)

B1b = Source(
    name="B1-b",
    astromol_name="B1b",
    type="Dark Cloud",
    ra="03:33:20.8",
    dec="+31:07:34",
    simbad_url="http://simbad.u-strasbg.fr/simbad/sim-id?Ident=%5BHKM99%5D+B1-b",
)

CRL2688 = Source(
    name="CRL 2688",
    astromol_name="CRL2688",
    type="Carbon Star",
    ra="21:02:18.27",
    dec="+36:41:37.0",
    simbad_url="http://simbad.u-strasbg.fr/simbad/sim-id?Ident=CRL+2688",
)

CRL618 = Source(
    name="CRL 618",
    astromol_name="CRL618",
    type="Carbon Star",
    ra="04:42:53.62",
    dec="+36:06:53.40",
    simbad_url="http://simbad.u-strasbg.fr/simbad/sim-id?Ident=CRL+618",
)

CasALOS = Source(
    name="Cas A LOS",
    astromol_name="CasALOS",
    type="LOS Cloud",
    ra="23:23:24.00",
    dec="+58:48:54.0",
    simbad_url="http://simbad.u-strasbg.fr/simbad/sim-id?Ident=Cas+A",
)

CrabNebula = Source(
    name="Crab Nebula",
    astromol_name="CrabNebula",
    type="SNR",
    ra="05:34:31.94",
    dec="+22:00:52.2",
    simbad_url="http://simbad.u-strasbg.fr/simbad/sim-id?Ident=Crab+Nebula",
)

CygnusOB212LOS = Source(
    name="Cygnus OB2 - 12",
    astromol_name="CygnusOB212LOS",
    type="LOS Cloud",
    ra="20:32:40.96",
    dec="+41:04:13.2",
    simbad_url="http://simbad.u-strasbg.fr/simbad/sim-id?Ident=%4011680932&Name=Schulte%2012",
)

DR21 = Source(
    name="DR 21",
    astromol_name="DR21",
    type="SFR",
    ra="20:39:01.6",
    dec="+42:19:38",
    simbad_url="http://simbad.u-strasbg.fr/simbad/sim-basic?Ident=DR+21",
)

DR21LOS = Source(
    name="DR 21 LOS",
    astromol_name="DR21LOS",
    type="LOS Cloud",
    ra="20:39:01.6",
    dec="+42:19:38",
    simbad_url="http://simbad.u-strasbg.fr/simbad/sim-basic?Ident=DR+21",
)

DR21OH = Source(
    name="DR 21(OH)",
    astromol_name="DR21OH",
    type="SFR",
    ra="20:39:01.01",
    dec="+42:22:50.22",
    simbad_url="http://simbad.u-strasbg.fr/simbad/sim-basic?Ident=DR+21%28OH%29",
)

G0693 = Source(
    name="G+0.693-0.027",
    astromol_name="G0693",
    type="Shock",
    ra="17:47:21.9",
    dec="-28:21:27.00",
    simbad_url="http://simbad.u-strasbg.fr/simbad/sim-id?Ident=%404986185&Name=GCM%20%2b0.693%20-0.027",
)

G327306LOS = Source(
    name="G327.3-0.6 LOS",
    astromol_name="G327306LOS",
    type="LOS Cloud",
    ra="15:53:05.0",
    dec="-54:35:24",
    simbad_url="http://simbad.u-strasbg.fr/simbad/sim-id?Ident=G327.3-0.6",
)

GL2136LOS = Source(
    name="GL2136 LOS",
    astromol_name="GL2136LOS",
    type="LOS Cloud",
    ra="18:27:18.43",
    dec="-25:04:02.84",
    simbad_url="http://simbad.u-strasbg.fr/simbad/sim-id?Ident=%402520798&Name=GJ%20%202136%20B",
)

GalacticCenter = Source(
    name="Galactic Center", 
    astromol_name="GalacticCenter",
    type="SFR",
)

HD124314LOS = Source(
    name="HD 124314 LOS",
    astromol_name="HD124314LOS",
    type="LOS Cloud",
    ra="14:15:01.61",
    dec="-61:42:24.38",
    simbad_url="http://simbad.u-strasbg.fr/simbad/sim-id?Ident=HD+124314",
)

HD27778LOS = Source(
    name="HD 27778 LOS",
    astromol_name="HD27778LOS",
    type="LOS Cloud",
    ra="04:23:59.78",
    dec="24:18:03.53",
    simbad_url="http://simbad.u-strasbg.fr/simbad/sim-id?Ident=HD+27778",
)

HorseheadPDR = Source(
    name="Horsehead PDR",
    astromol_name="HorseheadPDR",
    type="PDR",
    ra="05:40:53.936",
    dec="-02:28:00",
    simbad_url="http://simbad.u-strasbg.fr/simbad/sim-id?Ident=%40828287&Name=NAME%20Horsehead%20Nebula",
)

IC443G = Source(
    name="IC 443G",
    astromol_name="IC443G",
    type="SNR",
    ra="06:16:43.4",
    dec="+22:32:24",
    simbad_url="http://simbad.u-strasbg.fr/simbad/sim-id?Ident=IC+443G",
)

IRAS16293 = Source(
    name="IRAS 16293",
    astromol_name="IRAS16293",
    type="Protostar",
    ra="16:32:22.56",
    dec="-24:28:31.8",
    simbad_url="http://simbad.u-strasbg.fr/simbad/sim-id?Ident=IRAS+16293-2422",
)

IRC10216 = Source(
    name="IRC+10216",
    astromol_name="IRC10216",
    type="Carbon Star",
    ra="09:47:57.406",
    dec="+13:16:43.56",
    simbad_url="http://simbad.u-strasbg.fr/simbad/sim-id?Ident=IRC%2B10216",
)

K350 = Source(
    name="K3-50",
    astromol_name="K350",
    type="HII",
    ra="20:04:45.59",
    dec="33:32:42.0",
    simbad_url="http://simbad.u-strasbg.fr/simbad/sim-id?Ident=%402904320&Name=NAME%20K%203-50A",
)

L134 = Source(
    name="L134",
    astromol_name="L134",
    type="Dark Cloud",
    ra="15:53:36.3",
    dec="-04:35:26.0",
    simbad_url="http://simbad.u-strasbg.fr/simbad/sim-id?Ident=%402622026&Name=LDN%20%20134",
)

L1527 = Source(
    name="L1527",
    astromol_name="L1527",
    type="Dark Cloud",
    ra="04:39:53.0",
    dec="+25:45:00",
    simbad_url="http://simbad.u-strasbg.fr/simbad/sim-id?Ident=L1527",
)

L1544 = Source(
    name="L1544",
    astromol_name="L1544",
    type="Dark Cloud",
    ra="05:04:16.6",
    dec="+25:10:48",
    simbad_url="http://simbad.u-strasbg.fr/simbad/sim-id?Ident=L1544",
)

L183 = Source(
    name="L183",
    astromol_name="L183",
    type="Dark Cloud",
    ra="15:54:12.2",
    dec="-02:49:42",
    simbad_url="http://simbad.u-strasbg.fr/simbad/sim-basic?Ident=L183",
)

L483 = Source(
    name="L483",
    astromol_name="L483",
    type="Dark Cloud",
    ra="18:17:35.0",
    dec="-04:39:48",
    simbad_url="http://simbad.u-strasbg.fr/simbad/sim-id?Ident=L483",
)

LOSCloud = Source(
    name="LOS Cloud",
    astromol_name="LOSCloud",
    type="LOS Cloud",
)

Lupus1A = Source(
    name="Lupus-1A",
    astromol_name="Lupus1A",
    type="Dark Cloud",
    ra="15:42:52.4",
    dec="-34:07:53.5",
    simbad_url="http://simbad.u-strasbg.fr/simbad/sim-id?Ident=%405549180&Name=NAME%20Lupus-1A",
)

M17LOS = Source(
    name="M17 LOS",
    astromol_name="M17LOS",
    type="LOS Cloud",
    ra="18:20:47",
    dec="-16:10:18",
    simbad_url="http://simbad.u-strasbg.fr/simbad/sim-id?Ident=M17",
)

M17SW = Source(
    name="M17SW",
    astromol_name="M17SW",
    type="PDR",
    ra="18:20:23.1",
    dec="-16:11:43",
    simbad_url="http://simbad.u-strasbg.fr/simbad/sim-id?Ident=M17+SW",
)

M3LOS = Source(
    name="M3 LOS",
    astromol_name="M3LOS",
    type="LOS Cloud",
    ra="13:42:11.62",
    dec="+28:22:38.2",
    simbad_url="http://simbad.u-strasbg.fr/simbad/sim-id?Ident=M3",
)

NGC2024 = Source(
    name="NGC 2024",
    astromol_name="NGC2024",
    type="PDR",
    ra="05:41:43",
    dec="-01:50:30",
    simbad_url="http://simbad.u-strasbg.fr/simbad/sim-id?Ident=NGC+2024",
)

NGC2024LOS = Source(
    name="NGC 2024 LOS",
    astromol_name="NGC2024LOS",
    type="LOS Cloud",
    ra="05:41:43",
    dec="-01:50:30",
    simbad_url="http://simbad.u-strasbg.fr/simbad/sim-id?Ident=NGC+2024",
)

NGC2264 = Source(
    name="NGC 2264",
    astromol_name="NGC2264",
    type="YSO",
    ra="06:40:58",
    dec="+09:53:42",
    simbad_url="http://simbad.u-strasbg.fr/simbad/sim-id?Ident=NGC+2264",
)

NGC6334 = Source(
    name="NGC 6334",
    astromol_name="NGC6334",
    type="SFR",
    ra="17:20:53.3",
    dec="-35:46:59",
    simbad_url="http://simbad.u-strasbg.fr/simbad/sim-id?Ident=%402361705&Name=NAME%20NGC%206334-I",
)

NGC6334LOS = Source(
    name="NGC 6334 LOS",
    astromol_name="NGC6334LOS",
    type="LOS Cloud",
    ra="17:20:53.3",
    dec="-35:46:59",
    simbad_url="http://simbad.u-strasbg.fr/simbad/sim-id?Ident=%402361705&Name=NAME%20NGC%206334-I",
)

NGC7023 = Source(
    name="NGC 7023",
    astromol_name="NGC7023",
    type="PDR",
    ra="21:01:36.9",
    dec="+68:09:48",
    simbad_url="http://simbad.u-strasbg.fr/simbad/sim-id?Ident=NGC+7023",
)

NGC7027 = Source(
    name="NGC 7027",
    astromol_name="NGC7027",
    type="PN",
    ra="21:07:01.8",
    dec="+42:14:10",
    simbad_url="http://simbad.u-strasbg.fr/simbad/sim-id?Ident=NGC+7023",
)

NGC7538 = Source(
    name="NGC 7538",
    astromol_name="NGC7538",
    type="YSO",
    ra="23:13:37.2",
    dec="61:30:00",
    simbad_url="http://simbad.u-strasbg.fr/simbad/sim-id?Ident=NGC+7538",
)

NGC7538LOS = Source(
    name="NGC 7538 LOS",
    astromol_name="NGC7538LOS",
    type="LOS Cloud",
    ra="23:13:37.2",
    dec="61:30:00",
    simbad_url="http://simbad.u-strasbg.fr/simbad/sim-id?Ident=NGC+7538",
)

Orion = Source(
    name="Orion",
    astromol_name="Orion",
    type="SFR",
    ra="05:35:14.16",
    dec="-05:22:21.5",
    simbad_url="http://simbad.u-strasbg.fr/simbad/sim-id?Ident=Orion+KL",
)

OrionBar = Source(
    name="Orion Bar",
    astromol_name="OrionBar",
    type="PDR",
    ra="05:35:22.30",
    dec="-05:24:33.0",
    simbad_url="http://simbad.u-strasbg.fr/simbad/sim-id?Ident=Orion+Bar",
)

rhoOphA = Source(
    name="rho Ophiuchi A",
    astromol_name="rhoOphA",
    type="SFR",
    ra="16:26:27.20",
    dec="-24:24:04",
    simbad_url="http://simbad.u-strasbg.fr/simbad/sim-id?Ident=%402488602&Name=NAME%20rho%20Oph%20A%20SM%201",
)

SgrA = Source(
    name="Sgr A",
    astromol_name="SgrA",
    type="Sgr A",
    ra="17:45:40.0",
    dec="-29:00:28.2",
    simbad_url="http://simbad.u-strasbg.fr/simbad/sim-id?Ident=Sgr+A",
)

SgrALOS = Source(
    name="Sgr A LOS",
    astromol_name="SgrALOS",
    type="LOS Cloud",
    ra="17:45:40.0",
    dec="-29:00:28.2",
    simbad_url="http://simbad.u-strasbg.fr/simbad/sim-id?Ident=Sgr+A",
)

SgrB2 = Source(
    name="Sgr B2",
    astromol_name="SgrB2",
    type="SFR",
    ra="17:47:20.4",
    dec="-28:23:07",
    simbad_url="http://simbad.u-strasbg.fr/simbad/sim-id?Ident=Sgr+B2",
)

SgrB2LOS = Source(
    name="Sgr B2 LOS",
    astromol_name="SgrB2LOS",
    type="LOS Cloud",
    ra="17:47:20.4",
    dec="-28:23:07",
    simbad_url="http://simbad.u-strasbg.fr/simbad/sim-id?Ident=Sgr+B2",
)

TC1 = Source(
    name="TC 1",
    astromol_name="TC1",
    type="PN",
    ra="17:45:35.29",
    dec="-46:05:23.7",
    simbad_url="http://simbad.u-strasbg.fr/simbad/sim-id?Ident=PN%20Tc%201%20",
)

TMC1 = Source(
    name="TMC-1",
    astromol_name="TMC1",
    type="Dark Cloud",
    ra="04:41:45.9",
    dec="+25:41:27",
    simbad_url="http://simbad.u-strasbg.fr/simbad/sim-id?Ident=TMC-1",
)

VYCaMaj = Source(
    name="VY Ca Maj",
    astromol_name="VYCaMaj",
    type="Oxygen Star",
    ra="07:22:58.3",
    dec="-25:46:03.2",
    simbad_url="http://simbad.u-strasbg.fr/simbad/sim-id?Ident=VY+Canis+Majoris",
)

W3 = Source(
    name="W3",
    astromol_name="W3",
    type="LOS Cloud",
    ra="02:27:04.10",
    dec="+61:52:27.1",
    simbad_url="http://simbad.u-strasbg.fr/simbad/sim-id?Ident=W3",
)

W3OH = Source(
    name="W3(OH)",
    astromol_name="W3OH",
    type="SFR",
    ra="02:27:04.1",
    dec="+61:52:52",
    simbad_url="http://simbad.u-strasbg.fr/simbad/sim-id?Ident=W3+%28OH%29&NbIdent=1",
)

W31LOS = Source(
    name="W31 LOS",
    astromol_name="W31LOS",
    type="LOS Cloud",
    ra="18:10:28.6",
    dec="-19:55:51",
    simbad_url="http://simbad.u-strasbg.fr/simbad/sim-id?Ident=W31",
)

W33LOS = Source(
    name="W33 LOS",
    astromol_name="W33LOS",
    type="LOS Cloud",
    ra="18:14:14.0",
    dec="-17:55:50",
    simbad_url="http://simbad.u-strasbg.fr/simbad/sim-id?Ident=W33",
)

W43LOS = Source(
    name="W43 LOS",
    astromol_name="W43LOS",
    type="LOS Cloud",
    ra="18:47:32.4",
    dec="-01:56:31",
    simbad_url="http://simbad.u-strasbg.fr/simbad/sim-id?Ident=W43",
)

W44LOS = Source(
    name="W44 LOS",
    astromol_name="W44LOS",
    type="LOS Cloud",
    ra="18:56:10.65",
    dec="+01:13:21.30",
    simbad_url="http://simbad.u-strasbg.fr/simbad/sim-id?Ident=W43",
)

W49 = Source(
    name="W49",
    astromol_name="W49",
    type="SFR",
    ra="19:10:19.6",
    dec="+09:07:42",
    simbad_url="http://simbad.u-strasbg.fr/simbad/sim-id?Ident=W49",
)

W49LOS = Source(
    name="W49 LOS",
    astromol_name="W49LOS",
    type="LOS Cloud",
    ra="19:10:19.6",
    dec="+09:07:42",
    simbad_url="http://simbad.u-strasbg.fr/simbad/sim-id?Ident=W49",
)

W51 = Source(
    name="W51",
    astromol_name="W51",
    type="SFR",
    ra="19:23:50",
    dec="+14:06:0",
    simbad_url="http://simbad.u-strasbg.fr/simbad/sim-id?Ident=W51",
)

W51LOS = Source(
    name="W51 LOS",
    astromol_name="W51LOS",
    type="LOS Cloud",
    ra="19:23:50",
    dec="+14:06:0",
    simbad_url="http://simbad.u-strasbg.fr/simbad/sim-id?Ident=W51",
)

XiPerLOS = Source(
    name="Xi Per LOS",
    astromol_name="XiPerLOS",
    type="LOS Cloud",
    ra="03:58:57.9",
    dec="+35:47:27.74",
    simbad_url="http://simbad.u-strasbg.fr/simbad/sim-id?Ident=Xi+Per",
)

rhoOphA = Source(
    name="rho Oph A",
    astromol_name="rhoOphA",
    type="SFR",
    ra="16:25:35.14",
    dec="-23:26:49.9",
    simbad_url="http://simbad.u-strasbg.fr/simbad/sim-id?Ident=Rho+Oph+A",
)

all_sources = [
    AFGL890LOS,
    AFGL961LOS,
    AFGL989LOS,
    B1b,
    CRL2688,
    CRL618,
    CasALOS,
    CrabNebula,
    CygnusOB212LOS,
    DR21,
    DR21LOS,
    DR21OH,
    G0693,
    G327306LOS,
    GL2136LOS,
    GalacticCenter,
    HD124314LOS,
    HD27778LOS,
    HorseheadPDR,
    IC443G,
    IRAS16293,
    IRC10216,
    K350,
    L134,
    L1527,
    L1544,
    L183,
    L483,
    LOSCloud,
    Lupus1A,
    M17LOS,
    M17SW,
    M3LOS,
    NGC2024,
    NGC2024LOS,
    NGC2264,
    NGC6334,
    NGC6334LOS,
    NGC7023,
    NGC7027,
    NGC7538,
    NGC7538LOS,
    Orion,
    OrionBar,
    rhoOphA,
    SgrA,
    SgrALOS,
    SgrB2,
    SgrB2LOS,
    TC1,
    TMC1,
    VYCaMaj,
    W3,
    W3OH,
    W31LOS,
    W33LOS,
    W43LOS,
    W44LOS,
    W49,
    W49LOS,
    W51,
    W51LOS,
    XiPerLOS,
    rhoOphA,
]
