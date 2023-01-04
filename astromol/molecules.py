from astromol.sources import *
from astromol.telescopes import *
from astromol.references import *
import bibtexparser as btp
import re
import numpy as np
#from rdkit import Chem
#from rdkit.Chem.Draw import IPythonConsole

#IPythonConsole.ipython_useSVG = True


class Molecule(object):
    """
    A class used to represent a molecule detected in space

    ...

    Attributes
    ----------
    name : str
        The full name of the molecule
    formula : str
        The chemical formula for the molecule.  Should be interpretable
        by the latex mhchem package.  Does not include isomer or structural
        identifiers such as 'l' or 'c' for linear or cyclic.
    smiles : str
        The smiles string for the molecule. Not currently used.
    structure_params : dict
        A dictionary containing default values to pass to the structure function.
        keys and types are:
            showcarbon : bool (default is False) 
    table_formula : str
        The same as formula, but can include isomer or structural identifiers
        for listing in the Tables of Molecules.
    label : str
        The label that will be used for hypertext linking in LaTeX.  Labels
        must be unique among the database.    
    astromol_name : str
        A string representation of the variable representing this molecule
        in astromol

    year : int
        The year the molecule was first detected in space

    sources : list
        A list of source objects in which the molecule is detected
    telescopes : list
        A list of telescope objects that were used to perform the detections
    wavelengths : list
        A list of wavelength strings describing the wavelengths at which
        the detection was made.  Choose only from among this list:
            "cm"
            "mm"
            "sub-mm"
            "IR"
            "Vis"
            "UV"

    cyclic : bool
        A flag used to indicate if the molecule contains a cyclic moiety
        within its structure (default is False)
    fullerene : bool
        A flag used to indicate if the molecule is a fullerene (default is False). 
        Fullerene and pah should not both be True. 
    pah : bool
        A flag used to indicate if the molecule is a polycyclic aromatic
        hydrocarbon (default is False).  Fullerene and pah should not both be True. 

    Acon : int
        The A rotational constant of the molecule in units of MHz
    Bcon : int
        The B rotational constant of the molecule in units of MHz
    Ccon : int
        The C rotational constant of the molecule in units of MHz

    mua : float
        The A-component of the dipole moment of the molecule in units of Debye
    mub : float
        The B-component of the dipole moment of the molecule in units of Debye
    muc : float
        The C-component of the dipole moment of the molecule in units of Debye

    d_refs : str
        The literature references which describe the ISM/CSM detection of the molecule.
        Multiple references should be separated by a semi-colon (';'), which should
        not otherwise be used in the string.
    l_refs : str
        The literature references which describe the laboratory work that enabled
        the ISM/CSM detection of the molecule.  Multiple references should be separated
        by a semi-colon (';'), which should not otherwise be used in the string.        
    d_ref_bib_ids : list
        A list of strings, with each string a citekey for a bibtex entry corresponding
        to literature references which describe the ISM/CSM detection of the molecule.
        Not currently used.
    l_ref_bib_ids : list
        A list of strings, with each string a citekey for a bibtex entry corresponding
        to literature references which describe the laboratory work that enabled
        the ISM/CSM detection of the molecule.  Not currently used.

    isotopologues : str
        A comma-separated list of isotopically substituted species that have been detected.
        This will be updated to be a list of molecule objects instead in a future version.
    isos_d_refs : str
        Should start with the isotopologue formula in [], followed by the literature references
        describing the ISM/CSM detection of the molecule.
    isos_l_refs : str
        Should start with the isotopologue formula in [], followed by the literature references
        describing the laboratory work that enabled the ISM/CSM detection of the molecule.

    ice : bool
        Flag to indicate if the molecule has been detected in an ice (default is False)
    ice_sources : list
        A list of source objects in which the molecule has been detected in an ice.
    ice_telescopes : list
        A list of telescope objects that were used to perform the detections of the
        molecule in an ice.
    ice_wavelengths : list
        A list of wavelength strings describing the wavelengths at which
        the detection in an ice was made.  Choose only from among this list:
            "cm"
            "mm"
            "sub-mm"
            "IR"
            "Vis"
            "UV"
    ice_d_refs : str
        The literature references which describe the detection of the molecule in an ice.
        Multiple references should be separated by a semi-colon (';'), which should
        not otherwise be used in the string.
    ice_l_refs : str
        The literature references which describe the laboratory work that enabled the 
        detection of the molecule in an ice.  Multiple references should be separated by 
        a semi-colon (';'), which should not otherwise be used in the string.
    ice_d_bib_ids : list
        A list of strings, with each string a citekey for a bibtex entry corresponding
        to literature references which describe the detection of the molecule in an ice.
        Not currently used.
    ice_l_bib_ids : list
        A list of strings, with each string a citekey for a bibtex entry corresponding
        to literature references which describe the laboratory work that enabled
        the detection of the molecule in an ice.  Not currently used.
    ice_isos : list
        A list of molecule objects, each one describing an isotopically substituted species
        detected in an ice.     

    ppd : bool
        Flag to indicate if the molecule has been detected in a Class II 
        protoplanetary disk (default is False)
    ppd_sources : list
        A list of source objects in which the molecule has been detected in a Class II 
        protoplanetary disk.
    ppd_telescopes : list
        A list of telescope objects that were used to perform the detections of the
        molecule in a Class II protoplanetary disk.
    ppd_wavelengths : list
        A list of wavelength strings describing the wavelengths at which
        the detection in a Class II protoplanetary disk was made.  
        Choose only from among this list:
            "cm"
            "mm"
            "sub-mm"
            "IR"
            "Vis"
            "UV"
    ppd_d_refs : str
        The literature references which describe the detection of the molecule in a 
        Class II protoplanetary disk.  Multiple references should be separated by a 
        semi-colon (';'), which should not otherwise be used in the string.
    ppd_l_refs : str
        The literature references which describe the laboratory work that enabled the 
        detection of the molecule in a Class II protoplanetary disk.  Multiple references 
        should be separated by a semi-colon (';'), which should not otherwise be used 
        in the string.
    ppd_d_bib_ids : list
        A list of strings, with each string a citekey for a bibtex entry corresponding
        to literature references which describe the detection of the molecule in a Class II 
        protoplanetary disk.
    ppd_l_bib_ids : list
        A list of strings, with each string a citekey for a bibtex entry corresponding
        to literature references which describe the laboratory work that enabled
        the detection of the molecule in a Class II protoplanetary disk.  Not currently used.
    ppd_isos : list
        A list of molecule objects, each one describing an isotopically substituted species
        detected in a Class II protoplanetary disk.
    ppd_isos_refs : str
        Should start with the isotopologue formula in [], followed by the literature references
        describing the detection of the isotopically substituted species in a Class II 
        protoplanetary disk.

    exgal : bool
        Flag to indicate if the molecule has been detected in an external galaxy 
        (default is False)
    exgal_sources : list
        A list of source objects in which the molecule has been detected in an 
        external galaxy.
    exgal_telescopes : list
        A list of telescope objects that were used to perform the detections of the
        molecule in an external galaxy.
    exgal_wavelengths : list
        A list of wavelength strings describing the wavelengths at which
        the detection in an external galaxy was made.  
        Choose only from among this list:
            "cm"
            "mm"
            "sub-mm"
            "IR"
            "Vis"
            "UV"
    exgal_d_refs : str
        The literature references which describe the detection of the molecule in an 
        external galaxy.  Multiple references should be separated by a 
        semi-colon (';'), which should not otherwise be used in the string.
    exgal_l_refs : str
        The literature references which describe the laboratory work that enabled the 
        detection of the molecule in an external galaxy.  Multiple references 
        should be separated by a semi-colon (';'), which should not otherwise be used 
        in the string.
    exgal_d_bib_ids : list
        A list of strings, with each string a citekey for a bibtex entry corresponding
        to literature references which describe the detection of the molecule in an 
        external galaxy.
    exgal_l_bib_ids : list
        A list of strings, with each string a citekey for a bibtex entry corresponding
        to literature references which describe the laboratory work that enabled
        the detection of the molecule in an external galaxy.  Not currently used.
    exgal_isos : list
        A list of molecule objects, each one describing an isotopically substituted species
        detected in an external galaxy. 

    exo : bool
        Flag to indicate if the molecule has been detected in an exoplanetary atmosphere 
        (default is False)
    exo_sources : list
        A list of source objects in which the molecule has been detected in an 
        exoplanetary atmosphere.
    exo_telescopes : list
        A list of telescope objects that were used to perform the detections of the
        molecule in an exoplanetary atmosphere.
    exo_wavelengths : list
        A list of wavelength strings describing the wavelengths at which
        the detection in an exoplanetary atmosphere was made.  
        Choose only from among this list:
            "cm"
            "mm"
            "sub-mm"
            "IR"
            "Vis"
            "UV"
    exo_d_refs : str
        The literature references which describe the detection of the molecule in an 
        exoplanetary atmosphere.  Multiple references should be separated by a 
        semi-colon (';'), which should not otherwise be used in the string.
    exo_l_refs : str
        The literature references which describe the laboratory work that enabled the 
        detection of the molecule in an exoplanetary atmosphere.  Multiple references 
        should be separated by a semi-colon (';'), which should not otherwise be used 
        in the string.
    exo_d_bib_ids : list
        A list of strings, with each string a citekey for a bibtex entry corresponding
        to literature references which describe the detection of the molecule in an 
        exoplanetary atmosphere.
    exo_l_bib_ids : list
        A list of strings, with each string a citekey for a bibtex entry corresponding
        to literature references which describe the laboratory work that enabled
        the detection of the molecule in an exoplanetary atmosphere.  Not currently used.
    exo_isos : list
        A list of molecule objects, each one describing an isotopically substituted species
        detected in an exoplanetary atmosphere. 

    notes : str
        Any miscellaneous notes about the molecule
    census_version : str
        The version of the census when this entry was last updated.  
        The format is very specific: XXXX.Y.Z
            XXXX    :   the year of the most recent published census
            Y       :   reset to 0 with each published census and incremented by 1 whenever 
                        a molecule or molecules are added to an incremental update 
                        to astromol between published census releases
            Z       :   reset to 0 with each molecule addition and incremented by 1 whenever 
                        an update is made to astromol that does not involve a new 
                        molecule addition.
	change_log : dict
		A dictionary with keys corresponding to census_versions and entries therein
		describing updates to the entry.                        

    Properties
    ----------
    atoms
        The elemental composition of the molecule
    charge
        The charge of the molecule
    cation
        Whether the molecule is a cation or not
    anion
        Whether the molecule is an anion or not
    neutral
        Whether the molecule is neutral or not
    natoms
        The total number of atoms in the molecule
    mass
        The mass of the molecule in amu for the main isotopic species
    nelectrons
        The total number of electrons in the molecule
    radical
        Whether the molecule is a radical or not
    du
        The degree of unsaturation of the molecule
    maxdu
        The maximumd degree of unsaturation possible for a molecule with the
        number of heavy atoms this molecule has
    kappa
        The Ray's asymmetry parameter of this molecule
    
    Methods
    -------
    inspect()
        Prints out a list to the terminal all of the attributes for this molecule and their values.
    refs()
        Prints a nicely formatted list of references and notes to the terminal.
    summary()
        Prints a summary of the information in the database to the terminal

    """



    def __init__(
        self,
        name=None,
        formula=None,
        smiles=None,
        structure_params=None,
        table_formula=None,
        label=None,
        astromol_name=None,  

        year=None,
        
        sources=None,
        telescopes=None,
        wavelengths=None,

        cyclic=False,
        fullerene=False,
        pah=False,

        Acon=None,
        Bcon=None,
        Ccon=None,

        mua=None,
        mub=None,
        muc=None,

        d_refs=None,
        l_refs=None,
        d_ref_bib_ids=None,
        l_ref_bib_ids=None, 

        isotopologues=None,
        isos_d_refs=None,
        isos_l_refs=None,

        ice=False,
        ice_sources=None,
        ice_telescopes=None,
        ice_wavelengths=None,
        ice_d_refs=None,
        ice_l_refs=None,
        ice_d_bib_ids=None,
        ice_l_bib_ids=None, 
        ice_isos=None,       

        ppd=None,
        ppd_sources=None,
        ppd_telescopes=None,
        ppd_wavelengths=None,
        ppd_d_refs=None,
        ppd_l_refs=None,
        ppd_d_bib_ids=None,
        ppd_l_bib_ids=None,
        ppd_isos=None,
        ppd_isos_refs=None,

        exgal=None,
        exgal_sources=None,        
        exgal_telescopes=None,
        exgal_wavelengths=None,
        exgal_d_refs=None,
        exgal_l_refs=None,
        exgal_d_bib_ids=None,
        exgal_l_bib_ids=None,
        exgal_isos=None,

        exo=None,
        exo_sources=None,
        exo_telescopes=None,
        exo_wavelengths=None,
        exo_d_refs=None,
        exo_l_refs=None,
        exo_d_bib_ids=None,
        exo_l_bib_ids=None,
        exo_isos=None,

        notes=None,
        census_version=None,
        change_log=None,
    ):
        """
        Parameters
        ----------
        name : str
            The full name of the molecule
        formula : str
            The chemical formula for the molecule.  Should be interpretable
            by the latex mhchem package.  Does not include isomer or structural
            identifiers such as 'l' or 'c' for linear or cyclic.
        smiles : str
            The smiles string for the molecule. Not currently used.
        structure_params : dict
            A dictionary containing default values to pass to the structure function.
            keys and defaults are:
                showcarbon : bool (default is False)           
        table_formula : str
            The same as formula, but can include isomer or structural identifiers
            for listing in the Tables of Molecules.
        label : str
            The label that will be used for hypertext linking in LaTeX.  Labels
            must be unique among the database.    
        astromol_name : str
            A string representation of the variable representing this molecule
            in astromol

        year : int
            The year the molecule was first detected in space

        sources : list
            A list of source objects in which the molecule is detected
        telescopes : list
            A list of telescope objects that were used to perform the detections
        wavelengths : list
            A list of wavelength strings describing the wavelengths at which
            the detection was made.  Choose only from among this list:
                "cm"
                "mm"
                "sub-mm"
                "IR"
                "Vis"
                "UV"

        cyclic : bool
            A flag used to indicate if the molecule contains a cyclic moiety
            within its structure (default is False)
        fullerene : bool
            A flag used to indicate if the molecule is a fullerene (default is False). 
            Fullerene and pah should not both be True. 
        pah : bool
            A flag used to indicate if the molecule is a polycyclic aromatic
            hydrocarbon (default is False).  Fullerene and pah should not both be True. 

        Acon : int
            The A rotational constant of the molecule in units of MHz
        Bcon : int
            The B rotational constant of the molecule in units of MHz
        Ccon : int
            The C rotational constant of the molecule in units of MHz

        mua : float
            The A-component of the dipole moment of the molecule in units of Debye
        mub : float
            The B-component of the dipole moment of the molecule in units of Debye
        muc : float
            The C-component of the dipole moment of the molecule in units of Debye

        d_refs : str
            The literature references which describe the ISM/CSM detection of the molecule.
            Multiple references should be separated by a semi-colon (';'), which should
            not otherwise be used in the string.
        l_refs : str
            The literature references which describe the laboratory work that enabled
            the ISM/CSM detection of the molecule.  Multiple references should be separated
            by a semi-colon (';'), which should not otherwise be used in the string.        
        d_ref_bib_ids : list
            A list of strings, with each string a citekey for a bibtex entry corresponding
            to literature references which describe the ISM/CSM detection of the molecule.
            Not currently used.
        l_ref_bib_ids : list
            A list of strings, with each string a citekey for a bibtex entry corresponding
            to literature references which describe the laboratory work that enabled
            the ISM/CSM detection of the molecule.  Not currently used.

        isotopologues : str
            A comma-separated list of isotopically substituted species that have been detected.
            This will be updated to be a list of molecule objects instead in a future version.
        isos_d_refs : str
            Should start with the isotopologue formula in [], followed by the literature references
            describing the ISM/CSM detection of the molecule.
        isos_l_refs : str
            Should start with the isotopologue formula in [], followed by the literature references
            describing the laboratory work that enabled the ISM/CSM detection of the molecule.

        ice : bool
            Flag to indicate if the molecule has been detected in an ice (default is False)
        ice_sources : list
            A list of source objects in which the molecule has been detected in an ice.
        ice_telescopes : list
            A list of telescope objects that were used to perform the detections of the
            molecule in an ice.
        ice_wavelengths : list
            A list of wavelength strings describing the wavelengths at which
            the detection in an ice was made.  Choose only from among this list:
                "cm"
                "mm"
                "sub-mm"
                "IR"
                "Vis"
                "UV"
        ice_d_refs : str
            The literature references which describe the detection of the molecule in an ice.
            Multiple references should be separated by a semi-colon (';'), which should
            not otherwise be used in the string.
        ice_l_refs : str
            The literature references which describe the laboratory work that enabled the 
            detection of the molecule in an ice.  Multiple references should be separated by 
            a semi-colon (';'), which should not otherwise be used in the string.
        ice_d_bib_ids : list
            A list of strings, with each string a citekey for a bibtex entry corresponding
            to literature references which describe the detection of the molecule in an ice.
            Not currently used.
        ice_l_bib_ids : list
            A list of strings, with each string a citekey for a bibtex entry corresponding
            to literature references which describe the laboratory work that enabled
            the detection of the molecule in an ice.  Not currently used.
        ice_isos : list
            A list of molecule objects, each one describing an isotopically substituted species
            detected in an ice.     

        ppd : bool
            Flag to indicate if the molecule has been detected in a Class II 
            protoplanetary disk (default is False)
        ppd_sources : list
            A list of source objects in which the molecule has been detected in a Class II 
            protoplanetary disk.
        ppd_telescopes : list
            A list of telescope objects that were used to perform the detections of the
            molecule in a Class II protoplanetary disk.
        ppd_wavelengths : list
            A list of wavelength strings describing the wavelengths at which
            the detection in a Class II protoplanetary disk was made.  
            Choose only from among this list:
                "cm"
                "mm"
                "sub-mm"
                "IR"
                "Vis"
                "UV"
        ppd_d_refs : str
            The literature references which describe the detection of the molecule in a 
            Class II protoplanetary disk.  Multiple references should be separated by a 
            semi-colon (';'), which should not otherwise be used in the string.
        ppd_l_refs : str
            The literature references which describe the laboratory work that enabled the 
            detection of the molecule in a Class II protoplanetary disk.  Multiple references 
            should be separated by a semi-colon (';'), which should not otherwise be used 
            in the string.
        ppd_d_bib_ids : list
            A list of strings, with each string a citekey for a bibtex entry corresponding
            to literature references which describe the detection of the molecule in a Class II 
            protoplanetary disk.
        ppd_l_bib_ids : list
            A list of strings, with each string a citekey for a bibtex entry corresponding
            to literature references which describe the laboratory work that enabled
            the detection of the molecule in a Class II protoplanetary disk.  Not currently used.
        ppd_isos : list
            A list of molecule objects, each one describing an isotopically substituted species
            detected in a Class II protoplanetary disk.
        ppd_isos_refs : str
            Should start with the isotopologue formula in [], followed by the literature references
            describing the detection of the isotopically substituted species in a Class II 
            protoplanetary disk.

        exgal : bool
            Flag to indicate if the molecule has been detected in an external galaxy 
            (default is False)
        exgal_sources : list
            A list of source objects in which the molecule has been detected in an 
            external galaxy.
        exgal_telescopes : list
            A list of telescope objects that were used to perform the detections of the
            molecule in an external galaxy.
        exgal_wavelengths : list
            A list of wavelength strings describing the wavelengths at which
            the detection in an external galaxy was made.  
            Choose only from among this list:
                "cm"
                "mm"
                "sub-mm"
                "IR"
                "Vis"
                "UV"
        exgal_d_refs : str
            The literature references which describe the detection of the molecule in an 
            external galaxy.  Multiple references should be separated by a 
            semi-colon (';'), which should not otherwise be used in the string.
        exgal_l_refs : str
            The literature references which describe the laboratory work that enabled the 
            detection of the molecule in an external galaxy.  Multiple references 
            should be separated by a semi-colon (';'), which should not otherwise be used 
            in the string.
        exgal_d_bib_ids : list
            A list of strings, with each string a citekey for a bibtex entry corresponding
            to literature references which describe the detection of the molecule in an 
            external galaxy.
        exgal_l_bib_ids : list
            A list of strings, with each string a citekey for a bibtex entry corresponding
            to literature references which describe the laboratory work that enabled
            the detection of the molecule in an external galaxy.  Not currently used.
        exgal_isos : list
            A list of molecule objects, each one describing an isotopically substituted species
            detected in an external galaxy. 

        exo : bool
            Flag to indicate if the molecule has been detected in an exoplanetary atmosphere 
            (default is False)
        exo_sources : list
            A list of source objects in which the molecule has been detected in an 
            exoplanetary atmosphere.
        exo_telescopes : list
            A list of telescope objects that were used to perform the detections of the
            molecule in an exoplanetary atmosphere.
        exo_wavelengths : list
            A list of wavelength strings describing the wavelengths at which
            the detection in an exoplanetary atmosphere was made.  
            Choose only from among this list:
                "cm"
                "mm"
                "sub-mm"
                "IR"
                "Vis"
                "UV"
        exo_d_refs : str
            The literature references which describe the detection of the molecule in an 
            exoplanetary atmosphere.  Multiple references should be separated by a 
            semi-colon (';'), which should not otherwise be used in the string.
        exo_l_refs : str
            The literature references which describe the laboratory work that enabled the 
            detection of the molecule in an exoplanetary atmosphere.  Multiple references 
            should be separated by a semi-colon (';'), which should not otherwise be used 
            in the string.
        exo_d_bib_ids : list
            A list of strings, with each string a citekey for a bibtex entry corresponding
            to literature references which describe the detection of the molecule in an 
            exoplanetary atmosphere.
        exo_l_bib_ids : list
            A list of strings, with each string a citekey for a bibtex entry corresponding
            to literature references which describe the laboratory work that enabled
            the detection of the molecule in an exoplanetary atmosphere.  Not currently used.
        exo_isos : list
            A list of molecule objects, each one describing an isotopically substituted species
            detected in an exoplanetary atmosphere. 

        notes : str
            Any miscellaneous notes about the molecule
        census_version : str
            The version of the census when this entry was last updated.  
            The format is very specific: XXXX.Y.Z
                XXXX    :   the year of the most recent published census
                Y       :   reset to 0 with each published census and incremented by 1 whenever 
                            a molecule or molecules are added to an incremental update 
                            to astromol between published census releases
                Z       :   reset to 0 with each molecule addition and incremented by 1 whenever 
                            an update is made to astromol that does not involve a new 
                            molecule addition.
    	change_log : dict
    		A dictionary with keys corresponding to census_versions and entries therein
    		describing updates to the entry.
        """

        self.name = name
        self.formula = formula
        self.smiles = smiles
        self.structure_params = structure_params
        self.table_formula = table_formula
        self.label = label
        self.astromol_name = astromol_name

        self.year = year
        
        self.sources = sources
        self.telescopes = telescopes
        self.wavelengths = wavelengths

        self.cyclic = cyclic
        self.fullerene = fullerene
        self.pah = pah

        self.Acon = Acon
        self.Bcon = Bcon
        self.Ccon = Ccon

        self.mua = mua
        self.mub = mub
        self.muc = muc

        self.d_refs = d_refs
        self.l_refs = l_refs
        self.d_ref_bib_ids = d_ref_bib_ids
        self.l_ref_bib_ids = l_ref_bib_ids

        self.isotopologues = isotopologues
        self.isos_d_refs = isos_d_refs
        self.isos_l_refs = isos_l_refs

        self.ice = ice
        self.ice_sources = ice_sources
        self.ice_telescopes = ice_telescopes
        self.ice_wavelengths = ice_wavelengths
        self.ice_d_refs = ice_d_refs
        self.ice_l_refs = ice_l_refs
        self.ice_d_bib_ids = ice_d_bib_ids
        self.ice_l_bib_ids = ice_l_bib_ids
        self.ice_isos = ice_isos

        self.ppd = ppd
        self.ppd_sources = ppd_sources
        self.ppd_telescopes = ppd_telescopes
        self.ppd_wavelengths = ppd_wavelengths
        self.ppd_d_refs = ppd_d_refs
        self.ppd_l_refs = ppd_l_refs
        self.ppd_d_bib_ids = ppd_d_bib_ids
        self.ppd_l_bib_ids = ppd_l_bib_ids
        self.ppd_isos = ppd_isos 
        self.ppd_isos_refs = ppd_isos_refs       

        self.exgal = exgal
        self.exgal_sources = exgal_sources
        self.exgal_telescopes = exgal_telescopes
        self.exgal_wavelengths = exgal_wavelengths
        self.exgal_d_refs = exgal_d_refs
        self.exgal_l_refs = exgal_l_refs
        self.exgal_d_bib_ids = exgal_d_bib_ids
        self.exgal_l_bib_ids = exgal_l_bib_ids
        self.exgal_isos = exgal_isos

        self.exo = exo
        self.exo_sources = exo_sources
        self.exo_telescopes = exo_telescopes
        self.exo_wavelengths = exo_wavelengths
        self.exo_d_refs = exo_d_refs
        self.exo_l_refs = exo_l_refs
        self.exo_d_bib_ids = exo_d_bib_ids
        self.exo_l_bib_ids = exo_l_bib_ids
        self.exo_isos = exo_isos       

        self.notes = notes
        self.census_version = census_version
        self.change_log = change_log

        #self._set_structure_defaults()

    @property
    def atoms(self):
        """
        The elemental composition of the molecule

        Returns
        -------
        dict : the elemental composition of the molecule
        """
        # remove any charge from the molecule if present
        if self.formula[-1] == "+" or self.formula[-1] == "-":
            formula = self.formula[:-1]
        else:
            formula = self.formula
        # split the elements list on capital letters
        split_formula = re.findall("[A-Z][^A-Z]*", formula)
        # fill these into a dictionary
        atoms = {}
        for x in split_formula:
            # break up the element and the number (if present)
            y = [i for i in re.split(r"([a-zA-Z]+)([0-9]+)", x) if i]
            # if no number present, add a "1", then convert to an integer
            if len(y) == 1:
                y.append("1")
            y[1] = int(y[1])
            # add to the dictionary count, creating the element if it doesn't exist already
            if y[0] in atoms:
                atoms[y[0]] += y[1]
            else:
                atoms[y[0]] = y[1]
        return atoms

    @property
    def charge(self):
        """
        The charge of the molecule

        Returns
        -------
        int : the charge of the molecule
        """
        charge_dict = {"+": 1, "-": -1}
        if self.formula[-1] in charge_dict:
            return charge_dict[self.formula[-1]]
        else:
            return 0

    @property
    def cation(self):
        """
        Whether the molecule is a cation

        Returns
        -------
        bool : whether the molecule is a cation
        """        
        if self.charge > 0:
            return True
        else:
            return False

    @property
    def anion(self):
        """
        Whether the molecule is an anion

        Returns
        -------
        bool : whether the molecule is an anion
        """    
        if self.charge < 0:
            return True
        else:
            return False

    @property
    def neutral(self):
        """
        Whether the molecule is neutral

        Returns
        -------
        bool : whether the molecule is neutral
        """   
        if self.charge == 0:
            return True
        else:
            return False

    @property
    def natoms(self):
        """
        The number of atoms in the molecule

        Returns
        -------
        int : the number of atoms in the molecule
        """         
        return np.sum([self.atoms[x] for x in self.atoms])

    @property
    def mass(self):
        """
        The mass of the molecule in amu for the main isotopologue

        Returns
        -------
        int : The mass of the molecule in amu for the main isotopologue
        """   
        return np.sum([masses[x] * self.atoms[x] for x in self.atoms])

    @property
    def nelectrons(self):
        """
        The number of electrons in the molecule

        Returns
        -------
        int : the number of electrons in the molecule
        """   
        return np.sum([electrons[x] * self.atoms[x] for x in self.atoms]) - self.charge

    @property
    def radical(self):
        """
        Whether the molecule is a radical

        Determined by determining if there's an odd number of electrons in the molecule.
        Not capable of identifying biradicals.

        Returns
        -------
        bool : whether the molecule is a radical
        """         
        if (self.nelectrons % 2) == 0:
            return False
        else:
            return True

    @property
    def du(self):
        """
        The degree of unsaturation of the molecule

        Only calculated for molecules containing H, C, N, O, Cl, S, D, S, or F and no
        other elements.

        Returns
        -------
        float : the degree of unsaturation of the molecule
        """         
        # We only calculate a du if the molecule contains exclusively H, C, N, O, Cl, S, D, S, or F
        inclusion = [
            "D",
            "C",
            "H",
            "N",
            "O",
            "Cl",
            "F",
            "S",
        ]
        atoms = list(self.atoms.keys())
        remaining = list(filter(lambda x: x not in inclusion, atoms))
        if len(remaining) != 0:
            pass
        else:
            # assuming we're still here, we then calculate and return du
            return 1 + 0.5 * (
                (self.atoms["H"] if "H" in self.atoms else 0) * -1
                + (self.atoms["C"] if "C" in self.atoms else 0) * 2
                + (self.atoms["N"] if "N" in self.atoms else 0) * 1
                + (self.atoms["Cl"] if "Cl" in self.atoms else 0) * -1
                + (self.atoms["F"] if "F" in self.atoms else 0) * -1
            )

    @property
    def maxdu(self):
        """
        The maximum degree of unsaturation of a molecule with the same heavy atoms as
        the molecule

        Only calculated for molecules containing H, C, N, O, Cl, S, D, S, or F and no
        other elements.

        Returns
        -------
        float : the maximum degree of unsaturation of a molecule with the same heavy atoms as
                the molecule
        """           
        # We only calculate a maxdu if the molecule contains exclusively H, C, N, O, Cl, D, S, or F
        inclusion = [
            "D",
            "C",
            "H",
            "N",
            "O",
            "Cl",
            "F",
            "S",
        ]
        atoms = list(self.atoms.keys())
        remaining = list(filter(lambda x: x not in inclusion, atoms))
        if len(remaining) != 0:
            return None

        # assuming we're still here, we then calculate and return maxdu
        return 1 + 0.5 * (
            (self.atoms["C"] if "C" in self.atoms else 0) * 2 + (self.atoms["N"] if "N" in self.atoms else 0) * 1
        )

    @property
    def kappa(self):
        """
        The Ray's asymmetry parameter of the molecule

        Returns
        -------
        float : the Ray's asymmetry parameter of the molecule
        """         
        # only calculate if at least one of the rotational constants is not None
        if any(x is not None for x in [self.Acon, self.Bcon, self.Ccon]):
            if self.Acon == None and self.Ccon == None:
                return -1
            else:
                return (2 * self.Bcon - self.Acon - self.Ccon) / (self.Acon - self.Ccon)
        else:
            return None

    def inspect(self):
        """
        Prints out a list to the terminal all of the attributes for this molecule and their values.

        Note that this will be a largely unformatted output, and is really only intended to be used
        when a detailed look at what information is in the database for a given species is desired.
        """        

        # Get the attributes and values as a dictionary
        attr_dict = self.__dict__
        # We'll do something special for these so they aren't printed as objects
        special_cases = ['sources','telescopes','ppd_isos']

        for attr in attr_dict:
            if attr not in special_cases:
                # only print out values that aren't None, to keep table easy to read
                print(f"{attr:20}: {attr_dict[attr] if attr_dict[attr] is not None else ''}")
            elif attr == 'sources':
                print(f"{'sources':20}: {[x.name for x in attr_dict[attr]] if attr_dict[attr] is not None else ''}")
            elif attr == 'telescopes':
                print(f"{'telescopes':20}: {[x.shortname for x in attr_dict[attr]] if attr_dict[attr] is not None else ''}")
            elif attr == 'ppd_isos':
                print(f"{'ppd_isos':20}: {[x.formula for x in attr_dict[attr]] if attr_dict[attr] is not None else ''}")        
    
    def refs(self):
        """
        Prints a nicely formatted list of references and notes to the terminal.
        """

        l_refs = self.l_refs.split(";")
        d_refs = self.d_refs.split(";")

        if self.notes != None:
            notes = self.notes.strip("*")

        print("Detection Reference(s)")
        for x in range(len(d_refs)):
            print("[{}] {}".format(x + 1, d_refs[x].strip()))

        print("\nLaboratory Reference(s)")
        for x in range(len(l_refs)):
            print("[{}] {}".format(x + 1, l_refs[x].strip()))

        if self.notes != None:
            print("\nNotes")
            print(notes)

        if self.isotopologues != None:
            iso_d_refs = self.isos_d_refs.split("[")
            del iso_d_refs[0]

            print("\nIsotopologue Detection Reference(s)")
            for x in iso_d_refs:
                print("[" + x.strip())

        if self.ice == True or self.ice == "Tentative":
            print("\nIce Reference(s)")
            print("[Det] {}".format(self.ice_d_refs))
            print("[Lab] {}".format(self.ice_l_refs))

        if self.ppd == True or self.ppd == "Tentative":
            print("\nProtoplanetary Disks Reference(s)")
            print("[{}] {}".format(self.formula, self.ppd_d_refs))
            if self.ppd_isos != None:
                for x in self.ppd_isos:
                    print(f"[{x.formula}] {x.ppd_d_refs}")

        if self.exgal == True or self.exgal == "Tentative":
            print("\nExternal Galaxies Reference(s)")
            print("[{}] {}".format(self.formula, self.exgal_d_refs))

        if self.exo == True or self.exo == "Tentative":
            print("\nExoplanetary Atmospheres Reference(s)")
            print("[{}] {}".format(self.formula, self.exo_d_refs))
            if self.exo_isos != None:
            	for x in self.exo_isos:
            		print(f"[{x.formula}] {x.exo_d_refs}")

    def summary(self):
        """
        Prints a summary of the information in the database for the molecule to the terminal.
        """

        n_dash = len(self.name) + len(self.formula) + 3
        dashes = "-" * n_dash
        print("\n" + dashes)
        print("{} ({})".format(self.name, self.formula))
        print(dashes + "\n")
        print("Atoms:\t{}".format(self.natoms))
        print("Mass:\t{} amu".format(self.mass))
        print("Year Detected:\t{}".format(self.year))
        sources = [x.name for x in self.sources]
        sources_str = ", ".join(sources)
        print("Source(s):\t{}".format(sources_str))
        scopes = [x.shortname for x in self.telescopes]
        scopes_str = ", ".join(scopes)
        print("Telescope(s) Used:\t{}".format(scopes_str))

        attr_str = ""

        if self.neutral == True:
            attr_str += "Neutral, "
        if self.cation == True:
            attr_str += "Cation, "
        if self.anion == True:
            attr_str += "Anion, "
        if self.cyclic == True:
            attr_str += "Cyclic"
        if self.radical == True:
            attr_str += "Radical, "

        attr_str = attr_str.strip(" ").strip(",")

        print("Attributes:\t{}\n".format(attr_str))

        if self.isotopologues != None:
            print("Known Isotopologues:\t{}\n".format(self.isotopologues))

        other_envs = [self.ice, self.ppd, self.exgal, self.exo]

        if any(other_envs) == True:
            other_str = ""
            if self.ice == True:
                other_str += "Ices, "
            if self.ice == "Tentative":
                other_str += "Ices (Tentative), "
            if self.ppd == True:
                other_str += "Protoplanetary Disks, "
            if self.ppd == "Tentative":
                other_str += "Protoplanetary Disks (Tentative), "
            if self.exgal == True:
                other_str += "External Galaxies, "
            if self.exgal == "Tentative":
                other_str += "External Galaxies (Tentative), "
            if self.exo == True:
                other_str += "Exoplanetary Atmospheres, "
            if self.exo == "Tentative":
                other_str += "Exoplanetary Atmospheres (Tentative), "

            other_str = other_str.strip().strip(",") + "\n"

            print("Also Detected In:\t{}".format(other_str))

        if self.exgal == True or self.exgal == "Tentative":
            print("Sources of External Galaxy Detections:\t{}\n".format(self.exgal_sources))

        if self.ppd_isos != None:
            print("Isotopologues Also Detected in Protoplanetary Disks:\t{}\n".format(', '.join([x.formula for x in self.ppd_isos])))

        if self.exo_isos != None:
            print("Isotopologues Also Detected in Exoplanetary Atmospheres:\t{}\n".format(', '.join([x.formula for x in self.exo_isos])))		

        self.refs()

#     def _set_structure_defaults(self):
# 
#         defaults = {
#             "showcarbon" : False,
#         }
# 
#         if self.structure_params is None:
#             self.structure_params = {}
#         for param in defaults:
#             if param not in self.structure_params:
#                 self.structure_params[param] = defaults[param]
#     
#     def structure(
#         self,
#         showcarbon=None,
#         ):
#         if showcarbon is None:
#             showcarbon = self.structure_params["showcarbon"]
#         if self.smiles:
#             mol = Chem.MolFromSmiles(self.smiles)
#             if showcarbon is True:
#                 for atom in mol.GetAtoms():
#                     atom.SetProp('atomLabel', atom.GetSymbol())
#             return mol
#         else:
#             pass

    # for potentially labeling explicit things
    # for atom in mol.GetAtoms():
    # symbol = atom.GetSymbol()
    # if atom.GetNumImplicitHs() > 0:
    #     nhs = atom.GetNumImplicitHs()
    # if nhs == 1:
    #     my_str = symbol + "H"
    # elif nhs > 1:
    #     my_str = symbol + "H" + "<sub>" + str(nhs) + "</sub>"
    # else:
    #     my_str = symbol
    # atom.SetProp('atomLabel', my_str)
            
electrons = {
    "H": 1,
    "D": 1,
    "He": 2,
    "C": 6,
    "N": 7,
    "O": 8,
    "F": 9,
    "Na": 11,
    "Mg": 12,
    "Al": 13,
    "Si": 14,
    "P": 15,
    "S": 16,
    "Cl": 17,
    "Ar": 18,
    "K": 19,
    "Ca": 20,
    "Ti": 22,
    "V": 23,
    "Fe": 26,
}

masses = {
    "H": 1,
    "D": 2,
    "He": 4,
    "C": 12,
    "N": 14,
    "O": 16,
    "F": 19,
    "Na": 23,
    "Mg": 24,
    "Al": 27,
    "Si": 28,
    "P": 31,
    "S": 32,
    "Cl": 35,
    "Ar": 36,
    "K": 39,
    "Ca": 40,
    "Ti": 48,
    "V": 51,
    "Fe": 56,
}

######################################################################
#                            Two Atoms                               #
######################################################################

CH = Molecule(
    name="methylidyne",
    formula="CH",
    smiles="[C]-[H]",
    year=1937,
    label="CH",
    astromol_name="CH",
    sources=[LOSCloud],
    telescopes=[MtWilson],
    wavelengths=["UV", "Vis"],
    d_refs="Dunham 1937 PASP 49, 26; Swings & Rosenfeld 1937 ApJ 86, 483; McKellar 1940 PASP 52, 187",
    l_refs="Jevons 1932 Phys Soc. pp 177-179; Brazier & Brown 1983 JCP 78, 1608",
    notes="*First radio in Rydbeck et al. 1973 Nature 246, 466",
    exgal=True,
    exgal_d_refs="Whiteoak et al. 1980 MNRAS 190, 17p",
    exgal_d_bib_ids=["1980MNRAS.190P..17W"],
    exgal_sources="LMC, NGC 4945, NGC 5128",
    Bcon=425476,
    mua=1.5,
    census_version='2018.0.0',
)
CN = Molecule(
    name="cyano radical",
    formula="CN",
    smiles="[C]#[N]",
    year=1940,
    label="CN",
    astromol_name="CN",
    sources=[LOSCloud],
    telescopes=[MtWilson],
    wavelengths=["UV"],
    d_refs="McKellar 1940 PASP 52, 187",
    l_refs="Poletto and Rigutti 1965 Il Nuovo Cimento 39, 519; Dixon & Woods 1977 JCP 67, 3956; Thomas & Dalby 1968 Can. J. Phys. 46, 2815",
    notes="*First radio in Jefferts et al. 1970 ApJ 161, L87",
    ppd=True,
    ppd_d_refs="Kastner et al. 1997 Science 277, 67; Dutrey et al. 1997 A&A 317, L55",
    ppd_d_bib_ids = ["1997Sci...277...67K", "1997A&A...317L..55D"],    
    ppd_isos= [
                Molecule(
                        formula="C15N",
                        table_formula=r"C^{15}N",
                        ppd_d_refs = "Hily-Blant et al. 2017 A&A 603, L6",
                        ppd_d_bib_ids = ["2017A&A...603L...6H"],
                )
            ],
    ppd_isos_refs='"[C15N]" Hily-Blant et al. 2017 A&A 603, L6', # soon to be deprecated; here for legacy refs functionality
    exgal=True,
    exgal_d_refs="Henkel et al. 1988 A&A 201, L23",
    exgal_d_bib_ids=["1988A&A...201L..23H"],
    exgal_sources="M82, NGC 253, IC 342",
    Bcon=56693,
    mua=1.5,
    census_version='2018.0.0',
)
CHp = Molecule(
    name="methylidyne cation",
    formula="CH+",
    smiles="[C+][H]",
    year=1941,
    label="CH+",
    astromol_name="CHp",
    sources=[LOSCloud],
    telescopes=[MtWilson],
    wavelengths=["UV", "Vis"],
    d_refs="Douglas & Herzberg 1941 ApJ 94, 381; Dunham 1937 PASP 49, 26",
    l_refs="Douglas & Herzberg 1941 ApJ 94, 381",
    notes=None,
    ppd=True,
    ppd_d_refs="Thi et al. 2011 A&A 530, L2",
    ppd_d_bib_ids=["2011A&A...530L...2T"],
    exgal=True,
    exgal_d_refs="Magain & Gillet 1987 A&A 184, L5",
    exgal_d_bib_ids=["1987A&A...184L...5M"],
    exgal_sources="LMC",
    Bcon=417617,
    mua=1.7,
    census_version='2018.0.0',
)
OH = Molecule(
    name="hydroxyl radical",
    formula="OH",
    smiles="[O][H]",
    year=1963,
    label="OH",
    astromol_name="OH",
    sources=[CasALOS],
    telescopes=[Millstone],
    wavelengths=["cm"],
    d_refs="Weinreb et al. 1963 Nature 200, 829",
    l_refs="Ehrenstein et al. 1959 PRL 3, 40",
    notes=None,
    ppd=True,
    ppd_d_refs="Mandell et al. 2008 ApJ 681, L25; Salyk et al. 2008 ApJ 676, L49",
    ppd_d_bib_ids=["2008ApJ...681L..25M", "2008ApJ...676L..49S"],
    exgal=True,
    exgal_d_refs="Weliachew 1971 ApJ 167, L47",
    exgal_d_bib_ids=["1971ApJ...167L..47W"],
    exgal_sources="M82, NGC 253",
    exo=True,
    exo_d_bib_ids=["2021ApJ...910L...9N"],
    exo_d_refs="Nugroho et al. 2021 ApJ 910, L9",
    Bcon=556174,
    mua=1.7,
    census_version='2018.0.0',
)
CO = Molecule(
    name="carbon monoxide",
    formula="CO",
    smiles="[C]=[O]",
    structure_params={"showcarbon" : True},
    year=1970,
    label="CO",
    astromol_name="CO",
    sources=[Orion],
    telescopes=[NRAO36],
    wavelengths=["mm"],
    d_refs="Wilson et al. 1970 ApJ 161, L43",
    l_refs="Cord et al. 1968 Microwave Spectral Tables V5",
    notes=None,
    ice=True,
    ice_d_refs="Soifer et al. 1979 ApJ 232, L53",
    ice_l_refs="Mantz et al. 1975 JMS 57, 155",
    ice_d_bib_ids=["1979ApJ...232L..53S"],
    ppd=True,
    ppd_d_refs="Beckwith et al. 1986 ApJ 309, 755",
    ppd_d_bib_ids=["1986ApJ...309..755B"],
    ppd_isos=[
                Molecule(
                            formula="13CO",
                            table_formula=r'^{13}CO',
                            ppd_d_refs="Sargent & Beckwith 1987 ApJ 323, 294",
                            ppd_d_bib_ids=["1987ApJ...323..294S"],
                ),
                Molecule(
                            formula="C18O",
                            table_formula=r'C^{18}O',
                            ppd_d_refs="Dutrey et al. 1994 A&A 286, 149",
                            ppd_d_bib_ids=["1994A&A...286..149D"],
                ),
                Molecule(
                            formula="C17O",
                            table_formula=r'C^{17}O',
                            ppd_d_refs="Smith et al. 2009 ApJ 701, 163; Guilloteau et al. 2013 A&A 549, A92",
                            ppd_d_bib_ids=["2009ApJ...701..163S", "2013A&A...549A..92G"],
                ),
    ],
    ppd_isos_refs='"[13CO]" Sargent & Beckwith 1987 ApJ 323, 294 "[C18O]" Dutrey et al. 1994 A&A 286, 149 "[C17O]" Smith et al. 2009 ApJ 701, 163; Guilloteau et al. 2013 A&A 549, A92',
    exo=True,
    exo_d_refs="Madhusudhan et al. 2011 Nature 469, 64; Barman et al. 2011 ApJ 733, 65; Lanotte et al. 2014 A&A 572, A73; Barman et al. 2015 ApJ 804, 61",
    exo_d_bib_ids=["2011Natur.469...64M","2011ApJ...733...65B","2014A&A...572A..73L","2015ApJ...804...61B"],
    exo_isos = [
    			Molecule(
                            formula="13CO",
                            table_formula=r'^{13}CO',
                            exo_d_refs="Zhang et al. 2021 Nature 595, 370",
                            exo_d_bib_ids=["2021Natur.595..370Z"],
                            exo_sources= ["TYC 8998-760-1 b"] #to be updated later to full objects
                ),    
    ],
    exgal=True,
    exgal_d_refs="Rickard et al. 1975 ApJ 199, L75",
    exgal_d_bib_ids=["1975ApJ...199L..75R"],
    exgal_sources="M82, NGC 253",
    Bcon=57636,
    mua=0.1,
    census_version='2021.1.0',
    change_log = {
    				'2018.0.0' : 'First entry',
    				'2021.1.0' : 'Added 13CO detected in exoplanet atmosphere.'
    			}
)
H2 = Molecule(
    name="hydrogen",
    formula="H2",
    year=1970,
    label="H2",
    astromol_name="H2",
    sources=[XiPerLOS],
    telescopes=[Aerobee],
    wavelengths=["UV"],
    d_refs="Carruthers 1970 ApJ 161, L81",
    l_refs="Carruthers 1970 ApJ 161, L81",
    notes=None,
    ppd=True,
    ppd_d_refs="Thi et al. 1999 ApJ 521, L63",
    ppd_d_bib_ids=["1999ApJ...521L..63T"],
    ppd_isos=[Molecule(
                        formula="HD",
                        ppd_d_refs="Bergin et al. 2013 Nature 493, 644",
                        ppd_d_bib_ids=["2013Natur.493..644B"],
                    ),
    ],
    ppd_isos_refs='"[HD]" Bergin et al. 2013 Nature 493, 644',
    exgal=True,
    exgal_d_refs="Thompson et al. 1978 ApJ 222, L49",
    exgal_d_bib_ids=["1978ApJ...222L..49T"],
    exgal_sources="NGC 1068",
    mua=0.0,
    census_version='2018.0.0',
)
SiO = Molecule(
    name="silicon monoxide",
    formula="SiO",
    year=1971,
    label="SiO",
    astromol_name="SiO",
    sources=[SgrB2],
    telescopes=[NRAO36],
    wavelengths=["mm"],
    d_refs="Wilson et al. 1971 ApJ 167, L97",
    l_refs="Törring 1968 Z. Naturforschung 23A, 777; Raymonda et al. 1970 JCP 52, 3458",
    notes=None,
    exgal=True,
    exgal_d_refs="Mauersberger & Henkel 1991 A&A 245, 457",
    exgal_d_bib_ids=["1991A&A...245..457M"],
    exgal_sources="NGC 253",
    Bcon=21712,
    mua=3.1,
    census_version='2018.0.0',
)
CS = Molecule(
    name="carbon monosulfide",
    formula="CS",
    year=1971,
    label="CS",
    astromol_name="CS",
    sources=[Orion, W51, IRC10216, DR21],
    telescopes=[NRAO36],
    wavelengths=["mm"],
    d_refs="Penzias et al. 1971 ApJ 168, L53",
    l_refs="Mockler & Bird 1955 Phys Rev 98, 1837",
    notes=None,
    ppd=True,
    ppd_d_refs="Ohashi et al. 1991 AJ 102, 2054; Blake et al. 1992 ApJ 391, L99; Guilloteau et al. 2012 A&A 548, A70",
    ppd_d_bib_ids=["1991AJ....102.2054O", "1992ApJ...391L..99B", "2012A&A...548A..70G"],
    ppd_isos=[
                Molecule(
                    formula = 'C34S',
                    table_formula = r"C^{34}S",
                    ppd_d_refs="Le Gal et al. 2019 ApJ 876, 72; Loomis et al. 2020 ApJ 893, 101",
                    ppd_d_bib_ids=["2020ApJ...893..101L", "2019ApJ...876...72L"],
                ),
                Molecule(
                    formula = '13CS',
                    table_formula = r"^{13}CS",
                    ppd_d_refs="Le Gal et al. 2019 ApJ 876, 72; Loomis et al. 2020 ApJ 893, 101",
                    ppd_d_bib_ids=["2020ApJ...893..101L", "2019ApJ...876...72L"],
                ),
    ],
    ppd_isos_refs="[C34S] Le Gal et al. 2019 ApJ 876, 72; Loomis et al. 2020 ApJ 893, 101",
    exgal=True,
    exgal_d_refs="Henkel & Bally 1985 A&A 150, L25",
    exgal_d_bib_ids=["1985A&A...150L..25H"],
    exgal_sources="M82, IC 342",
    Bcon=24496,
    mua=2.0,
    census_version='2018.0.0',
)
SO = Molecule(
    name="sulfur monoxide",
    formula="SO",
    year=1973,
    label="SO",
    astromol_name="SO",
    sources=[Orion],
    telescopes=[NRAO36],
    wavelengths=["mm"],
    d_refs="Gottlieb & Ball 1973 ApJ 184, L59",
    l_refs="Winnewisser et al. 1964 JCP 41, 1687",
    notes=None,
    ppd=True,
    ppd_d_refs="Fuente et al. 2010 A&A 524, A19",
    ppd_d_bib_ids=["2010A&A...524A..19F"],
    exgal=True,
    exgal_d_refs="Johansson 1991 Proc. IAU Symposium 146, 1; Petuchowski & Bennett 1992 ApJ 391, 137",
    exgal_d_bib_ids=["1991IAUS..146....1J","1992ApJ...391..137P"],
    exgal_sources="M82, NGC 253",
    Bcon=21524,
    mua=1.5,
    census_version='2018.0.0',
)
SiS = Molecule(
    name="silicon monosulfide",
    formula="SiS",
    year=1975,
    label="SiS",
    astromol_name="SiS",
    sources=[IRC10216],
    telescopes=[NRAO36],
    wavelengths=["mm"],
    d_refs="Morris et al. 1975 ApJ 199, L47",
    l_refs="Hoeft 1965 Z. fur Naturforschung A, 20, 1327",
    notes=None,
    Bcon=9077,
    mua=1.7,
    census_version='2018.0.0',
)
NS = Molecule(
    name="nitrogen monosulfide",
    formula="NS",
    year=1975,
    label="NS",
    astromol_name="NS",
    sources=[SgrB2],
    telescopes=[NRAO36],
    wavelengths=["mm"],
    d_refs="Gottlieb et al. 1975 ApJ 200, L147; Kuiper et al. 1975 ApJ 200, L151",
    l_refs="Amano et al. 1969 JMS 32, 97",
    notes=None,
    exgal=True,
    exgal_d_refs="Martin et al. 2003 A&A 411, L465",
    exgal_d_bib_ids=["2003A&A...411L.465M"],
    exgal_sources="NGC 253",
    Bcon=23155,
    mua=1.8,
    census_version='2018.0.0',
)
C2 = Molecule(
    name="dicarbon",
    formula="C2",
    year=1977,
    label="C2",
    astromol_name="C2",
    sources=[CygnusOB212LOS],
    telescopes=[MtHopkins],
    wavelengths=["IR"],
    d_refs="Souza and Lutz 1977 ApJ 216, L49",
    l_refs="Phillips 1948 ApJ 107, 389",
    exgal=True,
    exgal_d_refs="Welty et al. 2013 MNRAS 428, 1107",
    exgal_d_bib_ids=["2013MNRAS.428.1107W"],
    exgal_sources="SMC",
    notes=None,
    mua=0.0,
    census_version='2018.0.0',
)
NO = Molecule(
    name="nitric oxide",
    formula="NO",
    year=1978,
    label="NO",
    astromol_name="NO",
    sources=[SgrB2],
    telescopes=[NRAO36],
    wavelengths=["mm"],
    d_refs="Liszt and Turner 1978 ApJ 224, L73",
    l_refs="Gallagher & Johnson 1956 Phys Rev 103, 1727",
    notes=None,
    exgal=True,
    exgal_d_refs="Martin et al. 2003 A&A 411, L465",
    exgal_d_bib_ids=["2003A&A...411L.465M"],
    exgal_sources="NGC 253",
    Bcon=50849,
    mua=0.2,
    census_version='2018.0.0',
)
HCl = Molecule(
    name="hydrogen chloride",
    formula="HCl",
    year=1985,
    label="HCl",
    astromol_name="HCl",
    sources=[Orion],
    telescopes=[Kuiper],
    wavelengths=["sub-mm"],
    d_refs="Blake et al. 1985 ApJ 295, 501",
    l_refs="de Lucia et al. 1971 Phys Rev A 3, 1849",
    exgal=True,
    exgal_d_refs="Wallstrom et al. 2019 A&A 629, A128",
    exgal_d_bib_ids=["2019A&A...629A.128W"],
    exgal_sources="PKS 1830-211 LOS",
    notes=None,
    Bcon=312989,
    mua=1.1,
    census_version='2018.0.0',
)
NaCl = Molecule(
    name="sodium chloride",
    formula="NaCl",
    year=1987,
    label="NaCl",
    astromol_name="NaCl",
    sources=[IRC10216],
    telescopes=[IRAM30],
    wavelengths=["mm"],
    d_refs="Cernicharo & Guélin 1987 A&A 183, L10",
    l_refs="Lovas & Tiemann 1974 J Phys Chem Ref Data 3, 609",
    notes=None,
    Bcon=6513,
    mua=9.0,
    census_version='2018.0.0',
)
AlCl = Molecule(
    name="aluminum chloride",
    formula="AlCl",
    year=1987,
    label="AlCl",
    astromol_name="AlCl",
    sources=[IRC10216],
    telescopes=[IRAM30],
    wavelengths=["mm"],
    d_refs="Cernicharo & Guélin 1987 A&A 183, L10",
    l_refs="Lovas & Tiemann 1974 J Phys Chem Ref Data 3, 609",
    notes=None,
    Bcon=7289,
    mua="*",
    census_version='2018.0.0',
)
KCl = Molecule(
    name="potassium chloride",
    formula="KCl",
    year=1987,
    label="KCl",
    astromol_name="KCl",
    sources=[IRC10216],
    telescopes=[IRAM30],
    wavelengths=["mm"],
    d_refs="Cernicharo & Guélin 1987 A&A 183, L10",
    l_refs="Lovas & Tiemann 1974 J Phys Chem Ref Data 3, 609",
    notes=None,
    Bcon=3845,
    mua=10.3,
    census_version='2018.0.0',
)
AlF = Molecule(
    name="aluminum fluoride",
    formula="AlF",
    year=1987,
    label="AlF",
    astromol_name="AlF",
    sources=[IRC10216],
    telescopes=[IRAM30],
    wavelengths=["mm"],
    d_refs="Cernicharo & Guélin 1987 A&A 183, L10",
    l_refs="Lovas & Tiemann 1974 J Phys Chem Ref Data 3, 609",
    notes="*Confirmed in 1994 ApJ 433, 729",
    isotopologues="26AlF",
    isos_d_refs='[26AlF] Kamiński et al. 2018 Nature Astronomy 2, 778',
    Bcon=16488,
    mua=1.5,
    census_version='2018.0.0',
)
PN = Molecule(
    name="phosphorous mononitride",
    formula="PN",
    year=1987,
    label="PN",
    astromol_name="PN",
    sources=[TMC1, Orion, W51, SgrB2],
    telescopes=[NRAOARO12, FCRAO14m, OVRO],
    wavelengths=["mm"],
    d_refs="Sutton et al. 1985 ApJS 58, 341",
    l_refs="Wyse et al. 1972 JCP 57, 1106",
    notes="*Confirmed in Turner & Bally 1987 ApJ 321, L75 and Ziurys 1987 ApJ 321 L81",
    Bcon=23495,
    mua=2.7,
    census_version='2018.0.0',
)
SiC = Molecule(
    name="silicon carbide",
    formula="SiC",
    year=1989,
    label="SiC",
    astromol_name="SiC",
    sources=[IRC10216],
    telescopes=[IRAM30],
    wavelengths=["mm"],
    d_refs="Cernicharo et al. 1989 ApJ 341, L25",
    l_refs="Cernicharo et al. 1989 ApJ 341, L25",
    notes=None,
    Bcon=20298,
    mua=1.7,
    census_version='2018.0.0',
)
CP = Molecule(
    name="carbon monophosphide",
    formula="CP",
    year=1990,
    label="CP",
    astromol_name="CP",
    sources=[IRC10216],
    telescopes=[IRAM30],
    wavelengths=["mm"],
    d_refs="Guélin et al. 1990 A&A 230, L9",
    l_refs="Saito et al. 1989 ApJ 341, 1114",
    notes=None,
    Bcon=23860,
    mua="*",
    census_version='2018.0.0',
)
NH = Molecule(
    name="imidogen radical",
    formula="NH",
    year=1991,
    label="NH",
    astromol_name="NH",
    sources=[XiPerLOS, HD27778LOS],
    telescopes=[KPNO4m, IRAM30],
    wavelengths=["UV"],
    d_refs="Meyer & Roth 1991 ApJ 376, L49",
    l_refs="Dixon 1959 Can J. Phys. 37, 1171 and Klaus et al. 1997 A&A 322, L1",
    notes="*First radio in Cernicharo et al. 2000 ApJ 534, L199",
    exgal=True,
    exgal_d_refs="Gonzalez-Alfonso et al. 2004 ApJ 613, 247",
    exgal_d_bib_ids=["2004ApJ...613..247G"],
    exgal_sources="Arp 220",
    Bcon=489959,
    mua=1.4,
    census_version='2018.0.0',
)
SiN = Molecule(
    name="silicon nitride ",
    formula="SiN",
    year=1992,
    label="SiN",
    astromol_name="SiN",
    sources=[IRC10216],
    telescopes=[NRAOARO12],
    wavelengths=["mm"],
    d_refs="Turner 1992 ApJ 388, L35",
    l_refs="Saito et al. 1983 JCP 78, 6447",
    notes=None,
    Bcon=21828,
    mua=2.6,
    census_version='2018.0.0',
)
SOp = Molecule(
    name="sulfur monoxide cation",
    formula="SO+",
    year=1992,
    label="SO+",
    astromol_name="SOp",
    sources=[IC443G],
    telescopes=[NRAOARO12],
    wavelengths=["mm"],
    d_refs="Turner 1992 ApJ 396, L107",
    l_refs="Amano et al. 1991 JMS 146, 519",
    notes=None,
    exgal=True,
    exgal_d_refs="Muller et al. 2011 A&A 535, A103",
    exgal_d_bib_ids=["2011A&A...535A.103M"],
    exgal_sources="PKS 1830-211 LOS",
    Bcon=23249,
    mua="*",
    census_version='2018.0.0',
)
COp = Molecule(
    name="carbon monoxide cation",
    formula="CO+",
    year=1993,
    label="CO+",
    astromol_name="COp",
    sources=[M17SW, NGC7027],
    telescopes=[NRAOARO12],
    wavelengths=["mm"],
    d_refs="Latter et al. 1993 ApJ 419, L97",
    l_refs="Sastry et al. 1981 ApJ 250, L91",
    notes=None,
    exgal=True,
    exgal_d_refs="Fuente et al. 2006 ApJ 641, L105",
    exgal_d_bib_ids=["2006ApJ...641L.105F"],
    exgal_sources="M82",
    Bcon=58983,
    mua=2.6,
    census_version='2018.0.0',
)
HF = Molecule(
    name="hydrogen fluoride",
    formula="HF",
    year=1997,
    label="HF",
    astromol_name="HF",
    sources=[SgrB2LOS],
    telescopes=[ISO],
    wavelengths=["IR"],
    d_refs="Neufeld et al. 1997 ApJ 488, L141",
    l_refs="Nolt et al. 1987 JMS 125, 274",
    notes=None,
    exgal=True,
    exgal_d_refs="van der Werf et al. 2010 A&A 518, L42; Rangwala et al. 2011 ApJ 743, 94; Monje et al. 2011 ApJL 742, L21",
    exgal_d_bib_ids=["2010A&A...518L..42V", "2011ApJ...743...94R", "2011ApJ...742L..21M"],
    exgal_sources="QSO Mrk 231, Arp 220, Cloverleaf LOS",
    Bcon=616365,
    mua=1.8,
    census_version='2018.0.0',
)
N2 = Molecule(
    name="nitrogen",
    formula="N2",
    year=2004,
    label="N2",
    astromol_name="N2",
    sources=[HD124314LOS],
    telescopes=[FUSE],
    wavelengths=["UV"],
    d_refs="Knauth et al. 2004 Nature 409, 636",
    l_refs="Stark et al. 2000 ApJ 531, 321",
    notes=None,
    mua=0.0,
    census_version='2018.0.0',
)
CFp = Molecule(
    name="fluoromethylidynium cation",
    formula="CF+",
    year=2006,
    label="CF+",
    astromol_name="CFp",
    sources=[OrionBar],
    telescopes=[IRAM30, APEX],
    wavelengths=["mm"],
    d_refs="Neufeld et al. 2006 A&A 454, L37",
    l_refs="Plummer et al. 1986 JCP 84, 2427",
    notes=None,
    exgal=True,
    exgal_d_refs="Muller et al. 2016 A&A 589, L5",
    exgal_d_bib_ids=["2016A&A...589L...5M"],
    exgal_sources="PKS 1830-211 LOS",
    Bcon=51294,
    mua=1.1,
    census_version='2018.0.0',
)
PO = Molecule(
    name="phosphorous monoxide",
    formula="PO",
    year=2007,
    label="PO",
    astromol_name="PO",
    sources=[VYCaMaj],
    telescopes=[SMT10],
    wavelengths=["mm"],
    d_refs="Tenenbaum et al. 2007 ApJ 666, L29",
    l_refs="Bailleux et al. 2002 JMS 216, 465",
    notes=None,
    Bcon=21900,
    mua=1.9,
    census_version='2018.0.0',
)
O2 = Molecule(
    name="oxygen",
    formula="O2",
    year=2007,
    label="O2",
    astromol_name="O2",
    sources=[Orion, rhoOphA],
    telescopes=[Odin, Herschel],
    wavelengths=["mm", "sub-mm"],
    d_refs="Larsson et al. 2007 A&A 466, 999",
    l_refs="Endo & Mizushima 1982 Jpn J Appl Phys 21, L379; Drouin et al. 2010 J Quant Spec Rad Transf 111, 1167",
    notes="*Also Larsson et al. 2007 A&A 466, 999; Tentative in Goldsmith 2002 ApJ 576, 814",
    exgal=True,
    exgal_d_refs="Wang et al. 2020 ApJ 889, 129",
    exgal_d_bib_ids=["2020ApJ...889..129W"],
    exgal_sources="QSO Mrk 231",    
    mua=0.0,
    census_version='2018.0.0',
)
AlO = Molecule(
    name="aluminum monoxide",
    formula="AlO",
    year=2009,
    label="AlO",
    astromol_name="AlO",
    sources=[VYCaMaj],
    telescopes=[SMT10],
    wavelengths=["mm"],
    d_refs="Tenenbaum & Ziurys 2009 ApJ 693, L59",
    l_refs="Yamada et al. 1990, JCP 92, 2146",
    notes=None,
    Bcon=19142,
    mua=4.6,
    census_version='2018.0.0',
)
CNm = Molecule(
    name="cyanide anion",
    formula="CN-",
    year=2010,
    label="CN-",
    astromol_name="CNm",
    sources=[IRC10216],
    telescopes=[IRAM30],
    wavelengths=["mm"],
    d_refs="Agúndez et al. 2010 A&A 517, L2",
    l_refs="Amano 2008 JCP 129, 244305",
    notes=None,
    Bcon=56133,
    mua=0.7,
    census_version='2018.0.0',
)
OHp = Molecule(
    name="hydroxyl cation",
    formula="OH+",
    year=2010,
    label="OH+",
    astromol_name="OHp",
    sources=[SgrB2LOS],
    telescopes=[APEX],
    wavelengths=["sub-mm"],
    d_refs="Wyrowski et al. 2010 A&A 518, A26; Gerin et al. 2010 A&A 518, L110; Benz et al. 2010 A&A 521, L35",
    l_refs="Bekooy et al. 1985 JCP 82, 3868",
    notes=None,
    exgal=True,
    exgal_d_refs="van der Werf et al. 2010 A&A 518, L42; Rangwala et al. 2011 ApJ 743, 94; Gonzalez-Alfonso et al. 2013 A&A 550, A25",
    exgal_d_bib_ids=["2010A&A...518L..42V", "2011ApJ...743...94R", "2013A&A...550A..25G"],
    exgal_sources="QSO Mrk 231, Arp 220, NGC 4418",
    Bcon=492346,
    mua=2.3,
    census_version='2018.0.0',
)
SHp = Molecule(
    name="sulfanylium cation",
    formula="SH+",
    year=2011,
    label="SH+",
    astromol_name="SHp",
    sources=[SgrB2],
    telescopes=[Herschel],
    wavelengths=["sub-mm"],
    d_refs="Benz et al. 2010 A&A 521, L35",
    l_refs="Brown et al. 2009 JMS 255, 68",
    notes="*Also in Menten et al. 2011 A&A 525, A77",
    exgal=True,
    exgal_d_refs="Muller et al. 2017 A&A 606, A109",
    exgal_d_bib_ids=["2017A&A...606A.109M"],
    exgal_sources="PKS 1830-211 LOS",
    Bcon=273810,
    mua=1.3,
    census_version='2018.0.0',
)
HClp = Molecule(
    name="hydrogen chloride cation",
    formula="HCl+",
    year=2012,
    label="HCl+",
    astromol_name="HClp",
    sources=[W31LOS, W49LOS],
    telescopes=[Herschel],
    wavelengths=["sub-mm"],
    d_refs="de Luca et al. 2012 ApJ 751, L37",
    l_refs="Gupta et al. 2012 ApJ 751, L38",
    notes=None,
    Bcon=293444,
    mua=1.8,
    census_version='2018.0.0',
)
SH = Molecule(
    name="mercapto radical",
    formula="SH",
    year=2012,
    label="SH",
    astromol_name="SH",
    sources=[W49LOS],
    telescopes=[SOFIA],
    wavelengths=["sub-mm"],
    d_refs="Neufeld et al. 2012 A&A 542, L6",
    l_refs="Morino & Kawaguchi 1995 JMS 170, 172; Klisch et al. 1996 ApJ 473, 1118",
    notes=None,
    Bcon=283588,
    mua=0.8,
    census_version='2018.0.0',
)
TiO = Molecule(
    name="titanium monoxide",
    formula="TiO",
    year=2013,
    label="TiO",
    astromol_name="TiO",
    sources=[VYCaMaj],
    telescopes=[SMA],
    wavelengths=["mm"],
    d_refs="Kamiński et al. 2013 A&A 551, A113",
    l_refs="Nakimi et al. 1998 JMS 191, 176",
    notes=None,
    exo=True,
    exo_d_refs="Haynes et al. 2015 ApJ 806, 146; Sedaghati et al. 2017 Nature 549, 238; Nugroho et al. 2017 ApJ 154, 221",
    exo_d_bib_ids=["2015ApJ...806..146H","2017Natur.549..238S","2017AJ....154..221N"],
    Bcon=16004,
    mua=3.3,
    census_version='2018.0.0',
)
ArHp = Molecule(
    name="argonium",
    formula="ArH+",
    year=2013,
    label="ArH+",
    astromol_name="ArHp",
    sources=[CrabNebula],
    telescopes=[Herschel],
    wavelengths=["sub-mm"],
    d_refs="Barlow et al. 2013 Science 342, 1343",
    l_refs="Barlow et al. 2013 Science 342, 1343",
    notes=None,
    exgal=True,
    exgal_d_refs="Müller et al. 2015 A&A 582, L4",
    exgal_d_bib_ids=["2015A&A...582L...4M"],
    exgal_sources="PKS 1830-211 LOS",
    Bcon=307966,
    mua=2.2,
    census_version='2018.0.0',
)
NSp = Molecule(
    name="nitrogen sulfide cation",
    formula="NS+",
    year=2018,
    label="NS+",
    astromol_name="NSp",
    sources=[B1b, TMC1, L483],
    telescopes=[IRAM30],
    wavelengths=["mm"],
    d_refs="Cernicharo et al. 2018 ApJL 853, L22",
    l_refs="Cernicharo et al. 2018 ApJL 853, L22",
    notes=None,
    Bcon=25050,
    mua=2.2,
    census_version='2018.0.0',
)
HeHp = Molecule(
    name="helium hydride cation",
    formula="HeH+",
    year=2019,
    label="HeH+",
    astromol_name="HeHp",
    sources=[NGC7027],
    telescopes=[SOFIA],
    wavelengths=["sub-mm"],
    d_refs="Gusten et al. 2019 Nature 568, 357",
    l_refs="Perry et al. 2014 JCP 141, 101101",
    notes="Dipole moment from Engel et al. 2005 MNRAS 357, 471 using equilibrium internuclear separation value of 1.45 from Peyerimhoff 1965 JCP 43, 998",
    Bcon=1006063,
    mua=1.7,
    census_version='2021.0.0',
)
VO = Molecule(
    name="vanadium oxide",
    formula="VO",
    year=2019,
    label="VO",
    astromol_name="VO",
    sources=[VYCaMaj],
    telescopes=[Hubble],
    wavelengths=["IR"],
    d_refs="Humphreys et al. 2019 ApJL 874, L26",
    l_refs="Adam et al. 1995 JMS 170, 94; Cheung et al. 1994 JMS 163, 443",
    notes="Some VO transitions originally observed in the source by Wallerstein & Gonzalez 2001 PASP 113, 954, Wallerstein 1971 ApJ 169, 195, and Wallerstein 1986 ApJ 164, 101, but not assigned as circumstellar until now.  B constant from Cheung et al. 1982 JMS 91, 165.  Dipole moment from Suenram et al. 1991 JMS 148, 114.",
    Bcon=16381,
    mua=3.355,
    census_version='2021.0.0',
)
POp = Molecule(
    name="phosphorus monoxide cation",
    formula="PO+",
    table_formula="PO+",
    year=2022,
    label="PO+",
    astromol_name="POp",
    sources=[G0693],
    telescopes=[IRAM30,Yebes40],
    wavelengths=["cm","mm"],
    d_refs="Rivilla et al. 2022 Frontiers Astron. Space Sci. 9, 829288",
    d_ref_bib_ids = ["2022FrASS...9.9288R"],
    l_refs="Petrmichl et al. 1991 J Chem Phys 94, 3504",
    l_ref_bib_ids = ["1991JChPh..94.3504P"],
    notes = None,
    Bcon=23571,
    mua=3.13,
    census_version='2021.4.0',
    change_log = {
    				'2021.4.0' : 'Initial entry',
    			}
)
SiP = Molecule(
    name="silicon phosphide",
    formula="SiP",
    year=2022,
    label="SiP",
    astromol_name="SiP",
    sources=[IRC10216],
    telescopes=[SMT10,NRAOARO12],
    wavelengths=["mm","sub-mm"],
    d_refs="Koelemay et al. 2022 ApJL 940, L11",
    d_ref_bib_ids = ["Koelemay:2022:L11"],
    l_refs="Koelemay et al. 2022 ApJL 940, L11",
    l_ref_bib_ids = ["Koelemay:2022:L11"],
    notes = "Dipole moment from Ornellas et al. 2000 ApJ 538, 675",
    Bcon=7967,
    mua=0.9,
    census_version='2021.6.0',
    change_log = {
    				'2021.6.0' : 'Initial entry',
    			}
)
######################################################################
#                          Three Atoms                               #
######################################################################


H2O = Molecule(
    name="water",
    formula="H2O",
    year=1969,
    label="H2O",
    astromol_name="H2O",
    sources=[SgrB2, Orion, W49],
    telescopes=[HatCreek],
    wavelengths=["cm"],
    d_refs="Cheung et al. 1969 Nature 221, 626",
    l_refs="Golden et al. 1948 Phys Rev 73, 92",
    notes=None,
    ice=True,
    ice_d_refs="Gillett & Forrest 1973 ApJ 179, 483",
    ice_l_refs="Irvine & Pollack 1968 Icarus 8, 324",
    ice_d_bib_ids=["1973ApJ...179..483G"],
    isotopologues="HDO",
    isos_d_refs='[HDO] Turner et al. 1975 ApJ 198, L125',
    isos_l_refs='[HDO] de Lucia et al. 1974 J Phys Chem Ref Data 3, 211; Erlandsson & Cox 1956 J Chem Phys 25, 778',
    ppd=True,
    ppd_d_refs="Carr et al. 2004 ApJ 603, 213; Hogerheijde et al. 2011 Science 334, 338; Salyk et al. 2008 ApJL 676, L49",
    ppd_d_bib_ids=["2004ApJ...603..213C","2011Sci...334..338H","2008ApJ...676L..49S"],
    exo=True,
    exo_d_refs="Tinetti et al. 2007 Nature 448, 169; Deming et al. 2013 ApJ 774, 95; Kreidberg et al. 2014 ApJL 793, L27; Kreidberg et al. 2015 ApJ 814, 66; Lockwood et al. 2014 ApJ 783, L29",
    exo_d_bib_ids=["2007Natur.448..169T","2013ApJ...774...95D","2014ApJ...793L..27K","2015ApJ...814...66K","2014ApJ...783L..29L"],
    exgal=True,
    exgal_d_refs="Churchwell et al. 1977 A&A 54, 969",
    exgal_d_bib_ids=["1977A&A....54..969C"],
    exgal_sources="M33",
    Acon=835840,
    Bcon=435352,
    Ccon=278139,
    mub=1.9,
    census_version='2018.0.0',
)
HCOp = Molecule(
    name="formylium cation",
    formula="HCO+",
    year=1970,
    label="HCO+",
    astromol_name="HCOp",
    sources=[W3OH, Orion, L134, SgrA, W51],
    telescopes=[NRAO36],
    wavelengths=["mm"],
    d_refs="Buhl & Snyder 1970 Nature 228, 267",
    l_refs="Woods et al. 1975 PRL 35, 1269",
    notes=None,
    ppd=True,
    ppd_d_refs="Kastner et al. 1997 Science 277, 67; Dutrey et al. 1997 A&A 317, L55",
    ppd_d_bib_ids = ["1997Sci...277...67K", "1997A&A...317L..55D"],
    ppd_isos= [
                Molecule (
                            formula = "DCO+",
                            ppd_d_refs = "van Dishoeck et al. 2003 A&A 400, L1",
                            ppd_d_bib_ids = ["2003A&A...400L...1V"],
                ),

                Molecule(
                            formula='H13CO+',
                            table_formula = r"H^{13}CO+",
                            ppd_d_refs = "van Zadelhoff et al. 2001 A&A 377, 566; van Dishoeck et al. 2003 A&A 400, L1",
                            ppd_d_bib_ids = ["2001A&A...377..566V", "2003A&A...400L...1V"]
                ),
        
        
        ],
    ppd_isos_refs='[DCO+] van Dishoeck et al. 2003 A&A 400, L1 [H13CO+] van Zadelhoff et al. 2001 A&A 377, 566; van Dishoeck et al. 2003 A&A 400, L1',
    exgal=True,
    exgal_d_refs="Stark et al. 1979 ApJ 229, 118",
    exgal_d_bib_ids=["1979ApJ...229..118S"],
    exgal_sources="M82",
    Bcon=44594,
    mua=3.9,
    census_version='2018.0.0',
)
HCN = Molecule(
    name="hydrogen cyanide",
    formula="HCN",
    year=1971,
    label="HCN",
    astromol_name="HCN",
    sources=[W3OH, Orion, SgrA, W49, W51, DR21],
    telescopes=[NRAO36],
    wavelengths=["mm"],
    d_refs="Snyder et al. 1971 ApJ 163, L47",
    l_refs="de Lucia & Gordy 1969 Phys Rev 187, 58",
    notes=None,
    ppd=True,
    ppd_d_refs="Kastner et al. 1997 Science 277, 67; Dutrey et al. 1997 A&A 317, L55",
    ppd_d_bib_ids=["1997Sci...277...67K", "1997A&A...317L..55D"],
    ppd_isos=[
                Molecule(
                            formula='DCN',
                            ppd_d_refs="Qi et al. 2008 ApJ 681, 1396",
                            ppd_d_bib_ids=["2008ApJ...681.1396Q"],
                ),
                Molecule(
                            formula='H13CN',
                            table_formula=r'H^{13}CN',
                            ppd_d_refs="Guzman et al. 2015 ApJ 814, 53",
                            ppd_d_bib_ids=["2015ApJ...814...53G"],
                ),
                Molecule(
                            formula='HC15N',
                            table_formula=r'H^{15}CN',
                            ppd_d_refs="Guzman et al. 2015 ApJ 814, 53",
                            ppd_d_bib_ids=["2015ApJ...814...53G"],
                ),
    ],
    ppd_isos_refs='[DCN] Qi et al. 2008 ApJ 681, 1396 [H13CN] Guzman et al. 2015 ApJ 814, 53 [HC15N] Guzman et al. 2015 ApJ 814, 53',
    exgal=True,
    exgal_d_refs="Rickard et al. 1977 ApJ 214, 390",
    exgal_d_bib_ids=["1977ApJ...214..390R"],
    exgal_sources="NGC 253, M82",
    exo=True,
    exo_d_refs="Hawker et al. 2018 ApJL 863, L11",
    exo_d_bib_ids=["2018ApJ...863L..11H"],
    Bcon=44316,
    mua=3.0,
    census_version='2018.0.0',
)
OCS = Molecule(
    name="carbonyl sulfide",
    formula="OCS",
    year=1971,
    label="OCS",
    astromol_name="OCS",
    sources=[SgrB2],
    telescopes=[NRAO36],
    wavelengths=["mm"],
    d_refs="Jefferts et al. 1971 ApJ 168, L111",
    l_refs="King & Gordy 1954 Phys Rev 93, 407",
    notes=None,
    ice=True,
    ice_d_refs="Palumbo et al. 1995 ApJ 449, 674; Palumbo et al. 1997 ApJ 479, 839",
    ice_l_refs="Palumbo et al. 1995 ApJ 449, 674",
    ice_d_bib_ids=["1995ApJ...449..674P","1997ApJ...479..839P"],
    exgal=True,
    exgal_d_refs="Mauersberger et al. 1995 A&A 294, 23",
    exgal_d_bib_ids=["1995A&A...294...23M"],
    exgal_sources="NGC 253",
    Bcon=6081,
    mua=0.7,
    census_version='2018.0.0',
)
HNC = Molecule(
    name="hydrogen isocyanide",
    formula="HNC",
    year=1972,
    label="HNC",
    astromol_name="HNC",
    sources=[W51, NGC2264],
    telescopes=[NRAO36],
    wavelengths=["mm"],
    d_refs="Snyder & Buhl 1972 Annals of the New York Academy of Science 194, 17; Zuckerman et al. 1972 ApJ 173, L125",
    l_refs="Blackman et al. 1976 Nature 261, 395",
    notes=None,
    ppd=True,
    ppd_d_refs="Dutrey et al. 1997 A&A 317, L55",
    ppd_d_bib_ids=["1997A&A...317L..55D"],
    ppd_isos=[
                Molecule(
                    formula = 'DNC',
                    ppd_d_refs="Loomis et al. 2020 ApJ 893, 101",
                    ppd_d_bib_ids=["2020ApJ...893..101L"],
                ),
    ],    
    exgal=True,
    exgal_d_refs="Henkel et al. 1988 A&A 201, L23",
    exgal_d_bib_ids=["1988A&A...201L..23H"],
    exgal_sources="IC 342",
    Bcon=45332,
    mua=3.1,
    census_version='2018.0.0',
)
H2S = Molecule(
    name="hydrogen sulfide",
    formula="H2S",
    year=1972,
    label="H2S",
    astromol_name="H2S",
    sources=[W3, W3OH, Orion, NGC2264, SgrB2, W51, DR21OH, NGC7538],
    telescopes=[NRAO36],
    wavelengths=["mm"],
    d_refs="Thaddeus et al. 1972 ApJ 176, L73",
    l_refs="Cupp et al. 1968 Phys Rev 171, 60",
    notes=None,
    exgal=True,
    exgal_d_refs="Hekkila et al. 1999 A&A 344, 817",
    exgal_d_bib_ids=["1999A&A...344..817H"],
    exgal_sources="LMC",
    ppd=True,
    ppd_d_refs="Phuong et al. 2018 A&A 616, L5",
    ppd_d_bib_ids=["2018A&A...616L...5P"],
    Acon=310584,
    Bcon=270368,
    Ccon=141820,
    mub=1.0,
    census_version='2018.0.0',
)
N2Hp = Molecule(
    name="protonated nitrogen",
    formula="N2H+",
    year=1974,
    label="N2H+",
    astromol_name="N2Hp",
    sources=[SgrB2, DR21, NGC6334, NGC2264],
    telescopes=[NRAO36],
    wavelengths=["mm"],
    d_refs="Turner 1974 ApJ 193, L83; Green et al. 1974 ApJ 193, L89; Thaddues & Turner 1975 ApJ 201, L25",
    l_refs="Saykally et al. 1976 ApJ 205, L101",
    notes=None,
    ppd=True,
    ppd_d_refs="Qi et al. 2003 ApJ 597, 986; Dutrey et al. 2007 A&A 464, 615",
    ppd_d_bib_ids=["2003ApJ...597..986Q", "2007A&A...464..615D"],
    ppd_isos=[
                Molecule(
                            formula="N2D+",
                            ppd_d_refs="Huang et al. 2015 ApJL 809, L26",
                            ppd_d_bib_ids=["2015ApJ...809L..26H"],
                ),
    ],
    ppd_isos_refs='[N2D+] Huang et al. 2015 ApJL 809, L26',
    exgal=True,
    exgal_d_refs="Mauersberger & Henkel 1991 A&A 245, 457",
    exgal_d_bib_ids=["1991A&A...245..457M"],
    exgal_sources="NGC 253, Maffei 2, IC 342, M82, NGC 6946",
    Bcon=46587,
    mua=3.4,
    census_version='2018.0.0',
)
C2H = Molecule(
    name="ethynyl radical",
    formula="C2H",
    year=1974,
    label="C2H",
    astromol_name="C2H",
    sources=[Orion],
    telescopes=[NRAO36],
    wavelengths=["mm"],
    d_refs="Tucker et al. 1974 ApJ 193, L115",
    l_refs="Sastry et al. 1981 ApJ 251, L119",
    notes=None,
    ppd=True,
    ppd_d_refs="Dutrey et al. 1997 A&A 317, L55",
    ppd_d_bib_ids=["1997A&A...317L..55D"],
    ppd_isos=[
                Molecule(
                    formula = 'C2D',
                    ppd_d_refs="Loomis et al. 2020 ApJ 893, 101",
                    ppd_d_bib_ids=["2020ApJ...893..101L"],
                ),
    ],    
    exgal=True,
    exgal_d_refs="Henkel et al. 1988 A&A 201, L23",
    exgal_d_bib_ids=["1988A&A...201L..23H"],
    exgal_sources="M82",
    Bcon=43675,
    mua=0.8,
    census_version='2018.0.0',
)
SO2 = Molecule(
    name="sulfur dioxide",
    formula="SO2",
    year=1975,
    label="SO2",
    astromol_name="SO2",
    sources=[Orion, SgrB2],
    telescopes=[NRAO36],
    wavelengths=["mm"],
    d_refs="Snyder et al. 1975 ApJ 198, L81",
    l_refs="Steenbeckeliers 1968 Ann. Soc. Sci. Brux 82, 331",
    notes=None,
    exgal=True,
    exgal_d_refs="Martin et al. 2003 A&A 411, L465",
    exgal_d_bib_ids=["2003A&A...411L.465M"],
    exgal_sources="NGC 253",
    Acon=60779,
    Bcon=10318,
    Ccon=8800,
    mub=1.6,
    census_version='2018.0.0',
)
HCO = Molecule(
    name="formyl radical",
    formula="HCO",
    year=1976,
    label="HCO",
    astromol_name="HCO",
    sources=[W3, NGC2024, W51, K350],
    telescopes=[NRAO36],
    wavelengths=["mm"],
    d_refs="Snyder et al. 1976 ApJ 208, L91",
    l_refs="Saito 1972 ApJ 178, L95",
    notes=None,
    exgal=True,
    exgal_d_refs="Sage & Ziurys 1995 ApJ 447, 625; Garcia-Burillo et al. 2002 ApJ 575, L55",
    exgal_d_bib_ids=["1995ApJ...447..625S", "2002ApJ...575L..55G"],
    exgal_sources="M82",
    Acon=7829365,
    Bcon=44788,
    Ccon=41930,
    mua=1.4,
    mub=0.7,
    census_version='2018.0.0',
)
HNO = Molecule(
    name="nitroxyl radical",
    formula="HNO",
    year=1977,
    label="HNO",
    astromol_name="HNO",
    sources=[SgrB2, NGC2024],
    telescopes=[NRAO36],
    wavelengths=["mm"],
    d_refs="Ulich et al. 1977 ApJ 217, L105",
    l_refs="Saito & Takagi 1973 JMS 47, 99",
    notes=None,
    Acon=553899,
    Bcon=42313,
    Ccon=39165,
    mua=1.0,
    mub=1.3,
    census_version='2018.0.0',
)
HCSp = Molecule(
    name="protonated carbon monosulfide",
    formula="HCS+",
    year=1981,
    label="HCS+",
    astromol_name="HCSp",
    sources=[Orion, SgrB2],
    telescopes=[NRAO36, Bell7m],
    wavelengths=["mm"],
    d_refs="Thaddeus et al. 1981 ApJ 246, L41",
    l_refs="Gudeman et al. 1981 ApJ 246, L47",
    notes=None,
    exgal=True,
    exgal_d_refs="Muller et al. 2013 A&A 551, A109",
    exgal_d_bib_ids=["2013A&A...551A.109M"],
    exgal_sources="PKS 1830-211 LOS",
    Bcon=10691,
    mua=1.9,
    census_version='2018.0.0',
)
HOCp = Molecule(
    name="hydroxymethyliumylidene",
    formula="HOC+",
    year=1983,
    label="HOC+",
    astromol_name="HOCp",
    sources=[SgrB2],
    telescopes=[FCRAO14m, Onsala20m],
    wavelengths=["mm"],
    d_refs="Woods et al. 1983 ApJ 270, 583",
    l_refs="Gudeman et al. 1982 PRL 48, 1344",
    notes="*Confirmed in 1995 ApJ 455, L73",
    exgal=True,
    exgal_d_refs="Usero et al. 2004 A&A 419, 897",
    exgal_d_bib_ids=["2004A&A...419..897U"],
    exgal_sources="NGC 1068",
    Bcon=44744,
    mua=4.0,
    census_version='2018.0.0',
)
SiC2 = Molecule(
    name="silacyclopropynylidene",
    formula="SiC2",
    year=1984,
    label="SiC2",
    astromol_name="SiC2",
    sources=[IRC10216],
    telescopes=[NRAO36, Bell7m],
    wavelengths=["mm"],
    cyclic=True,
    d_refs="Thaddeus et al. 1984 ApJ 283, L45",
    l_refs="Michalopoulos et al. 1984 JCP 80, 3556",
    notes=None,
    Acon=52474,
    Bcon=13157,
    Ccon=10443,
    mua=2.4,
    census_version='2018.0.0',
)
C2S = Molecule(
    name="dicarbon sulfide",
    formula="C2S",
    year=1987,
    label="C2S",
    astromol_name="C2S",
    sources=[TMC1, IRC10216, SgrB2],
    telescopes=[Nobeyama45, IRAM30],
    wavelengths=["cm", "mm"],
    d_refs="Saito et al. 1987 ApJ 317, L115",
    l_refs="Saito et al. 1987 ApJ 317, L115",
    notes="*Also Cernicharo et al. 1987 A&A 181, L9",
    ppd = True,
    ppd_d_refs="Phuong et al. 2021 A&A 653, L5",
    ppd_d_bib_ids = ["2021A&A...653L...5P"],
    ppd_sources = ["GG Tau"], #not yet implemented widely
    exgal=True,
    exgal_d_refs="Martin et al. 2006 ApJS 164, 450",
    exgal_d_bib_ids=["2006ApJS..164..450M"],
    exgal_sources="NGC 253",    
    Bcon=6478,
    mua=2.9,
    census_version='2021.1.0',
    change_log = {
    				'2018.0.0' : 'Initial Entry',
    				'2021.1.0' : 'Detection in protoplanetary disk'
    }
)
C3 = Molecule(
    name="tricarbon",
    formula="C3",
    year=1988,
    label="C3",
    astromol_name="C3",
    sources=[IRC10216],
    telescopes=[KPNO4m],
    wavelengths=["IR"],
    d_refs="Hinkle et al. 1988 Science 241, 1319",
    l_refs="Gausset et al. 1965 ApJ 142, 45",
    exgal=True,
    exgal_d_refs="Welty et al. 2013 MNRAS 428, 1107",
    exgal_d_bib_ids=["2013MNRAS.428.1107W"],
    exgal_sources="SMC",
    isotopologues="13CCC, C13CC",
    isos_d_refs="https://arxiv.org/abs/1911.09751",
    notes=None,
    mua=0.0,
    census_version='2018.0.0',
)
CO2 = Molecule(
    name="carbon dioxide",
    formula="CO2",
    year=1989,
    label="CO2",
    astromol_name="CO2",
    sources=[AFGL961LOS, AFGL989LOS, AFGL890LOS],
    telescopes=[IRAS],
    wavelengths=["IR"],
    d_refs="d'Hendecourt & Jourdain de Muizon 1989 A&A 223, L5; van Dishoeck et al. 1996 A&A 315, L349",
    l_refs="d'Hendecourt & Allamandola 1986 A&A Sup. Ser. 64, 453; Paso et al. 1980 JMS 79, 236; Reichle & Young 1972 Can J Phys 50, 2662",
    notes="*First detected in ices, then in gas phase",
    ice=True,
    ice_d_refs="d'Hendecourt & Jourdain de Muizon 1989 A&A 223, L5",
    ice_l_refs="d'Hendecourt & Allamandola 1986 A&A Sup. Ser. 64, 453",
    ice_d_bib_ids=["1989A&A...223L...5D"],
    ppd=True,
    ppd_d_refs="Carr & Najita 2008 Science 319, 1504",
    ppd_d_bib_ids=["2008Sci...319.1504C"],
    exo=True,
    exo_d_refs="Stevenson et al. 2010 Nature 464, 1161; Madhusudhan et al. 2011 Nature 469, 64; Lanotte et al. 2014 A&A 572, A73",
    exo_d_bib_ids=["2010Natur.464.1161S","2011Natur.469...64M","2014A&A...572A..73L"],
    mua=0.0,
    census_version='2018.0.0',
)
CH2 = Molecule(
    name="methylene",
    formula="CH2",
    year=1989,
    label="CH2",
    astromol_name="CH2",
    sources=[Orion],
    telescopes=[NRAOARO12],
    wavelengths=["mm"],
    d_refs="Hollis et al. 1989 ApJ 346, 794",
    l_refs="Lovas et al. 1983 ApJ 267, L131",
    notes="*Confirmed in 1995 ApJ 438, 259",
    Acon=2211494,
    Bcon=253618,
    Ccon=215102,
    mub=0.6,
    census_version='2018.0.0',
)
C2O = Molecule(
    name="dicarbon monoxide",
    formula="C2O",
    year=1991,
    label="C2O",
    astromol_name="C2O",
    sources=[TMC1],
    telescopes=[Nobeyama45],
    wavelengths=["cm"],
    d_refs="Ohishi et al. 1991 ApJ 380, L39",
    l_refs="Yamada et al. 1985 ApJ 290, L65",
    notes=None,
    Bcon=11546,
    mua=1.3,
    census_version='2018.0.0',
)
MgNC = Molecule(
    name="magnesium isocyanide",
    formula="MgNC",
    year=1993,
    label="MgNC",
    astromol_name="MgNC",
    sources=[IRC10216],
    telescopes=[IRAM30],
    wavelengths=["mm"],
    d_refs="Guélin et al. 1986 A&A 157, L17",
    l_refs="Kawaguchi et al. 1993 ApJ 406, L39",
    notes="*Actually identified in Kawaguchi et al. 1993 ApJ 406, L39 and Guélin et al. 1993 A&A 280, L19",
    Bcon=5967,
    mua=5.2,
    census_version='2018.0.0',
)
NH2 = Molecule(
    name="amidogen",
    formula="NH2",
    year=1993,
    label="NH2",
    astromol_name="NH2",
    sources=[SgrB2LOS],
    telescopes=[CSO],
    wavelengths=["sub-mm"],
    d_refs="van Dishoeck et al. 1993 ApJ 416, L83",
    l_refs="Charo et al. 1981 ApJ 244, L111",
    notes=None,
    exgal=True,
    exgal_d_refs="Muller et al. 2014 A&A 566, A112",
    exgal_d_bib_ids=["2014A&A...566A.112M"],
    exgal_sources="PKS 1830-211 LOS",
    Acon=710302,
    Bcon=388289,
    Ccon=245014,
    mub=1.8,
    isotopologues="NHD, ND2",
    isos_d_refs="Melosso et al. 2020 A&A 641, A153",
    isos_l_refs="Martin-Drumel et al. 2014 JPCA 118, 1331; Melosso et al. 2017 ApJS 233, 1; Bizzocchi et al. 2020 ApJS 247, 59",
    census_version='2018.0.0',
)

NaCN = Molecule(
    name="sodium cyanide",
    formula="NaCN",
    year=1994,
    label="NaCN",
    astromol_name="NaCN",
    sources=[IRC10216],
    telescopes=[NRAOARO12],
    wavelengths=["mm"],
    d_refs="Turner et al. 1994 ApJ 426, L97",
    l_refs="van Vaals et al. 1984 Chem Phys 86, 147",
    notes=None,
    Acon=57922,
    Bcon=8368,
    Ccon=7272,
    mua=8.9,
    census_version='2018.0.0',
)
N2O = Molecule(
    name="nitrous oxide",
    formula="N2O",
    year=1994,
    label="N2O",
    astromol_name="N2O",
    sources=[SgrB2],
    telescopes=[NRAOARO12],
    wavelengths=["mm"],
    d_refs="Ziurys et al. 1994 ApJ 436, L181",
    l_refs="Lovas 1978 J Phys Chem Ref Data 7, 1445",
    notes=None,
    Bcon=12562,
    mua=0.2,
    census_version='2018.0.0',
)
MgCN = Molecule(
    name="magnesium cyanide",
    formula="MgCN",
    year=1995,
    label="MgCN",
    astromol_name="MgCN",
    sources=[IRC10216],
    telescopes=[NRAOARO12, IRAM30],
    wavelengths=["mm"],
    d_refs="Ziurys et al. 1995 ApJ 445, L47",
    l_refs="Anderson et al. 1994 ApJ 429, L41",
    notes=None,
    Bcon=5095,
    mua="*",
    census_version='2018.0.0',
)
H3p = Molecule(
    name="",
    formula="H3+",
    year=1996,
    label="H3+",
    astromol_name="H3p",
    sources=[GL2136LOS, W33LOS],
    telescopes=[UKIRT],
    wavelengths=["IR"],
    d_refs="Geballe & Oka 1996 Nature 384, 334",
    l_refs="Oka 1980 PRL 45, 531",
    notes=None,
    exgal=True,
    exgal_d_refs="Geballe et al. 2006 ApJ 644, 907",
    exgal_d_bib_ids=["2006ApJ...644..907G"],
    exgal_sources="IRAS 08572+3915",
    mua=0.0,
    census_version='2018.0.0',
)
SiCN = Molecule(
    name="silicon monocyanide radical",
    formula="SiCN",
    year=2000,
    label="SiCN",
    astromol_name="SiCN",
    sources=[IRC10216],
    telescopes=[IRAM30],
    wavelengths=["mm"],
    d_refs="Guélin et al. 2000 A&A 363, L9",
    l_refs="Apponi et al. 2000 ApJ 536, L55",
    notes=None,
    Bcon=5543,
    mua=2.9,
    census_version='2018.0.0',
)
AlNC = Molecule(
    name="aluminum isocyanide",
    formula="AlNC",
    year=2002,
    label="AlNC",
    astromol_name="AlNC",
    sources=[IRC10216],
    telescopes=[IRAM30],
    wavelengths=["mm"],
    d_refs="Ziurys et al. 2002 ApJ 564, L45",
    l_refs="Robinson et al. 1997 Chem Phys Lett 278, 1",
    notes=None,
    Bcon=5985,
    mua=3.1,
    census_version='2018.0.0',
)
SiNC = Molecule(
    name="silicon monoisocyanide",
    formula="SiNC",
    year=2004,
    label="SiNC",
    astromol_name="SiNC",
    sources=[IRC10216],
    telescopes=[IRAM30],
    wavelengths=["mm"],
    d_refs="Guélin et al. 2004 A&A 426, L49",
    l_refs="Apponi et al. 2000 ApJ 536, L55",
    notes=None,
    Bcon=6397,
    mua=2.0,
    census_version='2018.0.0',
)
HCP = Molecule(
    name="phosphaethyne",
    formula="HCP",
    year=2007,
    label="HCP",
    astromol_name="HCP",
    sources=[IRC10216],
    telescopes=[IRAM30],
    wavelengths=["mm"],
    d_refs="Agúndez et al. 2007 ApJ 662, L91",
    l_refs="Bizzocchi et al. 2001 JMS 205, 110",
    notes="*First attempt 1990 ApJ 365, 59. Confirmed 2008 ApJ 684, 618",
    Bcon=19976,
    mua=0.4,
    census_version='2018.0.0',
)
CCP = Molecule(
    name="dicarbon phosphide radical",
    formula="CCP",
    year=2008,
    label="CCP",
    astromol_name="CCP",
    sources=[IRC10216],
    telescopes=[NRAOARO12],
    wavelengths=["mm"],
    d_refs="Halfen et al. 2008 ApJ 677, L101",
    l_refs="Halfen et al. 2008 ApJ 677, L101",
    notes=None,
    Bcon=6373,
    mua=3.4,
    census_version='2018.0.0',
)
AlOH = Molecule(
    name="aluminum hydroxide",
    formula="AlOH",
    year=2010,
    label="AlOH",
    astromol_name="AlOH",
    sources=[VYCaMaj],
    telescopes=[NRAOARO12, SMT10],
    wavelengths=["mm"],
    d_refs="Tenenbaum & Ziurys 2010 ApJ 712, L93",
    l_refs="Apponi et al. 1993 ApJ 414, L129",
    notes=None,
    Bcon=15740,
    mua=1.0,
    census_version='2018.0.0',
)
H2Op = Molecule(
    name="oxidaniumyl",
    formula="H2O+",
    year=2010,
    label="H2O+",
    astromol_name="H2Op",
    sources=[SgrB2, SgrB2LOS, NGC6334, DR21],
    telescopes=[Herschel],
    wavelengths=["sub-mm"],
    d_refs="Ossenkopf et al. 2010 A&A 518, L111; Gerin et al. 2010 A&A 518, L110",
    l_refs="Strahan et al. 1986 JCP 85, 1252; Murtz et al. 1998 JCP 109, 9744 ",
    notes=None,
    exgal=True,
    exgal_d_refs="Weiss et al. 2010 A&A 521, L1",
    exgal_d_bib_ids=["2010A&A...521L...1W"],
    exgal_sources="M82",
    Acon=870579,
    Bcon=372365,
    Ccon=253878,
    mub=2.4,
    census_version='2018.0.0',
)
H2Clp = Molecule(
    name="chloronium",
    formula="H2Cl+",
    year=2010,
    label="H2Cl+",
    astromol_name="H2Clp",
    sources=[SgrB2, SgrB2LOS, NGC6334, NGC6334LOS],
    telescopes=[Herschel],
    wavelengths=["sub-mm"],
    d_refs="Lis et al. 2010 A&A 521, L9",
    l_refs="Araki et al. 2001 JMS 210, 132",
    notes=None,
    exgal=True,
    exgal_d_refs="Muller et al. 2014 A&A 566, L6",
    exgal_d_bib_ids=["2014A&A...566L...6M"],
    exgal_sources="PKS 1830-211 LOS",
    Acon=337352,
    Bcon=273588,
    Ccon=148101,
    mub=1.9,
    census_version='2018.0.0',
)
KCN = Molecule(
    name="potassium cyanide",
    formula="KCN",
    year=2010,
    label="KCN",
    astromol_name="KCN",
    sources=[IRC10216],
    telescopes=[NRAOARO12, SMT10, IRAM30],
    wavelengths=["mm"],
    d_refs="Pulliam et al. 2010 ApJ 727, L181",
    l_refs="Torring et al. 1980 JCP 73, 4875",
    notes=None,
    Acon=58266,
    Bcon=4940,
    Ccon=4536,
    mub=10.0,
    census_version='2018.0.0',
)
FeCN = Molecule(
    name="iron cyanide",
    formula="FeCN",
    year=2011,
    label="FeCN",
    astromol_name="FeCN",
    sources=[IRC10216],
    telescopes=[NRAOARO12],
    wavelengths=["mm"],
    d_refs="Zack et al. 2011 ApJ 733, L36",
    l_refs="Flory & Ziurys 2011 JCP 135, 184303",
    notes=None,
    Bcon=4080,
    mua=4.5,
    census_version='2018.0.0',
)
HO2 = Molecule(
    name="hydroperoxyl radical",
    formula="HO2",
    year=2012,
    label="HO2",
    astromol_name="HO2",
    sources=[rhoOphA],
    telescopes=[IRAM30, APEX],
    wavelengths=["mm"],
    d_refs="Parise et al. 2012 A&A 541, L11",
    l_refs="Beers & Howard 1975 JCP 63, 4212; Saito 1977 JMS 65, 229; Charo & de Lucia 1982 JMS 94, 426",
    notes=None,
    Acon=610273,
    Bcon=33514,
    Ccon=31672,
    mua=1.4,
    mub=1.5,
    census_version='2018.0.0',
)
TiO2 = Molecule(
    name="titanium dioxide",
    formula="TiO2",
    year=2013,
    label="TiO2",
    astromol_name="TiO2",
    sources=[VYCaMaj],
    telescopes=[SMA, PdBI],
    wavelengths=["mm"],
    d_refs="Kamiński et al. 2013 A&A 551, A113",
    l_refs="Brunken 2008 APJ 676, 1367; Kania et al. 2011 JMS 268, 173",
    notes=None,
    Acon=30521,
    Bcon=8472,
    Ccon=6614,
    mua=6.3,
    census_version='2018.0.0',
)
CCN = Molecule(
    name="cyanomethylidyne",
    formula="CCN",
    year=2014,
    label="CCN",
    astromol_name="CCN",
    sources=[IRC10216],
    telescopes=[NRAOARO12, SMT10],
    wavelengths=["mm"],
    d_refs="Anderson & Ziurys 2014 ApJ 795, L1",
    l_refs="Anderson et al. 2015 JMS 307, 1",
    notes=None,
    Bcon=11939,
    mua=0.4,
    census_version='2018.0.0',
)
SiCSi = Molecule(
    name="disilicon carbide",
    formula="SiCSi",
    year=2015,
    label="SiCSi",
    astromol_name="SiCSi",
    sources=[IRC10216],
    telescopes=[IRAM30],
    wavelengths=["mm"],
    d_refs="Cernicharo et al. 2015 ApJ 806, L3",
    l_refs="McCarthy 2015 JPC Lett 6, 2107",
    notes=None,
    Acon=64074,
    Bcon=4396,
    Ccon=4102,
    mub=0.9,
    census_version='2018.0.0',
)
S2H = Molecule(
    name="hydrogen disulfide",
    formula="S2H",
    year=2017,
    label="S2H",
    astromol_name="S2H",
    sources=[HorseheadPDR],
    telescopes=[IRAM30],
    wavelengths=["mm"],
    d_refs="Fuente et al. 2017 ApJ 851, L49",
    l_refs="Tanimoto et al. 2000 JMS 199, 73",
    notes=None,
    Acon=296979,
    Bcon=7996,
    Ccon=7777,
    mua=1.2,
    mub=0.9,
    census_version='2018.0.0',
)
HCS = Molecule(
    name="thioformyl",
    formula="HCS",
    year=2018,
    label="HCS",
    astromol_name="HCS",
    sources=[L483],
    telescopes=[IRAM30],
    wavelengths=["mm"],
    d_refs="Agúndez et al. 2018 A&A 611, L1",
    l_refs="Habara et al. 2002 JCP 116, 9232",
    notes=None,
    Acon=954000,
    Bcon=20359,
    Ccon=19970,
    mua=0.4,
    mub=0.9,
    census_version='2018.0.0',
)
HSC = Molecule(
    name="sulfhydryl carbide",
    formula="HSC",
    year=2018,
    label="HSC",
    astromol_name="HSC",
    sources=[L483],
    telescopes=[IRAM30],
    wavelengths=["mm"],
    d_refs="Agúndez et al. 2018 A&A 611, L1",
    l_refs="Habara 2000 JCP 112, 10905",
    notes=None,
    Acon=295039,
    Bcon=22036,
    Ccon=19564,
    mua=2.5,
    mub=1.0,
    census_version='2018.0.0',
)
NCO = Molecule(
    name="isocyanate radical",
    formula="NCO",
    year=2018,
    label="NCO",
    astromol_name="NCO",
    sources=[L483],
    telescopes=[IRAM30],
    wavelengths=["mm"],
    d_refs="Marcelino et al. 2018 A&A 612, L10",
    l_refs="Kawaguchi et al. 1985 Mol Phys 55, 341; Saito and Amano 1970 JMS34, 383",
    notes=None,
    Bcon=11677,
    mua=0.6,
    census_version='2018.0.0',
)
CaNC = Molecule(
    name="calcium isocyanide",
    formula="CaNC",
    year=2019,
    label="CaNC",
    astromol_name="CaNC",
    sources=[IRC10216],
    telescopes=[IRAM30],
    wavelengths=["mm"],
    d_refs="Cernicharo et al. 2019 A&A 627, L4",
    l_refs="Steimle et al. 1993 ApJ 410, L49; Scurlock et al. 1994 JCP 100, 3497",
    notes="Dipole moment from Steimle et al. 1992 JCP 97, 2909",
    Bcon=4048,
    mua=6.985,
    census_version='2021.0.0',
)
NCS = Molecule(
    name="thiocyanogen",
    formula="NCS",
    year=2021,
    label="NCS",
    astromol_name="NCS",
    sources=[TMC1],
    telescopes=[Yebes40],
    wavelengths=["cm"],
    d_refs="Cernicharo et al. 2021 A&A 648, L3",
    l_refs="Amano & Amano 1991 J Chem Phys 95, 2275; McCarthy et al. 2003 ApJS 144, 287; Maeda et al. 2007 Mol Phys 105, 477",
    notes="Dipole moment taken from CDMS as calculated by Holger Müller.  Amano & Amano point out that this calculation shoudl be quite challenging, and neither of the more recent laboratory papers appears to perform it, so this value should be viewed with caution.",
    Bcon=6106,
    mua=2.45,
    census_version='2021.0.0',
)
MgC2 = Molecule(
    name="magnesium dicarbide",
    formula="MgC2",
    year=2022,
    label="MgC2",
    astromol_name="MgC2",
    sources=[IRC10216],
    telescopes=[Yebes40,IRAM30],
    wavelengths=["cm","mm"],
    d_refs="Changala et al. 2022 ApJL 940, L42",
    d_ref_bib_ids = ["Changala:2022:L42"],
    l_refs="Changala et al. 2022 ApJL 940, L42",
    l_ref_bib_ids = ["Changala:2022:L42"],
    notes = "Dipole moment from Itono et al. 2000 ApJ 538, L163.  Detection is of 24MgC2 (main isotopologue) as well as 25MgC2 and 26MgC2.",
    Acon=51900,
    Bcon=11504,
    Ccon=9392,
    mua=7.9,
    census_version='2021.6.0',
    change_log = {
    				'2021.6.0' : 'Initial entry',
    			}
)

######################################################################
#                           Four Atoms                               #
######################################################################


NH3 = Molecule(
    name="ammonia",
    formula="NH3",
    year=1968,
    label="NH3",
    astromol_name="NH3",
    sources=[GalacticCenter],
    telescopes=[HatCreek],
    wavelengths=["cm"],
    d_refs="Cheung et al. 1968 PRL 25, 1701",
    l_refs="Cleeton & Williams 1934 Phys Rev 45, 234",
    notes=None,
    ice=True,
    ice_d_refs="Lacy et al. 1998 ApJ 501, L105",
    ice_l_refs="d'Hendecourt & Allamandola 1986 A&A Sup. Ser. 64, 453",
    ice_d_bib_ids=["1998ApJ...501L.105L"],
    ppd=True,
    ppd_d_refs="Salinas et al. 2016 A&A 591, A122",
    ppd_d_bib_ids=["2016A&A...591A.122S"],
    exgal=True,
    exgal_d_refs="Martin & Ho 1979 A&A 74, L7",
    exgal_d_bib_ids=["1979A&A....74L...7M"],
    exgal_sources="IC 342, NGC 253",
    exo=True,
    exo_d_refs="Giacobbe et al. 2021 Nature 592, 205",
    exo_d_bib_ids=["2021Natur.592..205G"],
    exo_sources="HD 209458b",
    Acon=298193,
    Bcon=298193,
    Ccon=286696,
    muc=1.5,
    census_version='2018.0.0',
)
H2CO = Molecule(
    name="formaldehyde",
    formula="H2CO",
    year=1969,
    label="H2CO",
    astromol_name="H2CO",
    sources=[
        M17LOS,
        M3LOS,
        W49LOS,
        NGC2024LOS,
        DR21LOS,
        W43LOS,
        W44LOS,
        W51LOS,
        SgrALOS,
        SgrB2LOS,
        W33LOS,
        NGC6334LOS,
        CasALOS,
    ],
    telescopes=[NRAO140],
    wavelengths=["cm"],
    d_refs="Snyder et al. 1969 PRL 22, 679",
    l_refs="Shinegari 1967 J Phys Soc Jpn 23, 404",
    notes=None,
    ice=True,
    ice_d_refs="Keane et al. 2001 A&A 376, 254",
    ice_l_refs="Schutte et al. 1993 Icarus 104, 118",
    ice_d_bib_ids=["2001A&A...376..254K"],
    ppd=True,
    ppd_d_refs="Dutrey et al. 1997 A&A 317, L55",
    ppd_d_bib_ids=["1997A&A...317L..55D"],
    exgal=True,
    exgal_d_refs="Gardner & Whiteoak 1974 Nature 247, 526",
    exgal_d_bib_ids=["1974Natur.247..526G"],
    exgal_sources="NGC 253, NGC 4945",
    Acon=281971,
    Bcon=38834,
    Ccon=34004,
    mua=2.3,
    census_version='2018.0.0',
)
HNCO = Molecule(
    name="isocyanic acid",
    formula="HNCO",
    year=1972,
    label="HNCO",
    astromol_name="HNCO",
    sources=[SgrB2],
    telescopes=[NRAO36],
    wavelengths=["cm", "mm"],
    d_refs="Snyder & Buhl 1972 ApJ 177, 619",
    l_refs="Kewley et al. 1963 JMS 10, 418",
    notes=None,
    exgal=True,
    exgal_d_refs="Nguyen-Q-Rieu et al. 1991 A&A 241, L33",
    exgal_d_bib_ids=["1991A&A...241L..33N"],
    exgal_sources="NGC 253, Maffei 2, IC 342",
    Acon=912711,
    Bcon=11071,
    Ccon=10911,
    mua=1.6,
    mub=1.4,
    census_version='2018.0.0',
)
H2CS = Molecule(
    name="thioformaldehyde",
    formula="H2CS",
    year=1973,
    label="H2CS",
    astromol_name="H2CS",
    sources=[SgrB2LOS],
    telescopes=[Parkes64],
    wavelengths=["cm"],
    d_refs="Sinclair et al. 1973 Aust. J. Phys. 26, 85",
    l_refs="Johnson & Powell 1970 Science 169, 679",
    notes=None,
    ppd=True,
    ppd_d_refs="Le Gal et al. 2019 ApJ 876, 72; Loomis et al. 2020 ApJ 893, 101",
    ppd_d_bib_ids=["2020ApJ...893..101L", "2019ApJ...876...72L"],    
    exgal=True,
    exgal_d_refs="Martin et al. 2006 ApJS 164, 450",
    exgal_d_bib_ids=["2006ApJS..164..450M"],
    exgal_sources="NGC 253",
    Acon=291292,
    Bcon=17700,
    Ccon=16652,
    mua=1.6,
    census_version='2018.0.0',
)
C2H2 = Molecule(
    name="acetylene",
    formula="C2H2",
    year=1976,
    label="C2H2",
    astromol_name="C2H2",
    sources=[IRC10216],
    telescopes=[KPNO4m],
    wavelengths=["IR"],
    d_refs="Ridgway et al. 1976 Nature 264, 345",
    l_refs="Baldacci et al. 1973 JMS 48, 600",
    notes=None,
    ppd=True,
    ppd_d_refs="Lahuis et al. 2006 ApJ 636, L145",
    ppd_d_bib_ids=["2006ApJ...636L.145L"],
    exgal=True,
    exgal_d_refs="Matsuura et al. 2002 ApJ 580, L133",
    exgal_d_bib_ids=["2002ApJ...580L.133M"],
    exgal_sources="LMC",
    exo=True,
    exo_d_refs="Giacobbe et al. 2021 Nature 592, 205",
    exo_d_bib_ids=["2021Natur.592..205G"],
    exo_sources="HD 209458b",
    mua=0.0,
    census_version='2018.0.0',
)
C3N = Molecule(
    name="cyanoethynyl radical",
    formula="C3N",
    year=1977,
    label="C3N",
    astromol_name="C3N",
    sources=[IRC10216, TMC1],
    telescopes=[NRAO36, Onsala20m],
    wavelengths=["cm", "mm"],
    d_refs="Guelin & Thaddeus 1977 ApJ 212, L81",
    l_refs="Gottlieb et al. 1983 ApJ 275, 916",
    notes="*Confirmed in Friberg et al. 1980 ApJ 241, L99",
    exgal=True,
    exgal_d_refs="Tercero et al. 2020 A&AL 636, L7",
    exgal_d_bib_ids=["2020A&A...636L...7T"],
    exgal_sources="PKS 1830-211 LOS",
    Bcon=4968,
    mua=2.9,
    census_version='2018.0.0',
)
HNCS = Molecule(
    name="isothiocyanic acid",
    formula="HNCS",
    year=1979,
    label="HNCS",
    astromol_name="HNCS",
    sources=[SgrB2],
    telescopes=[Bell7m, NRAO36],
    wavelengths=["mm"],
    d_refs="Frerking et al. 1979 ApJ 234, L143",
    l_refs="Kewley et al. 1963 JMS 10, 418",
    notes=None,
    Acon=1348662,
    Bcon=5883,
    Ccon=5847,
    mua=1.6,
    mub="*",
    census_version='2018.0.0',
)
HOCOp = Molecule(
    name="protonated carbon dioxide",
    formula="HOCO+",
    year=1981,
    label="HOCO+",
    astromol_name="HOCOp",
    sources=[SgrB2],
    telescopes=[Bell7m],
    wavelengths=["mm"],
    d_refs="Thaddeus et al. 1981 ApJ 246, L41",
    l_refs="Green et al. 1976 Chem Phys 17, 479; Bogey et al. 1984 A&A 138, L11",
    notes=None,
    exgal=True,
    exgal_d_refs="Aladro et al. 2015 A&A 579, A101; Martin et al. 2006 ApJS 164, 450",
    exgal_d_bib_ids=["2015A&A...579A.101A", "2006ApJS..164..450M"],
    exgal_sources="NGC 253",
    Acon=789951,
    Bcon=10774,
    Ccon=10609,
    mua=2.7,
    mub=1.8,
    census_version='2018.0.0',
)
C3O = Molecule(
    name="tricarbon monoxide",
    formula="C3O",
    year=1985,
    label="C3O",
    astromol_name="C3O",
    sources=[TMC1],
    telescopes=[NRAOARO12, FCRAO14m, Nobeyama45],
    wavelengths=["cm", "mm"],
    d_refs="Matthews et al. 1984 Nature 310, 125",
    l_refs="Brown et al. 1983 JACS 105, 6496",
    notes="*Confirmed in Brown et al. 1985 ApJ 297, 302 and Kaifu et al. 2004 PASJ 56, 69",
    Bcon=4811,
    mua=2.4,
    census_version='2018.0.0',
)
lC3H = Molecule(
    name="propynylidyne radical",
    formula="C3H",
    table_formula="l-C3H",
    year=1985,
    label="l-C3H",
    astromol_name="lC3H",
    sources=[TMC1, IRC10216],
    telescopes=[NRAO36, Bell7m, FCRAO14m, Onsala20m],
    wavelengths=["cm", "mm"],
    d_refs="Thaddeus et al. 1985 ApJ 294, L49",
    l_refs="Gottlieb et al. 1985 ApJ 294, L55",
    notes=None,
    exgal=True,
    exgal_d_refs="Muller et al. 2011 A&A 535, A103",
    exgal_d_bib_ids=["2011A&A...535A.103M"],
    exgal_sources="PKS 1830-211 LOS",
    Bcon=11189,
    mua=3.6,
    mub=0.5,
    census_version='2018.0.0',
)
HCNHp = Molecule(
    name="protonated hydrogen cyanide",
    formula="HCNH+",
    year=1986,
    label="HCNH+",
    astromol_name="HCNHp",
    sources=[SgrB2],
    telescopes=[NRAOARO12, MWO4m],
    wavelengths=["mm"],
    d_refs="Ziurys & Turner 1986 ApJ 302, L31",
    l_refs="Bogey et al. 1985 JCP 83, 3703; Altman et al. 1984 JCP 80, 3911",
    notes=None,
    Bcon=37056,
    mua=0.3,
    census_version='2018.0.0',
)
H3Op = Molecule(
    name="hydronium",
    formula="H3O+",
    year=1986,
    label="H3O+",
    astromol_name="H3Op",
    sources=[Orion, SgrB2],
    telescopes=[NRAOARO12, MWO4m],
    wavelengths=["mm"],
    d_refs="Wootten et al. 1986 A&A 166, L15; Hollis et al. 1986 Nature 322, 524",
    l_refs="Plummer et al. 1985 JCP 83, 1428; Bogey et al. 1985 A&A 148, L11; Liu & Oka 1985 PRL 54, 1787",
    notes="*Confirmed in Wootten et al. 1991 ApJ 390, L79",
    exgal=True,
    exgal_d_refs="van der Tak et al. 2008 A&A 477, L5",
    exgal_d_bib_ids=["2008A&A...477L...5V"],
    exgal_sources="M82, Arp 220",
    Acon=334405,
    Bcon=334405,
    Ccon=184725,
    muc=1.4,
    census_version='2018.0.0',
)
C3S = Molecule(
    name="tricarbon monosulfide",
    formula="C3S",
    year=1987,
    label="C3S",
    astromol_name="C3S",
    sources=[TMC1, IRC10216],
    telescopes=[Nobeyama45, IRAM30],
    wavelengths=["cm", "mm"],
    d_refs="Yamamoto et al. 1987 ApJ 317, L119",
    l_refs="Yamamoto et al. 1987 ApJ 317, L119",
    notes=None,
    Bcon=2890,
    mua=3.7,
    census_version='2018.0.0',
)
cC3H = Molecule(
    name="cyclopropenylidene radical",
    formula="C3H",
    table_formula="c-C3H",
    year=1987,
    label="c-C3H",
    astromol_name="cC3H",
    sources=[TMC1],
    telescopes=[Nobeyama45],
    wavelengths=["mm"],
    cyclic=True,
    d_refs="Yamamoto et al. 1987 ApJ 322, L55",
    l_refs="Yamamoto et al. 1987 ApJ 322, L55",
    notes=None,
    exgal="Tentative",
    exgal_d_refs="Martin et al. 2006 ApJS 164, 450",
    exgal_d_bib_ids=["2006ApJS..164..450M"],
    exgal_sources="NGC 253",
    Acon=44517,
    Bcon=34016,
    Ccon=19189,
    mua=2.4,
    census_version='2018.0.0',
)
HC2N = Molecule(
    name="cyanocarbene radical",
    formula="HC2N",
    year=1991,
    label="HC2N",
    astromol_name="HC2N",
    sources=[IRC10216],
    telescopes=[IRAM30],
    wavelengths=["mm"],
    d_refs="Guélin & Cernicharo 1991 A&A 244, L21",
    l_refs="Saito et al. 1984 JCP 80, 1427; Brown et al. 1990 JMS 143, 203",
    notes=None,
    Bcon=10986,
    mua=3.0,
    census_version='2018.0.0',
)
H2CN = Molecule(
    name="methylene amidogen radical",
    formula="H2CN",
    year=1994,
    label="H2CN",
    astromol_name="H2CN",
    sources=[TMC1],
    telescopes=[NRAOARO12],
    wavelengths=["cm"],
    d_refs="Ohishi et al. 1994 ApJ 427, L51",
    l_refs="Yamamoto & Saito 1992 JCP 96, 4157",
    notes=None,
    exgal=True,
    exgal_d_refs="Tercero et al. 2020 A&AL 636, L7",
    exgal_d_bib_ids=["2020A&A...636L...7T"],
    exgal_sources="PKS 1830-211 LOS",
    Acon=284343,
    Bcon=39158,
    Ccon=34246,
    mua=2.5,
    census_version='2018.0.0',
)
SiC3 = Molecule(
    name="silicon tricarbide",
    formula="SiC3",
    year=1999,
    label="SiC3",
    astromol_name="SiC3",
    sources=[IRC10216],
    telescopes=[NRAOARO12],
    wavelengths=["mm"],
    cyclic=True,
    d_refs="Apponi et al. 1999 ApJ 516, L103",
    l_refs="Apponi et al. 1999 JCP 111, 3911; McCarthy et al. JCP 110, 1064",
    notes=None,
    Acon=37944,
    Bcon=6283,
    Ccon=5387,
    mua=4.0,
    census_version='2018.0.0',
)
CH3 = Molecule(
    name="methyl radical",
    formula="CH3",
    year=2000,
    label="CH3",
    astromol_name="CH3",
    sources=[SgrALOS],
    telescopes=[ISO],
    wavelengths=["IR"],
    d_refs="Feuchtgruber et al. 2000 ApJ 535, L111",
    l_refs="Yamada et al. 1981 JCP 75, 5256",
    notes=None,
    mua=0.0,
    census_version='2018.0.0',
)
C3Nm = Molecule(
    name="cyanoethynyl anion",
    formula="C3N-",
    year=2008,
    label="C3N-",
    astromol_name="C3Nm",
    sources=[IRC10216],
    telescopes=[IRAM30],
    wavelengths=["mm"],
    d_refs="Thaddeus et al. 2008 ApJ 677, 1132",
    l_refs="Thaddeus et al. 2008 ApJ 677, 1132",
    notes=None,
    Bcon=4852,
    mua=3.1,
    census_version='2018.0.0',
)
PH3 = Molecule(
    name="phosphine",
    formula="PH3",
    year=2008,
    label="PH3",
    astromol_name="PH3",
    sources=[IRC10216, CRL2688],
    telescopes=[IRAM30, Herschel, SMT10, CSO],
    wavelengths=["mm", "sub-mm"],
    d_refs="Agúndez et al. 2008 A&A 485, L33",
    l_refs="Cazzoli & Puzzarini 2006 JMS 239, 64; Sousa-Silva et al. 2013 JMS 288, 28; Muller 2013 JQSRT 130, 335",
    notes="*Confirmed in Agúndez et al. 2014 ApJL 790, L27",
    Acon=133480,
    Bcon=133480,
    Ccon=117488,
    census_version='2018.0.0',
)
HCNO = Molecule(
    name="fulminic acid",
    formula="HCNO",
    year=2009,
    label="HCNO",
    astromol_name="HCNO",
    sources=[B1b, L1544, L183, L1527],
    telescopes=[IRAM30],
    wavelengths=["mm"],
    d_refs="Marcelino et al. 2009 ApJ 690, L27",
    l_refs="Winnewisser & Winnewisser 1971 Z Naturforsch 26, 128",
    notes=None,
    Bcon=11469,
    mua=3.1,
    census_version='2018.0.0',
)
HOCN = Molecule(
    name="cyanic acid",
    formula="HOCN",
    year=2009,
    label="HOCN",
    astromol_name="HOCN",
    sources=[SgrB2],
    telescopes=[Bell7m, NRAO36, NRAOARO12],
    wavelengths=["mm"],
    d_refs="Brünken et al. 2009 ApJ 697, 880",
    l_refs="Brünken et al. 2009 ApJ 697, 880",
    notes="*Confirmed in Brünken et al. 2010 A&A 516, A109",
    Acon=681000,
    Bcon=10577,
    Ccon=10398,
    mua=3.7,
    census_version='2018.0.0',
)
HSCN = Molecule(
    name="thiocyanic acid",
    formula="HSCN",
    year=2009,
    label="HSCN",
    astromol_name="HSCN",
    sources=[SgrB2],
    telescopes=[NRAOARO12],
    wavelengths=["mm"],
    d_refs="Halfen et al. 2009 ApJ 702, L124",
    l_refs="Brunken et al. 2009 ApJ 706, 1588",
    notes=None,
    Acon=289830,
    Bcon=5795,
    Ccon=5675,
    mua=3.3,
    census_version='2018.0.0',
)
HOOH = Molecule(
    name="hydrogen peroxide",
    formula="HOOH",
    year=2011,
    label="HOOH",
    astromol_name="HOOH",
    sources=[rhoOphA],
    telescopes=[APEX],
    wavelengths=["mm"],
    d_refs="Bergman et al. 2011 A&A 531, L8",
    l_refs="Petkie et al. 1995 JMS 171, 145; Helminger et al. 1981 JMS 85, 120",
    notes=None,
    Acon=301878,
    Bcon=26212,
    Ccon=25099,
    muc=1.6,
    census_version='2018.0.0',
)
lC3Hp = Molecule(
    name="propynylidyne cation",
    formula="C3H+",
    table_formula="l-C3H+",
    year=2012,
    label="l-C3H+",
    astromol_name="lC3Hp",
    sources=[HorseheadPDR],
    telescopes=[IRAM30],
    wavelengths=["mm"],
    d_refs="Pety et al. 2012 A&A 549, A68",
    l_refs="Brunken et al. 2014 ApJ 783, L4",
    notes=None,
    exgal=True,
    exgal_d_refs="Tercero et al. 2020 A&AL 636, L7",
    exgal_d_bib_ids=["2020A&A...636L...7T"],
    exgal_sources="PKS 1830-211 LOS",
    Bcon=11245,
    mua=3.0,
    census_version='2021.4.0',
    change_log = {
    				'2018.0.0' : 'Initial entry',
    				'2021.4.0' : 'Corrected error in name (was duplicating cC3Hp)',
    			}    
)
HMgNC = Molecule(
    name="hydromagnesium isocyanide",
    formula="HMgNC",
    year=2013,
    label="HMgNC",
    astromol_name="HMgNC",
    sources=[IRC10216],
    telescopes=[IRAM30],
    wavelengths=["mm"],
    d_refs="Cabezas et al. 2013 ApJ 75, 133",
    l_refs="Cabezas et al. 2013 ApJ 75, 133",
    notes=None,
    Bcon=5481,
    mua=3.5,
    census_version='2018.0.0',
)
HCCO = Molecule(
    name="ketenyl radical",
    formula="HCCO",
    year=2015,
    label="HCCO",
    astromol_name="HCCO",
    sources=[Lupus1A, L483],
    telescopes=[IRAM30],
    wavelengths=["mm"],
    d_refs="Agúndez et al. 2015 A&A 577, L5",
    l_refs="Endo & Hirota 1987 JCP 86, 4319; Oshima & Endo 1993 JMS 159, 458",
    notes=None,
    Bcon=10831,
    mua=1.6,
    census_version='2018.0.0',
)
CNCN = Molecule(
    name="isocyanogen",
    formula="CNCN",
    year=2018,
    label="CNCN",
    astromol_name="CNCN",
    sources=[L483],
    telescopes=[IRAM30],
    wavelengths=["mm"],
    d_refs="Agundez et al. 2018 ApJL 861, L22",
    l_refs="Gerry et al. 1990 JMS 140, 147; Winnewisser et al. 1992 JMS 153, 635",
    notes=None,
    Bcon=5174,
    mua=0.7,
    census_version='2018.0.0',
)
HONO = Molecule(
    name="nitrous acid",
    formula="HONO",
    year=2019,
    label="HONO",
    astromol_name="HONO",
    sources=[IRAS16293],
    telescopes=[ALMA],
    wavelengths=["sub-mm"],
    d_refs="Coutens et al. 2019 A&A 623, L13",
    l_refs="Guilmot et al. 1993 JMS 160, 387; Guilmot et al. 1993 JMS 160, 401; Dehayem-Kamadjeu et al. 2005 JMS 234, 182",
    notes="Only lines of trans-HONO are claimed as detected. As such, constants for this entry are for trans-HONO.",
    Acon=92892,
    Bcon=12525,
    Ccon=11017,
    mua=1.378,
    mub=1.242,
    census_version='2021.0.0',
)
MgCCH = Molecule(
    name="magnesium ethynyl radical",
    formula="MgCCH",
    year=2019,
    label="MgCCH",
    astromol_name="MgCCH",
    sources=[IRC10216],
    telescopes=[IRAM30],
    wavelengths=["mm"],
    d_refs="Agundez et al. 2014 A&A 570, A45 (tentative); Cernicharo et al. 2019 A&A 630, L2 (confirmation)",
    l_refs="Brewster et al. 1999 Chem. Phys. Lett. 310, 411",
    notes="Dipole from Woon 1996 ApJ 456, 602",
    Bcon=4965,
    mua=1.68,
    census_version='2121.0.0',
)
HCCS = Molecule(
    name="thioketenyl radical",
    formula="HCCS",
    year=2021,
    label="HCCS",
    astromol_name="HCCS",
    sources=[TMC1],
    telescopes=[Yebes40],
    wavelengths=["cm"],
    d_refs="Cernicharo et al. 2021 A&A 648, L3",
    l_refs="Kim et al. 2002 J Mol Spec 212, 83; Vrtilek et al. 1992 ApJ 398, L73",
    notes="Dipole moment is from Vrtilek et al. who in turn got it from J. D. Goddard 1992, private communication.",
    Bcon=5876,
    mua=1.2,
    census_version='2121.0.0',
)
HNCN = Molecule(
    name="cyanomidyl radical",
    formula="HNCN",
    year=2021,
    label="HNCN",
    astromol_name="HNCN",
    sources=[G0693],
    telescopes=[Yebes40, IRAM30],
    wavelengths=["cm","mm"],
    d_refs="Rivilla et al. 2021 MNRAS 506, L79",
    d_ref_bib_ids = ["2021MNRAS.506L..79R"],
    l_refs="Yamamoto & Saito 1994 J Chem Phys 101, 10350",
    l_ref_bib_ids = ["1994JChPh.10110350Y"],
    notes="Dipole moment from CDMS.",
    Acon=634900,
    Bcon=11087,
    Ccon=10882,
    mua=2.28,
    mub=1.,
    census_version='2021.1.0',
)

H2NC = Molecule(
    name="aminocarbyne",
    formula="H2NC",
    year=2021,
    label="H2NC",
    astromol_name="H2NC",
    sources=[L483],
    telescopes=[IRAM30],
    wavelengths=["mm"],
    d_refs="Cabezas et al. 2021 A&A 654, A45",
    d_ref_bib_ids = ["2021A&A...654A..45C"],
    l_refs="Cabezas et al. 2021 A&A 654, A45",
    l_ref_bib_ids = ["2021A&A...654A..45C"],
    notes="No laboratory work; purely theoretical constants.",
    exgal=True,
    exgal_d_refs="Cabezas et al. 2021 A&A 654, A45",
    exgal_d_bib_ids=["2021A&A...654A..45C"],
    exgal_sources="PKS 1830-211 LOS",
    Acon=339670,
    Bcon=38086,
    Ccon=34124,
    mua=3.83,
    census_version='2021.3.0',
    change_log = {
    				'2021.3.0' : 'Initial Entry',
    }
)
HCCSp = Molecule(
    name="thioketenylium",
    formula="HCCS+",
    table_formula="HCCS+",
    year=2022,
    label="HCCS+",
    astromol_name="HCCSp",
    sources=[TMC1],
    telescopes=[IRAM30,Yebes40],
    wavelengths=["cm","mm"],
    d_refs="Cabezas et al. 2022 A&A 657, L4",
    d_ref_bib_ids = ["2022A&A...657L...4C"],
    l_refs="Cabezas et al. 2022 A&A 657, L4",
    l_ref_bib_ids = ["2022A&A...657L...4C"],
    notes = "Identification performed based on quantum chemical calculations only.  Cabezas et al. do not state what dipole moment they assumed, but it was presumably taken from Puzzarini 2008 Chem Phys 346, 45",
    Bcon=6021,
    mua=2.3,
    census_version='2021.4.0',
    change_log = {
    				'2021.4.0' : 'Initial entry',
    			}
)

######################################################################
#                           Five Atoms                               #
######################################################################


HC3N = Molecule(
    name="cyanoacetylene",
    formula="HC3N",
    year=1971,
    label="HC3N",
    astromol_name="HC3N",
    sources=[SgrB2],
    telescopes=[NRAO140],
    wavelengths=["cm"],
    d_refs="Turner 1971 ApJ 163, L35",
    l_refs="Tyler & Sheridan 1963 Trans Faraday Soc 59, 2661",
    notes="*Confirmed in Dickinson 1972 AL 12, 235",
    ppd=True,
    ppd_d_refs="Chapillon et al. 2012 ApJ 756, 58",
    ppd_d_bib_ids=["2012ApJ...756...58C"],
    exgal=True,
    exgal_d_refs="Mauersberger et al. 1990 A&A 236, 63; Henkel et al. 1988 A&A 201, L23",
    exgal_d_bib_ids=["1990A&A...236...63M", "1988A&A...201L..23H"],
    exgal_sources="NGC 253",
    Bcon=4549,
    mua=3.7,
    census_version='2018.0.0',
)
HCOOH = Molecule(
    name="formic acid",
    formula="HCOOH",
    year=1971,
    label="HCOOH",
    astromol_name="HCOOH",
    sources=[SgrB2],
    telescopes=[NRAO140],
    wavelengths=["cm"],
    d_refs="Zukerman et al. 1971 ApJ 163, L41",
    l_refs="Zukerman et al. 1971 ApJ 163, L41; Bellet et al. 1971 J Mol Struct 9, 49; Bellet et al. 1971 J Mol Struct 9, 65",
    notes="*Confirmed in Winnewisser & Churchwell 1975 ApJ 200, L33",
    ice=True,
    ice_d_refs="Schutte et al. 1999 A&A 343, 966",
    ice_l_refs="Schutte et al. 1999 A&A 343, 966",
    ice_d_bib_ids=["1999A&A...343..966S"],
    ppd=True,
    ppd_d_refs="Favre et al. 2018 ApJL 862, L2",
    ppd_d_bib_ids=["2018ApJ...862L...2F"],
    exgal=True,
    exgal_d_refs="Tercero et al. 2020 A&AL 636, L7",
    exgal_d_bib_ids=["2020A&A...636L...7T"],
    exgal_sources="PKS 1830-211 LOS",
    Acon=77512,
    Bcon=12055,
    Ccon=10416,
    mua=1.4,
    mub=0.2,
    census_version='2018.0.0',
)
CH2NH = Molecule(
    name="methanimine",
    formula="CH2NH",
    year=1973,
    label="CH2NH",
    astromol_name="CH2NH",
    sources=[SgrB2],
    telescopes=[Parkes64],
    wavelengths=["cm"],
    d_refs="Godfrey et al. 1973 ApL 13, 119",
    l_refs="Godfrey et al. 1973 ApL 13, 119; Johnson & Lovas 1972 CPL 15, 65",
    notes=None,
    exgal=True,
    exgal_d_refs="Muller et al. 2011 A&A 535, A103",
    exgal_d_bib_ids=["2011A&A...535A.103M"],
    exgal_sources="PKS 1830-211 LOS",
    Acon=196211,
    Bcon=34532,
    Ccon=29352,
    mua=1.3,
    mub=1.5,
    census_version='2018.0.0',
)
NH2CN = Molecule(
    name="cyanamide",
    formula="NH2CN",
    year=1975,
    label="NH2CN",
    astromol_name="NH2CN",
    sources=[SgrB2],
    telescopes=[NRAO36],
    wavelengths=["mm"],
    d_refs="Turner et al. 1975 ApJ 201, L149",
    l_refs="Tyler et al. 1972 JMS 43, 248; Miller et al. 1962 JMS 8, 153; Lide 1962 JMS 8, 142; Johnson & Suenram 1976 ApJ 208, 245",
    notes=None,
    exgal=True,
    exgal_d_refs="Martin et al. 2006 ApJS 164, 450",
    exgal_d_bib_ids=["2006ApJS..164..450M"],
    exgal_sources="NGC 253",
    Acon=312142,
    Bcon=10130,
    Ccon=9866,
    mua=4.3,
    muc=1.0,
    census_version='2018.0.0',
)
H2CCO = Molecule(
    name="ketene",
    formula="H2CCO",
    year=1977,
    label="H2CCO",
    astromol_name="H2CCO",
    sources=[SgrB2],
    telescopes=[NRAO36],
    wavelengths=["mm"],
    d_refs="Turner 1977 ApJ 213, L75",
    l_refs="Johnson & Strandberg 1952 JCP 20, 687; Johns et al. 1972 JMS 42, 523",
    notes=None,
    exgal=True,
    exgal_d_refs="Muller et al. 2011 A&A 535, A103",
    exgal_d_bib_ids=["2011A&A...535A.103M"],
    exgal_sources="PKS 1830-211 LOS",
    Acon=282473,
    Bcon=10294,
    Ccon=9916,
    mua=1.4,
    census_version='2018.0.0',
)
C4H = Molecule(
    name="butadiynyl radical",
    formula="C4H",
    year=1978,
    label="C4H",
    astromol_name="C4H",
    sources=[IRC10216],
    telescopes=[NRAO36],
    wavelengths=["mm"],
    d_refs="Guélin et al. 1978 ApJ 224, L27",
    l_refs="Gottlieb et al. 1983 ApJ 275, 916",
    notes=None,
    exgal=True,
    exgal_d_refs="Muller et al. 2011 A&A 535, A103",
    exgal_d_bib_ids=["2011A&A...535A.103M"],
    exgal_sources="PKS 1830-211 LOS",
    Bcon=4759,
    mua=0.9,
    census_version='2018.0.0',
)
SiH4 = Molecule(
    name="silane",
    formula="SiH4",
    year=1984,
    label="SiH4",
    astromol_name="SiH4",
    sources=[IRC10216],
    telescopes=[IRTF],
    wavelengths=["IR"],
    d_refs="Goldhaber and Betz 1977 ApJ 279, L55",
    l_refs="Goldhaber and Betz 1977 ApJ 279, L55",
    notes=None,
    mua=0.0,
    census_version='2018.0.0',
)
cC3H2 = Molecule(
    name="cyclopropenylidene",
    formula="C3H2",
    table_formula="c-C3H2",
    year=1985,
    label="c-C3H2",
    astromol_name="cC3H2",
    sources=[SgrB2, Orion, TMC1],
    telescopes=[Bell7m],
    wavelengths=["cm", "mm"],
    cyclic=True,
    d_refs="Thaddeus et al. 1985 ApJ 299, L63",
    l_refs="Thaddeus et al. 1985 ApJ 299, L63",
    notes="*See also Vrtilek et al. 1987 ApJ 314, 716",
    ppd=True,
    ppd_d_refs="Qi et al. 2013 ApJL 765, L14",
    ppd_d_bib_ids=["2013ApJ...765L..14Q"],
    exgal=True,
    exgal_d_refs="Seaquist & Bell 1986 ApJ 303, L67",
    exgal_d_bib_ids=["1986ApJ...303L..67S"],
    exgal_sources="NGC 5128",
    Acon=35093,
    Bcon=32213,
    Ccon=16749,
    mub=3.4,
    census_version='2018.0.0',
)
CH2CN = Molecule(
    name="cyanomethyl radical",
    formula="CH2CN",
    year=1988,
    label="CH2CN",
    astromol_name="CH2CN",
    sources=[TMC1, SgrB2],
    telescopes=[FCRAO14m, NRAO140, Onsala20m, Nobeyama45],
    wavelengths=["cm"],
    d_refs="Irvine et al. 1988 ApJ 334, L107",
    l_refs="Saito et al. 1988 ApJ 334, L113",
    notes=None,
    exgal=True,
    exgal_d_refs="Muller et al. 2011 A&A 535, A103",
    exgal_d_bib_ids=["2011A&A...535A.103M"],
    exgal_sources="PKS 1830-211 LOS",
    isotopologues='HDCCN',
    isos_d_refs='[HDCCN] Cabezas et al. 2021 A&A 646, 1',
    isos_l_refs='[HDCCN] Cabezas et al. 2021 A&A 646, 1',
    Acon=285130,
    Bcon=10246,
    Ccon=9877,
    mua=1.6,
    census_version='2018.0.0',
)
C5 = Molecule(
    name="pentacarbon",
    formula="C5",
    year=1989,
    label="C5",
    astromol_name="C5",
    sources=[IRC10216],
    telescopes=[KPNO4m],
    wavelengths=["IR"],
    d_refs="Bernath et al. 1989 Science 244, 562",
    l_refs="Vala et al. 1989 JCP 90, 595",
    notes=None,
    mua=0.0,
    census_version='2018.0.0',
)
SiC4 = Molecule(
    name="silicon tetracarbide",
    formula="SiC4",
    year=1989,
    label="SiC4",
    astromol_name="SiC4",
    sources=[IRC10216],
    telescopes=[Nobeyama45],
    wavelengths=["cm", "mm"],
    d_refs="Ohishi et al. 1989 ApJ 345, L83",
    l_refs="Ohishi et al. 1989 ApJ 345, L83",
    notes=None,
    Bcon=1534,
    mua=6.4,
    census_version='2018.0.0',
)
H2CCC = Molecule(
    name="propadienylidene",
    formula="H2CCC",
    year=1991,
    label="H2CCC",
    astromol_name="H2CCC",
    sources=[TMC1],
    telescopes=[IRAM30, Effelsberg100],
    wavelengths=["cm", "mm"],
    d_refs="Cernicharo et al. 1991 ApJ 368, L39",
    l_refs="Vrtilek et al. 1990 ApJ 364, L53",
    notes=None,
    exgal=True,
    exgal_d_refs="Muller et al. 2011 A&A 535, A103",
    exgal_d_bib_ids=["2011A&A...535A.103M"],
    exgal_sources="PKS 1830-211 LOS",
    Acon=288775,
    Bcon=10589,
    Ccon=10204,
    mua=4.1,
    census_version='2018.0.0',
)
CH4 = Molecule(
    name="methane",
    formula="CH4",
    year=1991,
    label="CH4",
    astromol_name="CH4",
    sources=[NGC7538LOS],
    telescopes=[IRTF],
    wavelengths=["IR"],
    d_refs="Lacy et al. 1991 ApJ 376, 556",
    l_refs="Champion et al. 1989 JMS 133, 256; d'Hendecourt & Allamandola 1986 A&A Supp Ser. 64, 453 ",
    notes=None,
    ice=True,
    ice_d_refs="Lacy et al. 1991 ApJ 376, 556",
    ice_l_refs="d'Hendecourt & Allamandola 1986 A&A Sup. Ser. 64, 453",
    ice_d_bib_ids=["1991ApJ...376..556L"],
    ppd=True,
    ppd_d_refs="Gibb et al. 2013 ApJL 776, L28",
    ppd_d_bib_ids=["2013ApJ...776L..28G"],
    exo=True,
    exo_d_refs="Swain et al. 2008 Nature 452, 329; Barman et al. 2011 ApJ 733, 65; Stevenson et al. 2014 ApJ 791, 36; Barman et al. 2015 ApJ 804, 61",
    exo_d_bib_ids=["2008Natur.452..329S","2011ApJ...733...65B","2014ApJ...791...36S","2015ApJ...804...61B"],
    mua=0.0,
    census_version='2018.0.0',
)
HCCNC = Molecule(
    name="isocyanoacetylene",
    formula="HCCNC",
    year=1992,
    label="HCCNC",
    astromol_name="HCCNC",
    sources=[TMC1],
    telescopes=[Nobeyama45],
    wavelengths=["cm", "mm"],
    d_refs="Kawaguchi et al. 1992 ApJ 386, L51",
    l_refs="Kruger et al. 2010 Ang. Chem. 23, 1644",
    notes=None,
    Bcon=4968,
    mua=2.9,
    census_version='2018.0.0',
)
HNCCC = Molecule(
    name="",
    formula="HNCCC",
    year=1992,
    label="HNCCC",
    astromol_name="HNCCC",
    sources=[TMC1],
    telescopes=[Nobeyama45],
    wavelengths=["cm"],
    d_refs="Kawaguchi et al. 1992 ApJ 396, L49",
    l_refs="Kawaguchi et al. 1992 ApJ 396, L49",
    notes=None,
    Bcon=4668,
    mua=5.7,
    census_version='2018.0.0',
)
H2COHp = Molecule(
    name="protonated formaldehyde",
    formula="H2COH+",
    year=1996,
    label="H2COH+",
    astromol_name="H2COHp",
    sources=[SgrB2, Orion, W51],
    telescopes=[Nobeyama45, NRAOARO12],
    wavelengths=["cm", "mm"],
    d_refs="Ohishi et al. 1996 ApJ 471, L61",
    l_refs="Chomiak et al. 1994 Can J Phys 72, 1078",
    notes=None,
    Acon=197582,
    Bcon=34351,
    Ccon=29173,
    mua=1.4,
    mub=1.8,
    census_version='2018.0.0',
)
C4Hm = Molecule(
    name="butadiynyl anion",
    formula="C4H-",
    year=2007,
    label="C4H-",
    astromol_name="C4Hm",
    sources=[IRC10216],
    telescopes=[IRAM30],
    wavelengths=["mm"],
    d_refs="Cernicharo et al. 2007 A&A 467, L37",
    l_refs="Gupta et al. 2007 ApJ 655, L57",
    notes=None,
    Bcon=4655,
    mua=6.2,
    census_version='2018.0.0',
)
CNCHO = Molecule(
    name="cyanoformaldehyde",
    formula="CNCHO",
    year=2007,
    label="CNCHO",
    astromol_name="CNCHO",
    sources=[SgrB2],
    telescopes=[GBT],
    wavelengths=["cm"],
    d_refs="Remijan et al. 2008 ApJ 675, L85",
    l_refs="Bogey et al. 1988 CPL 146, 227; Bogey et al. 1995 JMS 172, 344",
    notes=None,
    Acon=67470,
    Bcon=5010,
    Ccon=4657,
    mua=0.8,
    mub=1.9,
    census_version='2018.0.0',
)
HNCNH = Molecule(
    name="carbodiimide",
    formula="HNCNH",
    year=2012,
    label="HNCNH",
    astromol_name="HNCNH",
    sources=[SgrB2],
    telescopes=[GBT],
    wavelengths=["cm"],
    d_refs="McGuire et al. 2012 ApJ 758, L33",
    l_refs="Birk et al. 1989 JMS 135, 402; Wagener et al. 1995 JMS 170, 323; Jabs et al. 1997 Chem Phys 225, 77",
    notes=None,
    Acon=379244,
    Bcon=10367,
    Ccon=10366,
    mub=1.9,
    census_version='2018.0.0',
)
CH3O = Molecule(
    name="methoxy radical",
    formula="CH3O",
    year=2012,
    label="CH3O",
    astromol_name="CH3O",
    sources=[B1b],
    telescopes=[IRAM30],
    wavelengths=["mm"],
    d_refs="Cernicharo et al. 2012 ApJ 759, L43",
    l_refs="Momose et al. 1988 JCP 88, 5338; Endo et al. 1984 JCP 81, 122",
    notes=None,
    Acon=513887,
    Bcon=27930,
    Ccon=27930,
    mua=2.1,
    census_version='2018.0.0',
)
NH3Dp = Molecule(
    name="ammonium ion",
    formula="NH3D+",
    year=2013,
    label="NH3D+",
    astromol_name="NH3Dp",
    sources=[Orion, B1b],
    telescopes=[IRAM30],
    wavelengths=["mm"],
    d_refs="Gupta et al. 2013 ApJ 778, L1",
    l_refs="Gupta et al. 2013 ApJ 778, L1",
    notes="*Confirmed in Marcelino et al. 2018 A&A 612, L10",
    Acon=175439,
    Bcon=131412,
    Ccon=131412,
    mua=0.3,
    census_version='2018.0.0',
)
H2NCOp = Molecule(
    name="protonated isocyanic acid",
    formula="H2NCO+",
    year=2013,
    label="H2NCO+",
    astromol_name="H2NCOp",
    sources=[SgrB2, L483],
    telescopes=[GBT],
    wavelengths=["cm"],
    d_refs="Cernicharo et al. 2013 ApJ 771, L10",
    l_refs="Cernicharo et al. 2013 ApJ 771, L10",
    notes="*See also Doménech et al. 2013 ApJ 77, L11",
    Acon=319800,
    Bcon=10279,
    Ccon=9949,
    mua=4.1,
    census_version='2018.0.0',
)
NCCNHp = Molecule(
    name="protonated cyanogen",
    formula="NCCNH+",
    year=2015,
    label="NCCNH+",
    astromol_name="NCCNHp",
    sources=[TMC1, L483],
    telescopes=[IRAM30, Yebes40],
    wavelengths=["cm", "mm"],
    d_refs="Agúndez et al. 2015 A&A 579, L10",
    l_refs="Amano & Scappini 1991 JCP 95, 2280; Gottlieb et al. 200 JCP 113, 1910",
    notes=None,
    Bcon=4438,
    mua=6.5,
    census_version='2018.0.0',
)
CH3Cl = Molecule(
    name="chloromethane",
    formula="CH3Cl",
    year=2017,
    label="CH3Cl",
    astromol_name="CH3Cl",
    sources=[IRAS16293],
    telescopes=[ALMA],
    wavelengths=["mm"],
    d_refs="Fayolle et al. 2017 Nature Astron. 1, 702",
    l_refs="Wlodarczak et al. 1986 JMS 116, 251",
    notes=None,
    Acon=156051,
    Bcon=13293,
    Ccon=13293,
    mua=1.9,
    census_version='2018.0.0',
)
MgC3N = Molecule(
    name="magnesium cyanoethynyl radical",
    formula="MgC3N",
    year=2019,
    label="MgC3N",
    astromol_name="MgC3N",
    sources=[IRC10216],
    telescopes=[IRAM30, Yebes40],
    wavelengths=["cm", "mm"],
    d_refs="Cernicharo et al. 2019 A&A 630, L2",
    l_refs="Cernicharo et al. 2019 A&A 630, L2",
    notes="Assigned based entirely on quantum chemistry; no lab work.",
    Bcon=1381,
    mua=6.3,
    census_version='2021.0.0',
)
HC3Op = Molecule(
    name="protonated tricarbon monoxide",
    formula="HC3O+",
    year=2020,
    label="HC3O+",
    astromol_name="HC3Op",
    sources=[TMC1],
    telescopes=[IRAM30, Yebes40],
    wavelengths=["cm", "mm"],
    d_refs="Cernicharo et al. 2020 A&A 642, L17",
    l_refs="Cernicharo et al. 2020 A&A 642, L17",
    notes=None,
    Bcon=4461,
    mua=3.4,
    census_version='2021.0.0',
)
NH2OH = Molecule(
    name="hydroxylamine",
    formula="NH2OH",
    year=2020,
    label="NH2OH",
    astromol_name="NH2OH",
    sources=[G0693],
    telescopes=[IRAM30],
    wavelengths=['mm'],
    d_refs="Rivilla et al. 2020 ApJL 899, L28",
    l_refs="Tsunekawa 1972 J Phys Soc Japn 33, 167; Morino et al. 2000 J Mol Struct 517-518, 367",
    notes=None,
    Acon=190976,
    Bcon=25219,
    Ccon=25157,
    mua=0.589,
    muc=0.060,
    census_version="2021.0.0",
)
HC3Sp = Molecule(
    name="protonated tricarbon monosulfide",
    formula="HC3S+",
    year=2021,
    label="HC3S+",
    astromol_name="HC3Sp",
    sources=[TMC1],
    telescopes=[Yebes40],
    wavelengths=["cm"],
    d_refs="Cernicharo et al. 2021 A&A 646, L3",
    l_refs="Cernicharo et al. 2021 A&A 646, L3",
    notes=None,
    Bcon=2735,
    mua=1.7,
    census_version='2021.0.0',
)
H2CCS = Molecule(
    name="thioketene",
    formula="H2CCS",
    year=2021,
    label="H2CCS",
    astromol_name="H2CCS",
    sources=[TMC1],
    telescopes=[Yebes40],
    wavelengths=["cm"],
    d_refs="Cernicharo et al. 2021 A&A 648, L3",
    l_refs="Georgiou et al. 1979 J Mol Spectrosc 77, 365; Winnewisser & Schäfer 1980 Z Natur Forsch A 35, 483; McNaughton et al. 1996 J Mol Spectrosc 175, 377",
    notes="Dipole is from Georgiou et al. 1979.",
    Acon=286616,
    Bcon=5663,
    Ccon=5548,
    mua=1.02,
    census_version='2021.0.0',
)
C4S = Molecule(
    name="tetracarbon monosulfide",
    formula="C4S",
    year=2021,
    label="C4S",
    astromol_name="C4S",
    sources=[TMC1],
    telescopes=[Yebes40],
    wavelengths=["cm"],
    d_refs="Cernicharo et al. 2021 A&A 648, L3",
    l_refs="Hirahara et al. 1993 ApJL 408, L113;  Gordon et al. 2001 ApJS 134, 311",
    notes="Dipole moment from Lee 1997 Chem Phys Lett 268, 69 and Pascoli & Lavendy 1998 Int J Mass Spectr 181, 11",
    Bcon=1519,
    mua=4.03,
    census_version='2021.0.0',
)
CHOSH = Molecule(
    name="monothioformic acid",
    formula="CHOSH",
    year=2021,
    label="CHOSH",
    astromol_name="CHOSH",
    sources=[G0693],
    telescopes=[Yebes40, IRAM30],
    wavelengths=["cm", "mm"],
    d_refs="Rodríguez-Almeida et al. 2021 ApJL 912, L11",
    l_refs="Hocking & Winnewisser 1976 Z Naturforsch A 31, 995",
    notes="Dataset was refit by H.S.P. Müller for CDMS, and those are the predictions used for the detection.",
    Acon=62036,
    Bcon=6125,
    Ccon=5570,
    mua=1.366,
    mub=0.702,
    census_version='2021.0.0',
)
HCSCN = Molecule(
    name="monothioformic acid",
    formula="HCSCN",
    year=2021,
    label="HCSCN",
    astromol_name="HCSCN",
    sources=[TMC1],
    telescopes=[Yebes40],
    wavelengths=["cm"],
    d_refs="Cernicharo et al. 2021 A&A 650, L14",
    d_ref_bib_ids=["2021A&A...650L..14C"],
    l_refs="Cernicharo et al. 2021 A&A 650, L14",
    l_ref_bib_ids=["2021A&A...650L..14C"],
    Acon=42910,
    Bcon=3117,
    Ccon=2901,
    mua=2.16,
    mub=3.09,
    census_version='2021.1.0',
)
HC3O = Molecule(
    name="tricarbon monoxide",
    formula="HC3O",
    year=2021,
    label="HC3O",
    astromol_name="HC3O",
    sources=[TMC1],
    telescopes=[Yebes40],
    wavelengths=["cm"],
    d_refs="Cernicharo et al. 2021 A&A 656, L21",
    d_ref_bib_ids=["2021A&A...656L..21C"],
    l_refs="Cooksy et al. 1992 ApJ 386, L27; Cooksy et al. 1992 JMS 153, 610; Chen et al. 1996 ApJ 462, 561",
    l_ref_bib_ids=["1992ApJ...386L..27C", "1992JMoSp.153..610C", "1996ApJ...462..561C"],
    notes= "Dipole from Cooksy et al. 1995 J Phys Chem 99, 11095",
    Acon=261120,
    Bcon=4577,
    Ccon=4489,
    mua=2.39,
    census_version='2021.4.0',
    change_log = {
    				'2021.4.0' : 'Initial entry',
    			}
)


######################################################################
#                           Six Atoms                               #
######################################################################


CH3OH = Molecule(
    name="methanol",
    formula="CH3OH",
    year=1970,
    label="CH3OH",
    astromol_name="CH3OH",
    sources=[SgrA, SgrB2],
    telescopes=[NRAO140],
    wavelengths=["cm"],
    d_refs="Ball et al. 1970 ApJ 162, L203",
    l_refs="Ball et al. 1970 ApJ 162, L203",
    notes=None,
    ice=True,
    ice_d_refs="Grim et al. 1991 A&A 243, 473",
    ice_l_refs="d'Hendecourt & Allamandola 1986 A&A Sup. Ser. 64, 453",
    ice_d_bib_ids=["1991A&A...243..473G"],
    ppd=True,
    ppd_d_refs="Walsh et al. 2016 ApJL 823, L10",
    ppd_d_bib_ids=["2016ApJ...823L..10W"],
    exgal=True,
    exgal_d_refs="Henkel et al. 1987 A&A 188, L1",
    exgal_d_bib_ids=["1987A&A...188L...1H"],
    exgal_sources="NGC 253, IC 342",
    Acon=127523,
    Bcon=24690,
    Ccon=23760,
    mua=0.9,
    mub=1.4,
    census_version='2018.0.0',
)
CH3CN = Molecule(
    name="methyl cyanide",
    formula="CH3CN",
    year=1971,
    label="CH3CN",
    astromol_name="CH3CN",
    sources=[SgrA, SgrB2],
    telescopes=[NRAO36],
    wavelengths=["mm"],
    d_refs="Solomon et al. 1971 ApJ 168, L107",
    l_refs="Cord et al. 1968 Microwave Spectral Tables V5; Kessler et al. Phys Rev 79, 54",
    notes=None,
    ppd=True,
    ppd_d_refs="Oberg et al. 2015 Nature 520, 198",
    ppd_d_bib_ids=["2015Natur.520..198O"],
    exgal=True,
    exgal_d_refs="Mauersberger et al. 1991 A&A 247, 307",
    exgal_d_bib_ids=["1991A&A...247..307M"],
    exgal_sources="NGC 253",
    Acon=158099,
    Bcon=9199,
    Ccon=9199,
    mua=3.9,
    census_version='2018.0.0',
)
NH2CHO = Molecule(
    name="formamide",
    formula="NH2CHO",
    year=1971,
    label="NH2CHO",
    astromol_name="NH2CHO",
    sources=[SgrB2],
    telescopes=[NRAO140],
    wavelengths=["cm"],
    d_refs="Rubin et al. 1971 ApJ 169, L39",
    l_refs="Rubin et al. 1971 ApJ 169, L39",
    notes=None,
    exgal=True,
    exgal_d_refs="Muller et al. 2013 A&A 551, A109",
    exgal_d_bib_ids=["2013A&A...551A.109M"],
    exgal_sources="PKS 1830-211 LOS",
    Acon=72717,
    Bcon=11373,
    Ccon=9834,
    mua=3.6,
    mub=0.9,
    census_version='2018.0.0',
)
CH3SH = Molecule(
    name="methyl mercaptan",
    formula="CH3SH",
    year=1979,
    label="CH3SH",
    astromol_name="CH3SH",
    sources=[SgrB2],
    telescopes=[Bell7m],
    wavelengths=["mm"],
    d_refs="Linke et al. 1979 ApJ 234, L139",
    l_refs="Kilb 1955 JCP 23, 1736",
    notes=None,
    exgal=True,
    exgal_d_refs="Tercero et al. 2020 A&AL 636, L7",
    exgal_d_bib_ids=["2020A&A...636L...7T"],
    exgal_sources="PKS 1830-211 LOS",
    Acon=102771,
    Bcon=12952,
    Ccon=12400,
    mua=1.3,
    mub=0.8,
    census_version='2018.0.0',
)
C2H4 = Molecule(
    name="ethylene",
    formula="C2H4",
    year=1981,
    label="C2H4",
    astromol_name="C2H4",
    sources=[IRC10216],
    telescopes=[McMath],
    wavelengths=["IR"],
    d_refs="Betz 1981 ApJ 244, L103",
    l_refs="Lambeau et al. 1980 JMS 81, 227",
    notes=None,
    mua=0.0,
    census_version='2018.0.0',
)
C5H = Molecule(
    name="pentynylidyne radical",
    formula="C5H",
    year=1986,
    label="C5H",
    astromol_name="C5H",
    sources=[IRC10216],
    telescopes=[IRAM30],
    wavelengths=["cm"],
    d_refs="Cernicharo et al. 1986 A&A 164, L1",
    l_refs="Gottlieb et al. 1986 A&A 164, L5",
    notes="*See also Cernicharo et al. 1986 A&A 167, L5 and Cernicharo et al. 1987 A&A 172, L5",
    Bcon=2395,
    mua=4.9,
    census_version='2018.0.0',
)
CH3NC = Molecule(
    name="methyl isocyanide",
    formula="CH3NC",
    year=1988,
    label="CH3NC",
    astromol_name="CH3NC",
    sources=[SgrB2],
    telescopes=[IRAM30],
    wavelengths=["mm"],
    d_refs="Cernicharo et al. 1988 A&A 189, L1",
    l_refs="Kukolich 1972 JCP 57, 869; Ring et al. 1947 Phys Rev 72, 1262",
    notes="*Confirmed in Remijan et al. 2005 ApJ 632, 333 and Gratier et al. 2013 557, A101",
    Acon=157151,
    Bcon=10053,
    Ccon=10053,
    mua=3.9,
    census_version='2018.0.0',
)
HC2CHO = Molecule(
    name="propynal",
    formula="HC2CHO",
    year=1988,
    label="HC2CHO",
    astromol_name="HC2CHO",
    sources=[TMC1],
    telescopes=[NRAO140, Nobeyama45],
    wavelengths=["cm"],
    d_refs="Irvine et al. 1988 ApJ 335, L89",
    l_refs="Winnewisser 1973 JMS 46, 16",
    notes=None,
    Acon=68035,
    Bcon=4826,
    Ccon=4500,
    mua=2.4,
    mub=0.6,
    census_version='2018.0.0',
)
H2C4 = Molecule(
    name="butatrienylidene",
    formula="H2C4",
    year=1991,
    label="H2C4",
    astromol_name="H2C4",
    sources=[IRC10216],
    telescopes=[IRAM30],
    wavelengths=["mm"],
    d_refs="Cernicharo et al. 1991 ApJ 368, L43",
    l_refs="Killian et al. 1990 ApJ 365, L89",
    notes=None,
    Acon=286234,
    Bcon=4503,
    Ccon=4429,
    mua=4.1,
    census_version='2018.0.0',
)
C5S = Molecule(
    name="pentacarbon monosulfide radical",
    formula="C5S",
    year=1993,
    label="C5S",
    astromol_name="C5S",
    sources=[IRC10216],
    telescopes=[NRAO140],
    wavelengths=["cm"],
    d_refs="Bell et al. 1993 ApJ 417, L37",
    l_refs="Kasai et al. 1993 ApJ 410, L45; Gordon et al. 2001 ApJS 134, 311",
    notes="*Confirmed in Agúndez et al. 2014 A&A 570, A45",
    Bcon=923,
    mua=5.1,
    census_version='2018.0.0',
)
HC3NHp = Molecule(
    name="protonated cyanoacetylene",
    formula="HC3NH+",
    year=1994,
    label="HC3NH+",
    astromol_name="HC3NHp",
    sources=[TMC1],
    telescopes=[Nobeyama45],
    wavelengths=["cm"],
    d_refs="Kawaguchi et al. 1994 ApJ 420, L95",
    l_refs="Lee & Amano 1987 ApJ 323",
    notes=None,
    Bcon=4329,
    mua=1.6,
    census_version='2018.0.0',
)
C5N = Molecule(
    name="cyanobutadiynyl radical",
    formula="C5N",
    year=1998,
    label="C5N",
    astromol_name="C5N",
    sources=[TMC1],
    telescopes=[IRAM30, Effelsberg100],
    wavelengths=["cm"],
    d_refs="Guélin et al. 1998 A&A 355, L1",
    l_refs="Kasai et al. 1997 ApJ 477, L65",
    notes=None,
    Bcon=1403,
    mua=3.4,
    census_version='2018.0.0',
)
HC4H = Molecule(
    name="diacetylene",
    formula="HC4H",
    year=2001,
    label="HC4H",
    astromol_name="HC4H",
    sources=[CRL618],
    telescopes=[ISO],
    wavelengths=["IR"],
    d_refs="Cernicharo et al. 2001 ApJ 546, L123",
    l_refs="Arie & Johns 1992 JMS 155, 195",
    notes="*Confirmed in 2018 ApJ 852, 80",
    exgal=True,
    exgal_d_refs="Bernard-Salas et al. 2006 ApJ 652, L29",
    exgal_d_bib_ids=["2006ApJ...652L..29B"],
    exgal_sources="SMP LMC 11",
    mua=0.0,
    census_version='2018.0.0',
)
HC4N = Molecule(
    name="",
    formula="HC4N",
    year=2004,
    label="HC4N",
    astromol_name="HC4N",
    sources=[IRC10216],
    telescopes=[IRAM30],
    wavelengths=["mm"],
    d_refs="Cernicharo et al. 2004 ApJ 615, L145",
    l_refs="Tang et al. 1999 CPL 315, 69",
    notes=None,
    Bcon=2302,
    mua=4.3,
    census_version='2018.0.0',
)
cH2C3O = Molecule(
    name="cyclopropenone",
    formula="H2C3O",
    table_formula="c-H2C3O",
    year=2006,
    label="c-H2C3O",
    astromol_name="cH2C3O",
    sources=[SgrB2],
    telescopes=[GBT],
    wavelengths=["cm"],
    cyclic=True,
    d_refs="Hollis et al. 2006 ApJ 642, 933",
    l_refs="Benson et al. 1973 JACS 95, 2772; Guillemin et al. 1990 JMS 140, 190",
    notes=None,
    Acon=32041,
    Bcon=7825,
    Ccon=6281,
    mua=4.4,
    census_version='2018.0.0',
)
CH2CNH = Molecule(
    name="ketenimine",
    formula="CH2CNH",
    year=2006,
    label="CH2CNH",
    astromol_name="CH2CNH",
    sources=[SgrB2],
    telescopes=[GBT],
    wavelengths=["cm"],
    d_refs="Lovas et al. 2006 ApJ 645, L137",
    l_refs="Rodler et al. 1984 CPL 110, 447; Rodler et al. 1986 JMS 118, 267",
    notes=None,
    Acon=201444,
    Bcon=9663,
    Ccon=9470,
    mua=0.4,
    mub=1.4,
    census_version='2018.0.0',
)
C5Nm = Molecule(
    name="cyanobutadiynyl anion",
    formula="C5N-",
    year=2008,
    label="C5N-",
    astromol_name="C5Nm",
    sources=[IRC10216],
    telescopes=[IRAM30],
    wavelengths=["mm"],
    d_refs="Cernicharo et al. 2008 ApJ 688, L83",
    l_refs="Botschwina & Oswald 2008 JCP 129, 044305",
    notes=None,
    Bcon=1389,
    mua=5.2,
    census_version='2018.0.0',
)
HNCHCN = Molecule(
    name="E-cyanomethanimine",
    formula="HNCHCN",
    year=2013,
    label="HNCHCN",
    astromol_name="HNCHCN",
    sources=[SgrB2],
    telescopes=[GBT],
    wavelengths=["cm"],
    d_refs="Zaleski et al. 2013 ApJ 765, L9",
    l_refs="Zaleski et al. 2013 ApJ 765, L9",
    notes=None,
    Bcon=1389,
    mua=3.3,
    mub=2.5,
    census_version='2018.0.0',
)
SiH3CN = Molecule(
    name="silyl cyanide",
    formula="SiH3CN",
    year=2014,
    label="SiH3CN",
    astromol_name="SiH3CN",
    sources=[IRC10216],
    telescopes=[IRAM30],
    wavelengths=["mm"],
    d_refs="Agúndez et al. 2014 A&A 570, A45",
    l_refs="Priem et al. 1998 JMS 191, 183",
    notes="*Confirmed in Cernicharo et al. 2017 A&A 606, L5",
    Acon=62695,
    Bcon=4972,
    Ccon=4600,
    mua=3.4,
    census_version='2018.0.0',
)
MgC4H = Molecule(
    name="magnesium butadiynyl radical",
    formula="MgC4H",
    year=2019,
    label="MgC4H",
    astromol_name="MgC4H",
    sources=[IRC10216],
    telescopes=[IRAM30, Yebes40],
    wavelengths=["cm", "mm"],
    d_refs="Cernicharo et al. 2019 A&A 630, L2",
    l_refs="Forthomme et al. 2010 Chem. Phys. Lett. 488, 116",
    notes="Lab spectroscopy is electronic - no pure rotational spectra are available for this species.  Assignment was made based on quantum chemical calculations performed in Cernicharo et al. 2019.",
    Bcon=1381,
    mua=2.1,
    census_version='2021.0.0',
)
CH3COp = Molecule(
    name="acetyl cation",
    formula="CH3CO+",
    year=2021,
    label="CH3CO+",
    astromol_name="CH3COp",
    sources=[TMC1, L483, L1527, L1544],
    telescopes=[IRAM30, Yebes40],
    wavelengths=["cm", "mm"],
    d_refs="Cernicharo et al. 2021 A&A 646, L7",
    l_refs="Cernicharo et al. 2021 A&A 646, L7",
    notes="",
    Bcon=9134,
    mua=3.5,
    census_version='2021.0.0',
)
H2CCCS = Molecule(
    name="propadienthione",
    formula="H2CCCS",
    year=2021,
    label="H2CCCS",
    astromol_name="H2CCCS",
    sources=[TMC1],
    telescopes=[Yebes40],
    wavelengths=["cm"],
    d_refs="Cernicharo et al. 2021 A&A 648, L3",
    l_refs="Brown et al. 1988 J Am Chem Soc 110, 789",
    notes="",
    Acon=328500,
    Bcon=2539,
    Ccon=2516,
    mua=2.064,
    census_version='2021.0.0',
)
CH2CCH = Molecule(
    name="propargyl radical",
    formula="CH2CCH",
    year=2021,
    label="CH2CCH",
    astromol_name="CH2CCH",
    sources=[TMC1],
    telescopes=[Yebes40],
    wavelengths=["cm"],
    d_refs="Agundez et al. 2021 A&A 647, L10",
    l_refs="Tanaka 1997 J Chem Phys 107, 2728",
    Acon=288055,
    Bcon=9523,
    Ccon=9207,
    mua=0.14,
    notes="Dipole moment taken from Kupper et al. 2002 JCP 117, 647",
    census_version="2021.0.0",
)
HCSCCH = Molecule(
    name="propynethial",
    formula="HCSCCH",
    year=2021,
    label="HCSCCH",
    astromol_name="HCSCCH",
    sources=[TMC1],
    telescopes=[Yebes40],
    wavelengths=["cm"],
    d_refs="Cernicharo et al. 2021 A&A 650, L14",
    d_ref_bib_ids=["2021A&A...650L..14C"],
    l_refs="Brown et al. 1982 Aus. J. Chem. 35, 1747; Crabtree et al. 2016 JCP 144, 124201; Margules et al. 2020 A&A 642, A206",
    l_ref_bib_ids=["Brown:1982ur", "Crabtree:2016fj", "2020A&A...642A.206M"],
    Acon=42652,
    Bcon=3109,
    Ccon=2894,
    mua=1.763,
    mub=0.661,
    census_version='2021.1.0',
)
C5O = Molecule(
    name="pentacarbon monoxide",
    formula="C5O",
    table_formula="C5O",
    year=2021,
    label="C5O",
    astromol_name="C5O",
    sources=[TMC1],
    telescopes=[Yebes40],
    wavelengths=["cm"],
    d_refs="Cernicharo et al. 2021 A&A 656, L21",
    d_ref_bib_ids = ["2021A&A...656L..21C"],
    l_refs="Ogata et al. 1995 JACS 117, 3593",
    l_ref_bib_ids = ["Ogata:1995:3593"],
    notes = "Dipole moment from Botschwina et al. 1995 J Phys Chem 99, 9755",
    Bcon=1367,
    mua=4.06,
    census_version='2021.4.0',
    change_log = {
    				'2021.4.0' : 'Initial entry',
    			}
)
C5Hp = Molecule(
    name="pentynylidyne cation",
    formula="C5H+",
    table_formula="C5H+",
    year=2022,
    label="C5H+",
    astromol_name="C5Hp",
    sources=[TMC1],
    telescopes=[Yebes40],
    wavelengths=["cm"],
    d_refs="Cernicharo et al. 2022 A&A 657, L16",
    d_ref_bib_ids = ["2022A&A...657L..16C"],
    l_refs="Cernicharo et al. 2022 A&A 657, L16",
    l_ref_bib_ids = ["2022A&A...657L..16C"],
    notes = "Identification performed solely on the basis of quantum chemical calculations. Dipole moment from Botschwina 1991 J Chem Phys 95, 4360.",
    Bcon=2412,
    mua=2.88,
    census_version='2021.4.0',
    change_log = {
    				'2021.4.0' : 'Initial entry',
    			}
)
cC5H = Molecule(
    name="ethynylcyclopropenylidene",
    formula="C5H",
    table_formula="c-C5H",
    year=2022,
    label="cC5H",
    astromol_name="cC5H",
    cyclic=True,
    sources=[TMC1],
    telescopes=[Yebes40],
    wavelengths=["cm"],
    d_refs="Cabezas et al. 2022 A&A 663, L2",
    d_ref_bib_ids = ["Cabezas:2022:L2"],
    l_refs="Apponi et al. 2001 ApJL 547, L65; Cabezas et al. 2022 A&A 663, L2",
    l_ref_bib_ids = ["Apponi:2001:L65","Cabezas:2022:L2"],
    notes = "Not 100% certain on the name; Dipole moment from Woon 1995 Chem Phys Lett 244, 45",
    Acon=45018,
    Bcon=3504,
    Ccon=3247,
    mua=4.88,
    census_version='2021.5.0',
    change_log = {
    				'2021.5.0' : 'Initial entry',
    			}
)
HC4S = Molecule(
    name="butadiynethionyl",
    formula="HC4S",
    table_formula="HC4S",
    year=2022,
    label="HC4S",
    astromol_name="HC4S",
    sources=[TMC1],
    telescopes=[Yebes40],
    wavelengths=["cm"],
    d_refs="Fuentetaja et al. 2022 A&AL 667, L4",
    d_ref_bib_ids = ["Fuentetaja:2022:L4"],
    l_refs="Hirahara et al. 1994 JCP 101, 7342",
    l_ref_bib_ids = ["Hirahara:1994:7342"],
    notes = None,
    Bcon=1434,
    mua=1.45,
    census_version='2021.6.0',
    change_log = {
    				'2021.6.0' : 'Initial entry',
    			}
)

######################################################################
#                           Seven Atoms                              #
######################################################################

CH3CHO = Molecule(
    name="acetaldehyde",
    formula="CH3CHO",
    year=1973,
    label="CH3CHO",
    astromol_name="CH3CHO",
    sources=[SgrB2],
    telescopes=[NRAO140],
    wavelengths=["cm"],
    d_refs="Gottlieb 1973 Molecules in the Galactic Environment 181; Fourikis et al. 1974 Aust J Phys 27, 425; Gilmore et al. 1976 ApJ 204, 43",
    l_refs="Kilb et al. 1957 JCP 26, 1695; Souter & Wood 1970 JCP 52, 674",
    notes=None,
    ice="Tentative",
    ice_d_refs="Schutte et al. 1999 A&A 343, 966",
    ice_l_refs="Schutte et al. 1999 A&A 343, 966",
    exgal=True,
    exgal_d_refs="Muller et al. 2011 A&A 535, A103",
    exgal_d_bib_ids=["2011A&A...535A.103M"],
    exgal_sources="PKS 1830-211 LOS",
    Acon=56449,
    Bcon=10160,
    Ccon=9101,
    mua=2.4,
    mub=1.3,
    census_version='2018.0.0',
)
CH3CCH = Molecule(
    name="methylacetylene",
    formula="CH3CCH",
    year=1973,
    label="CH3CCH",
    astromol_name="CH3CCH",
    sources=[SgrB2],
    telescopes=[NRAO36],
    wavelengths=["mm"],
    d_refs="Buhl & Snyder 1973 Molecules in the Galactic Environment 187",
    l_refs="Trambarulo et al. 1950 JCP 18, 1613",
    notes=None,
    exgal=True,
    exgal_d_refs="Mauersberger et al. 1991 A&A 247, 307",
    exgal_d_bib_ids=["1991A&A...247..307M"],
    exgal_sources="NGC 253, M82",
    Acon=158590,
    Bcon=8546,
    Ccon=8546,
    mua=0.8,
    census_version='2018.0.0',
)
CH3NH2 = Molecule(
    name="methylamine",
    formula="CH3NH2",
    year=1974,
    label="CH3NH2",
    astromol_name="CH3NH2",
    sources=[SgrB2, Orion],
    telescopes=[Mitaka6, NRAO36, Parkes64],
    wavelengths=["cm", "mm"],
    d_refs="Fourikis et al. 1974 ApJ 191, L139; Kaifu et al. 1974 ApJ 191, L135",
    l_refs="Takagi & Kojima 1973 ApJ 181, L91",
    notes=None,
    exgal=True,
    exgal_d_refs="Muller et al. 2011 A&A 535, A103",
    exgal_d_bib_ids=["2011A&A...535A.103M"],
    exgal_sources="PKS 1830-211 LOS",
    Acon=103156,
    Bcon=22608,
    Ccon=21730,
    mua=0.3,
    mub=1.3,
    census_version='2018.0.0',
)
CH2CHCN = Molecule(
    name="vinyl cyanide",
    formula="CH2CHCN",
    year=1975,
    label="CH2CHCN",
    astromol_name="CH2CHCN",
    sources=[SgrB2],
    telescopes=[Parkes64],
    wavelengths=["cm"],
    d_refs="Gardner & Winnewisser 1975 ApJ 195, L127",
    l_refs="Gerry & Winnewisser 1973 JMS 48, 1",
    notes=None,
    exgal=True,
    exgal_d_refs="Tercero et al. 2020 A&AL 636, L7",
    exgal_d_bib_ids=["2020A&A...636L...7T"],
    exgal_sources="PKS 1830-211 LOS",
    Acon=49851,
    Bcon=4971,
    Ccon=4514,
    mua=3.8,
    mub=0.9,
    census_version='2018.0.0',
)
HC5N = Molecule(
    name="cyanodiacetylene",
    formula="HC5N",
    year=1976,
    label="HC5N",
    astromol_name="HC5N",
    sources=[SgrB2],
    telescopes=[Algonquin46],
    wavelengths=["cm"],
    d_refs="Broten et al. 1976 ApJ 209, L143; Avery et al. 1976 ApJ 205 L173",
    l_refs="Alexander et al. 1976 JMS 62, 175",
    notes=None,
    exgal="Tentative",
    exgal_d_refs="Aladro et al. 2015 A&A 579, A101",
    exgal_d_bib_ids=["2015A&A...579A.101A"],
    exgal_sources="NGC 253",
    Bcon=1331,
    mua=4.3,
    census_version='2018.0.0',
)
C6H = Molecule(
    name="hexatriynyl radical",
    formula="C6H",
    year=1986,
    label="C6H",
    astromol_name="C6H",
    sources=[TMC1],
    telescopes=[Nobeyama45],
    wavelengths=["cm"],
    d_refs="Suzuki et al. 1986 PASJ 38, 911",
    l_refs="Pearson et al. 1988 A&A 189, L13",
    notes=None,
    Bcon=1391,
    mua=5.5,
    census_version='2018.0.0',
)
cC2H4O = Molecule(
    name="ethylene oxide",
    formula="C2H4O",
    table_formula="c-C2H4O",
    year=1997,
    label="c-C2H4O",
    astromol_name="cC2H4O",
    sources=[SgrB2],
    telescopes=[Haystack37, Nobeyama45, SEST15],
    wavelengths=["cm", "mm"],
    cyclic=True,
    d_refs="Dickens et al. 1997 ApJ 489, 753",
    l_refs="Hirose 1974 ApJ 189, L145",
    notes=None,
    Acon=25484,
    Bcon=22121,
    Ccon=14098,
    mub=1.9,
    census_version='2018.0.0',
)
CH2CHOH = Molecule(
    name="vinyl alcohol",
    formula="CH2CHOH",
    year=2001,
    label="CH2CHOH",
    astromol_name="CH2CHOH",
    sources=[SgrB2],
    telescopes=[NRAOARO12],
    wavelengths=["mm"],
    d_refs="Turner & Apponi 2001 ApJ 561, L207",
    l_refs="Rodler 1985 JMS 114, 23; Kaushik 1977 CPL 49, 90",
    notes=None,
    Acon=62868,
    Bcon=10456,
    Ccon=8963,
    mua=0.5,
    mub=1.7,
    census_version='2018.0.0',
)
C6Hm = Molecule(
    name="hexatriynyl anion",
    formula="C6H-",
    year=2006,
    label="C6H-",
    astromol_name="C6Hm",
    sources=[TMC1, IRC10216],
    telescopes=[GBT],
    wavelengths=["cm"],
    d_refs="McCarthy et al. 2006 ApJ 652, L141",
    l_refs="McCarthy et al. 2006 ApJ 652, L141",
    notes="*First gas-phase molecular anion",
    Bcon=1377,
    mua=8.2,
    census_version='2018.0.0',
)
CH3NCO = Molecule(
    name="methyl isocyanate",
    formula="CH3NCO",
    year=2015,
    label="CH3NCO",
    astromol_name="CH3NCO",
    sources=[SgrB2, Orion],
    telescopes=[NRAOARO12, SMT10],
    wavelengths=["mm"],
    d_refs="Halfen et al. 2015 ApJ 812, L5",
    l_refs="Halfen et al. 2015 ApJ 812, L5",
    notes="*see also Cernicharo et al. 2016 A&A 587, L4",
    Acon=128400,
    Bcon=4415,
    Ccon=4257,
    mua=2.9,
    census_version='2018.0.0',
)
HC5O = Molecule(
    name="butadiynylformyl radical",
    formula="HC5O",
    year=2017,
    label="HC5O",
    astromol_name="HC5O",
    sources=[TMC1],
    telescopes=[GBT],
    wavelengths=["cm"],
    d_refs="McGuire et al. 2017 ApJ 843, L28",
    l_refs="Mohamed et al. 2005 JCP 123, 234301",
    notes=None,
    Bcon=1294,
    mua=2.2,
    census_version='2018.0.0',
)
HOCH2CN = Molecule(
    name="glycolonitrile",
    formula="HOCH2CN",
    year=2019,
    label="HOCH2CN",
    astromol_name="HOCH2CN",
    sources=[IRAS16293],
    telescopes=[ALMA],
    wavelengths=["mm"],
    d_refs="Zeng et al. 2019 MNRAS 484, L43",
    l_refs="Margules et al. 2017 A&A 601, A50",
    notes=None,
    Acon=33610,
    Bcon=4838,
    Ccon=4377,
    mua=2.32,
    mub=1.31,
    muc=1.23,
    census_version='2021.0.0',
)
HC4NC = Molecule(
    name="isocyanoacetylene",
    formula="HC4NC",
    year=2020,
    label="HC4NC",
    astromol_name="HC4NC",
    sources=[TMC1],
    telescopes=[GBT, Yebes40],
    wavelengths=["mm"],
    d_refs="Xue et al. 2020 ApJL 900, L9; Cernicharo et al. 2020 A&A 642, L8",
    l_refs="Botschwina et al. 1998 JCP 109, 3108",
    notes="Also known as isocyanobutadiyne",
    Bcon=1402,
    mua=3.24,
    census_version='2021.0.0',
)
HC3HNH = Molecule(
    name="propargylimine",
    formula="HC3HNH",
    year=2020,
    label="HC3HNH",
    astromol_name="HC3HNH",
    sources=[G0693],
    telescopes=[IRAM30],
    wavelengths=["mm"],
    d_refs="Bizzocchi et al. 2020 A&A 640, A98",
    l_refs="Bizzocchi et al. 2020 A&A 640, A98; Kroto et al. 1984 J. Chem. Soc. Chem. Comm. 993; Sugie et al. 1985 JMS 111, 83; McNaughton et al. 1988 J. Mol. Struct. 190, 195.",
    Acon=54640,
    Bcon=4862,
    Ccon=4458,
    mua=2.14,
    mub=0.17,
    census_version='2021.0.0',
)
C3HCCH = Molecule(
    name="ethynyl cyclopropenylidne",
    formula="C3HCCH",
    table_formula="c-C3HCCH",
    year=2021,
    label="C3HCCH",
    astromol_name="C3HCCH",
    sources=[TMC1],
    telescopes=[Yebes40],
    wavelengths=["cm"],
    d_refs="Cernicharo et al. 2021 A&A 649, L15",
    l_refs="Travers et al. 1997 ApJL 483, L135",
    cyclic=True,
    Acon=34639,
    Bcon=3425,
    Ccon=3114,
    mua=2.04,
    mub=2.89,
    census_version='2021.0.0',
)
MgC5N = Molecule(
    name="magnesium cyanobutadiynyl radical",
    formula="MgC5N",
    table_formula="MgC5N",
    year=2021,
    label="MgC5N",
    astromol_name="MgC5N",
    sources=[IRC10216],
    telescopes=[Yebes40],
    wavelengths=["cm"],
    d_refs="Pardo et al. 2021 A&A 652, L13",
    d_ref_bib_ids = ["2021A&A...652L..13P"],
    l_refs="Pardo et al. 2021 A&A 652, L13",
    l_ref_bib_ids = ["2021A&A...652L..13P"],
    notes = "Identification performed based on quantum chemical calculations only.",
    Bcon=576,
    mua=7.3,
    census_version='2021.1.0',
    change_log = {
    				'2021.1.0' : 'Initial entry',
    			}
)
CH2C3N = Molecule(
    name="3-cyano propargyl radical",
    formula="CH2C3N",
    table_formula="CH2C3N",
    year=2021,
    label="CH2C3N",
    astromol_name="CH2C3N",
    sources=[TMC1],
    telescopes=[Yebes40],
    wavelengths=["cm"],
    d_refs="Cabezas et al. 2021 A&A 654, L9",
    d_ref_bib_ids = ["2021A&A...654L...9C"],
    l_refs="Chen et al. 1998 ApJ 492, 849; Tang et al. 2001 ApJ 552, 409",
    l_ref_bib_ids = ["1998ApJ...492..849C", "2001ApJ...552..409T"],
    notes = None,
    Bcon=2186,
    mua=4.43,
    census_version='2021.4.0',
    change_log = {
    				'2021.4.0' : 'Initial entry',
    			}
)
NC4NHp = Molecule(
    name="protonated butynedinitrile",
    formula="NC4NH+",
    year=2023,
    label="NC4NH+",
    astromol_name="NC4NHp",
    sources=[TMC1],
    telescopes=[Yebes40],
    wavelengths=["cm"],
    d_refs="Agundez et al. 2023 A&AL 669, L1",
    d_ref_bib_ids = ["Agundez:2023:L1"],
    l_refs="Agundez et al. 2023 A&AL 669, L1",
    l_ref_bib_ids = ["Agundez:2023:L1"],
    notes = "No laboratory experiment - identification based on theory.  Dipole moment from Marcelino et al. 2020 A&AL 643, L6.",
    Bcon=1294,
    mua=9.1,
    census_version='2021.6.0',
    change_log = {
    				'2021.6.0' : 'Initial entry',
    			}
)

######################################################################
#                           Eight Atoms                              #
######################################################################

HCOOCH3 = Molecule(
    name="methyl formate",
    formula="HCOOCH3",
    year=1975,
    label="HCOOCH3",
    astromol_name="HCOOCH3",
    sources=[SgrB2],
    telescopes=[Parkes64, Effelsberg100],
    wavelengths=["cm"],
    d_refs="Churchwell & Winnewisser 1975 A&A 45, 229; Brown et al. 1975 ApJ 197, L29",
    l_refs="Brown et al. 1975 ApJ 197, L29",
    notes="The trans conformer was detected in Neill et al. 2012 ApJ 755, 143",
    exgal=True,
    exgal_d_refs="Sewiło et al. 2018 ApJL 853, L19",
    exgal_d_bib_ids=["2018ApJ...853L..19S"],
    exgal_sources="LMC",
    Acon=17630,
    Bcon=9243,
    Ccon=5318,
    mua=1.6,
    mub=0.7,
    census_version='2021.5.0',
    change_log = {
    				'2018.0.0' : 'Initial entry',
    				'2021.5.0' : 'Minor formatting update to notes text'
    			}
)
CH3C3N = Molecule(
    name="methylcyanoacetylene",
    formula="CH3C3N",
    year=1984,
    label="CH3C3N",
    astromol_name="CH3C3N",
    sources=[TMC1],
    telescopes=[NRAO140],
    wavelengths=["cm"],
    d_refs="Broten et al. 1984 ApJ 276, L25",
    l_refs="Moises et al. 1982 JMS 92, 497",
    notes=None,
    Acon=158099,
    Bcon=2066,
    Ccon=2066,
    mua=4.8,
    census_version='2018.0.0',
)
C7H = Molecule(
    name="heptatriynylidyne radical",
    formula="C7H",
    year=1997,
    label="C7H",
    astromol_name="C7H",
    sources=[IRC10216],
    telescopes=[IRAM30],
    wavelengths=["mm"],
    d_refs="Guélin et al. 1997 A&A 317, L1",
    l_refs="Travers et al. 1996 ApJ 465, L77",
    notes=None,
    Bcon=875,
    mua=5.9,
    census_version='2018.0.0',
)
CH3COOH = Molecule(
    name="acetic acid",
    formula="CH3COOH",
    year=1997,
    label="CH3COOH",
    astromol_name="CH3COOH",
    sources=[SgrB2],
    telescopes=[BIMA, OVRO],
    wavelengths=["mm"],
    d_refs="Mehringer et al. 1997 ApJ 480, L71",
    l_refs="Tabor 1957 JCP 27, 974",
    notes=None,
    Acon=11335,
    Bcon=9479,
    Ccon=5325,
    mua=2.9,
    mub=4.9,
    census_version='2018.0.0',
)
H2C6 = Molecule(
    name="hexapentaenylidene",
    formula="H2C6",
    year=1997,
    label="H2C6",
    astromol_name="H2C6",
    sources=[TMC1],
    telescopes=[Goldstone70],
    wavelengths=["cm"],
    d_refs="Langer et al. 1997 ApJ 480, L63",
    l_refs="McCarthy et al. 1997 Science 275, 518",
    notes=None,
    Acon=268400,
    Bcon=1348,
    Ccon=1341,
    mua=6.2,
    census_version='2018.0.0',
)
CH2OHCHO = Molecule(
    name="glycolaldehyde",
    formula="CH2OHCHO",
    year=2000,
    label="CH2OHCHO",
    astromol_name="CH2OHCHO",
    sources=[SgrB2],
    telescopes=[NRAOARO12],
    wavelengths=["mm"],
    d_refs="Hollis et al. 2000 ApJ 540, L107",
    l_refs="Marstokk & Mollendal 1973 J Mol Struct 16, 259",
    notes=None,
    Acon=18446,
    Bcon=6526,
    Ccon=4969,
    mua=0.3,
    mub=2.3,
    census_version='2018.0.0',
)
HC6H = Molecule(
    name="triacetylene",
    formula="HC6H",
    year=2001,
    label="HC6H",
    astromol_name="HC6H",
    sources=[CRL618],
    telescopes=[ISO],
    wavelengths=["IR"],
    d_refs="Cernicharo et al. 2001 ApJ 546, L123",
    l_refs="Haas etal. 1994 JMS 167, 176",
    notes=None,
    exgal=True,
    exgal_d_refs="Bernard-Salas et al. 2006 ApJ 652, L29",
    exgal_d_bib_ids=["2006ApJ...652L..29B"],
    exgal_sources="SMP LMC 11",
    mua=0.0,
    census_version='2018.0.0',
)
CH2CHCHO = Molecule(
    name="propenal",
    formula="CH2CHCHO",
    year=2004,
    label="CH2CHCHO",
    astromol_name="CH2CHCHO",
    sources=[SgrB2, G327306LOS],
    telescopes=[NRAOARO12, SEST15],
    wavelengths=["cm"],
    d_refs="Hollis et al. 2004 ApJ 610, L21",
    l_refs="Winnewisser et al. 1975 Z Naturforsch 30, 1001",
    notes=None,
    Acon=47354,
    Bcon=4660,
    Ccon=4243,
    mua=3.1,
    mub=0.6,
    census_version='2018.0.0',
)
CH2CCHCN = Molecule(
    name="cyanoallene",
    formula="CH2CCHCN",
    year=2006,
    label="CH2CCHCN",
    astromol_name="CH2CCHCN",
    sources=[TMC1],
    telescopes=[GBT],
    wavelengths=["cm"],
    d_refs="Lovas et al. 2006 ApJ 637, L37",
    l_refs="Bouche et al. 1973 J Mol Struct 18, 211",
    notes="*Also Chin et al. 2006 AIP Conf. Proc. 855, 149",
    Acon=25981,
    Bcon=2689,
    Ccon=2475,
    mua=4.1,
    mub=1.3,
    census_version='2018.0.0',
)
NH2CH2CN = Molecule(
    name="aminoacetonitrile",
    formula="NH2CH2CN",
    year=2008,
    label="NH2CH2CN",
    astromol_name="NH2CH2CN",
    sources=[SgrB2],
    telescopes=[IRAM30, PdBI, ATCA],
    wavelengths=["mm"],
    d_refs="Belloche et al. 2008 A&A 482, 179",
    l_refs="Bogey et al. 1990 JMS 143, 180",
    notes=None,
    Acon=30246,
    Bcon=4761,
    Ccon=4311,
    mua=2.6,
    mub=0.6,
    census_version='2018.0.0',
)
CH3CHNH = Molecule(
    name="ethanimine",
    formula="CH3CHNH",
    year=2013,
    label="CH3CHNH",
    astromol_name="CH3CHNH",
    sources=[SgrB2],
    telescopes=[GBT],
    wavelengths=["cm"],
    d_refs="Loomis et al. 2013 ApJL 765, L10",
    l_refs="Loomis et al. 2013 ApJL 765, L10",
    notes=None,
    Acon=49961,
    Bcon=9828,
    Ccon=8650,
    mua=0.8,
    mub=1.9,
    census_version='2018.0.0',
)
CH3SiH3 = Molecule(
    name="methyl silane",
    formula="CH3SiH3",
    year=2017,
    label="CH3SiH3",
    astromol_name="CH3SiH3",
    sources=[IRC10216],
    telescopes=[IRAM30],
    wavelengths=["mm"],
    d_refs="Cernicharo et al. 2017 A&A 606, L5",
    l_refs="Wong et al. 1983 JMS 102, 89",
    notes=None,
    Acon=56189,
    Bcon=10986,
    Ccon=10986,
    mua=0.7,
    census_version='2018.0.0',
)
NH2CONH2 = Molecule(
    name="urea",
    formula="NH2CONH2",
    year=2019,
    label="NH2CONH2",
    astromol_name="NH2CONH2",
    sources=[SgrB2],
    telescopes=[ALMA],
    wavelengths=["mm"],
    d_refs="Belloche et al. 2019 A&A 628, A10",
    l_refs="Brown et al. 1975 JMS 58, 445; Kasten & Dreizler 1986 Z. Naturforsch A. 41, 1173; Kretschmer et al. 1996 Mol. Phys. 87, 1159; Godfrey et al. 1997 J. Mol. Struct. 413-414, 405; Remijan et al. 2014 ApJ 783, 77; Additional work used in Belloche et al. 2019 A&A 628, A10 to be reported in Medvedev et al. in prep as of 9/16/2019.",
    notes="Evidence for the detection, but no claim, made in Remijan et al. 2014 ApJ 783, 77.  Dipole is from Brown et al. 1975 JMS 58, 445",
    Acon=11233,
    Bcon=10369,
    Ccon=5417,
    mub=3.83,
    census_version='2021.0.0',
)
HCCCH2CN = Molecule(
    name="propargyl cyanide",
    formula="HCCCH2CN",
    year=2020,
    label="HCCCH2CN",
    astromol_name="HCCCH2CN",
    sources=[TMC1],
    telescopes=[GBT],
    wavelengths=["cm"],
    d_refs="McGuire et al. 2020 ApJL 900, L10",
    l_refs="Jones & Sheridan 1982 J Mol Struct 78, 303; Demaison et al. 1985 JMS 114, 210; McNaughton et al. 1988 JMS 132, 407; Jager et al. 1990 JMS 143, 50; McGuire et al. 2020 ApJL 900, L10",
    notes="Also known as 3-butynenitrile and 1-cyanoprop-2-yne",
    Acon=19820,
    Bcon=2910,
    Ccon=2573,
    mua=3.23,
    mub=2.34,
    census_version='2021.0.0',
)
CH2CHCCH = Molecule(
    name="vinyl acetylene",
    formula="CH2CHCCH",
    year=2021,
    label="CH2CHCCH",
    astromol_name="CH2CHCCH",
    sources=[TMC1],
    telescopes=[Yebes40],
    wavelengths=["cm"],
    d_refs="Cernicharo et al. 2021 A&A 647, L2",
    l_refs="Thorwirth et al. 2004 J Mol Struct 695, 263; Thorwirth et al. 2003 A&A 398, L11; Sobolev 1961 Optics and Spectroscopy 12, 78",
    notes="Dipole moment is of some debate, especially mub.  See Thorwirth et al. 2003 attempts to improve upon Sobolev 1961.",
    Acon=50300,
    Bcon=4745,
    Ccon=4330,
    mua=0.43,
    mub=0.02,
    census_version='2021.0.0',
)
MgC6H = Molecule(
    name="magnesium hexatriynyl radical",
    formula="MgC6H",
    table_formula="MgC6H",
    year=2021,
    label="MgC6H",
    astromol_name="MgC6H",
    sources=[IRC10216],
    telescopes=[Yebes40],
    wavelengths=["cm"],
    d_refs="Pardo et al. 2021 A&A 652, L13",
    d_ref_bib_ids = ["2021A&A...652L..13P"],
    l_refs="Pardo et al. 2021 A&A 652, L13",
    l_ref_bib_ids = ["2021A&A...652L..13P"],
    notes = "Identification performed based on quantum chemical calculations only.",
    Bcon=580,
    mua=2.5,
    census_version='2021.1.0',
    change_log = {
    				'2021.1.0' : 'Initial entry',
    			}
)
C2H3NH2 = Molecule(
    name="vinylamine",
    formula="C2H3NH2",
    table_formula="C2H3NH2",
    year=2021,
    label="C2H3NH2",
    astromol_name="C2H3NH2",
    sources=[G0693],
    telescopes=[IRAM30,Yebes40],
    wavelengths=["cm","mm"],
    d_refs="Zeng et al. 2021 ApJL 920, L17",
    d_ref_bib_ids = ["2021ApJ...920L..27Z"],
    l_refs="Brown et al. 1990 JMS 142, 195; McNaughton et al. 1994 JMS 163, 80",
    l_ref_bib_ids = ["1990JMoSp.142..195B", "1994JMoSp.163...80M"],
    notes = "Also known as ethenamine",
    Acon=56320,
    Bcon=10035,
    Ccon=8565,
    mua=1.078,
    mub=0.19,
    census_version='2021.4.0',
    change_log = {
    				'2021.4.0' : 'Initial entry',
    			}
)
HOCHCHOH = Molecule(
    name="1,2-ethenediol",
    formula="HOCHCHOH",
    table_formula="HOCHCHOH",
    year=2022,
    label="HOCHCHOH",
    astromol_name="HOCHCHOH",
    sources=[G0693],
    telescopes=[Yebes40,IRAM30],
    wavelengths=["cm","mm"],
    d_refs="Rivilla et al. 2022 ApJL 929, L11",
    d_ref_bib_ids = ["2022ApJ...929L..11R"],
    l_refs="Melosso et al. 2022 Chem Comm 58, 2750",
    l_ref_bib_ids = ["Melosso:2022:2750"],
    notes = None,
    Acon=19507,
    Bcon=6312,
    Ccon=4772,
    mua=1.96,
    mub=0.62,
    census_version='2021.4.0',
    change_log = {
    				'2021.4.0' : 'Initial entry',
    			}
)
HCCCHCCC = Molecule(
    name="ethynylbutatrienyliden",
    formula="HCCCHCCC",
    table_formula="HCCCHCCC",
    year=2022,
    label="HCCCHCCC",
    astromol_name="HCCCHCCC",
    sources=[TMC1],
    telescopes=[Yebes40],
    wavelengths=["cm"],
    d_refs="Fuentetaja et al. 2022 A&AL 667, L4",
    d_ref_bib_ids = ["Fuentetaja:2022:L4"],
    l_refs="McCarthy & Thaddeus 2002 ApJL 569, L55",
    l_ref_bib_ids = ["McCarthy:2002:L55"],
    notes = "Dipole moments from Sattlemeyer & Stanton 2000 JACS 122, 8220",
    Acon=21094,
    Bcon=1676,
    Ccon=1550,
    mua=3.718,
    mub=1.325,
    census_version='2021.6.0',
    change_log = {
    				'2021.6.0' : 'Initial entry',
    			}
)

######################################################################
#                           Nine Atoms                               #
######################################################################

CH3OCH3 = Molecule(
    name="dimethyl ether",
    formula="CH3OCH3",
    year=1974,
    label="CH3OCH3",
    astromol_name="CH3OCH3",
    sources=[Orion],
    telescopes=[NRAO36, NRL85],
    wavelengths=["cm", "mm"],
    d_refs="Snyder et al. 1974 ApJ 191, L79",
    l_refs="Kasai & Myers JCP 30, 1096; Blukis et al. 1963 JCP 38, 2753",
    notes=None,
    exgal=True,
    exgal_d_refs="Qiu et al. 2018 A&A 613, A3; Sewiło et al. 2018 ApJL 853, L19",
    exgal_d_bib_ids=["2018A&A...613A...3Q","2018ApJ...853L..19S"],
    exgal_sources="NGC 1068, LMC",
    Acon=38788,
    Bcon=10057,
    Ccon=8887,
    mub=1.3,
    census_version='2018.0.0',
)
CH3CH2OH = Molecule(
    name="ethanol",
    formula="CH3CH2OH",
    year=1975,
    label="CH3CH2OH",
    astromol_name="CH3CH2OH",
    sources=[SgrB2],
    telescopes=[NRAO36],
    wavelengths=["mm"],
    d_refs="Zukerman et al. 1975 ApJ 196, L99",
    l_refs="Takano et al. 1986 JMS 26, 157",
    notes="*g-ethanol detected 1997 ApJ 480, 420",
    Acon=34892,
    Bcon=9351,
    Ccon=8135,
    mua=0.1,
    mub=1.4,
    census_version='2018.0.0',
)
CH3CH2CN = Molecule(
    name="ethyl cyanide",
    formula="CH3CH2CN",
    year=1977,
    label="CH3CH2CN",
    astromol_name="CH3CH2CN",
    sources=[SgrB2, Orion],
    telescopes=[NRAO36],
    wavelengths=["mm"],
    d_refs="Johnson et al. 1977 ApJ 218, 370",
    l_refs="Johnson et al. 1977 ApJ 218, 370",
    notes=None,
    Acon=27664,
    Bcon=4714,
    Ccon=4235,
    mua=3.9,
    mub=1.2,
    census_version='2018.0.0',
)
HC7N = Molecule(
    name="cyanotriacetylene",
    formula="HC7N",
    year=1977,
    label="HC7N",
    astromol_name="HC7N",
    sources=[TMC1],
    telescopes=[Algonquin46, Haystack37],
    wavelengths=["cm"],
    d_refs="Kroto et al. 1977 Bull. Am. As. Soc. 9, 303; Kroto et al. 1978 ApJ 219, L133",
    l_refs="Kirby et al. 1980 JMS 83, 261",
    notes=None,
    Bcon=564,
    mua=4.8,
    census_version='2021.5.1',
	change_log = {
				'2018.0.0' : 'Initial Entry',
				'2021.5.1' : 'Updated detection references with additional paper.',
    }
)
CH3C4H = Molecule(
    name="methyldiacetylene",
    formula="CH3C4H",
    year=1984,
    label="CH3C4H",
    astromol_name="CH3C4H",
    sources=[TMC1],
    telescopes=[Haystack37, NRAO140, Effelsberg100],
    wavelengths=["cm"],
    d_refs="Walmsley et al. 1984 A&A 134, L11",
    l_refs="Heath et al. 1955 Faraday Discuss. 19, 38",
    notes=None,
    Acon=159140,
    Bcon=2036,
    Ccon=2036,
    mua=1.2,
    census_version='2018.0.0',
)
C8H = Molecule(
    name="octatriynyl radical",
    formula="C8H",
    year=1996,
    label="C8H",
    astromol_name="C8H",
    sources=[IRC10216],
    telescopes=[IRAM30, Nobeyama45],
    wavelengths=["mm"],
    d_refs="Cernicharo & Guélin 1996 A&A 309, L27",
    l_refs="Pauzat et al. 1991 ApJ 369, L13",
    notes=None,
    Bcon=587,
    mua=6.5,
    census_version='2018.0.0',
)
CH3CONH2 = Molecule(
    name="acetamide",
    formula="CH3CONH2",
    year=2006,
    label="CH3CONH2",
    astromol_name="CH3CONH2",
    sources=[SgrB2],
    telescopes=[GBT],
    wavelengths=["cm"],
    d_refs="Hollis et al. 2006 ApJ 643, L25",
    l_refs="Suenram et al. 2001 JMS 208, 188",
    notes=None,
    Acon=10788,
    Bcon=9331,
    Ccon=5157,
    mua=1.1,
    mub=3.5,
    census_version='2018.0.0',
)
C8Hm = Molecule(
    name="octatriynyl anion",
    formula="C8H-",
    year=2007,
    label="C8H-",
    astromol_name="C8Hm",
    sources=[TMC1, IRC10216],
    telescopes=[GBT],
    wavelengths=["cm"],
    d_refs="Brünken et al. 2007 ApJ 664, L43; Remijan et al. 2007 ApJ 664, L47",
    l_refs="Gupta et al. 2007 ApJ 655, L57",
    notes=None,
    Bcon=583,
    mua=10.4,
    census_version='2018.0.0',
)
CH2CHCH3 = Molecule(
    name="propylene",
    formula="CH2CHCH3",
    year=2007,
    label="CH2CHCH3",
    astromol_name="CH2CHCH3",
    sources=[TMC1],
    telescopes=[IRAM30],
    wavelengths=["mm"],
    d_refs="Marcelino et al. 2007 ApJ 665, L127",
    l_refs="Pearson et al. 1994 JMS 166, 120; Wlodarczak et al. 1994 JMS 167, 239",
    notes=None,
    Acon=46281,
    Bcon=9308,
    Ccon=8130,
    mua=0.4,
    mub=0.1,
    census_version='2018.0.0',
)
CH3CH2SH = Molecule(
    name="ethyl mercaptan",
    formula="CH3CH2SH",
    year=2014,
    label="CH3CH2SH",
    astromol_name="CH3CH2SH",
    sources=[Orion],
    telescopes=[IRAM30],
    wavelengths=["mm"],
    d_refs="Kolesniková et al. 2014 ApJ 784, L7",
    l_refs="Kolesniková et al. 2014 ApJ 784, L7",
    notes="Confirmed in Rodríguez-Almeida et al. 2021 ApJL 912, L11",
    Acon=28747,
    Bcon=5295,
    Ccon=4846,
    mua=1.5,
    mub=0.2,
    muc=0.6,
    census_version='2021.0.0',
)
HC7O = Molecule(
    name="hexadiynylformyl radical",
    formula="HC7O",
    year=2017,
    label="HC7O",
    astromol_name="HC7O",
    sources=[TMC1],
    telescopes=[GBT],
    wavelengths=["cm"],
    d_refs="McGuire et al. 2017 ApJ 843, L28",
    l_refs="Mohamed et al. 2005 JCP 123, 234301",
    notes="*Confirmed in Cordiner et al. 2017 ApJ 850, 194",
    Bcon=549,
    mua=2.2,
    census_version='2018.0.0',
)
CH3NHCHO = Molecule(
    name="n-methyl formamide",
    formula="CH3NHCHO",
    year=2021,
    label="CH3NHCHO",
    astromol_name="CH3NHCHO",
    sources=[SgrB2],
    telescopes=[ALMA],
    wavelengths=["mm"],
    d_refs="Belloche et al. 2019 A&A 628, A10",
    l_refs="Belloche et al. 2017 A&A 601, A49",
    notes="Tentatively detected in Belloche et al. 2017 A&A 601, A49.  Rotational constants are for the A-state using an asymmetric top Hamiltonian.  See Kawashima et al. 2010 JMS 263, 11 for a more complete discussion.  Dipole moment is from Fontoni & Caminati 1996 J Chem Soc Faraday Trans 92, 343 who use an opposite naming convention than the others for cis/trans.",
    Acon=19987,
    Bcon=6405,
    Ccon=4902,
    mua=2.914,
    mub=2.41,
    census_version='2021.2.0',
    change_log = {
    				'2021.2.0' : 'Initial entry.'
    				}
)
H2CCCHCCH = Molecule(
    name="allenyl acetylene",
    formula="H2CCCHCCH",
    year=2021,
    label="H2CCCHCCH",
    astromol_name="H2CCCHCCH",
    sources=[TMC1],
    telescopes=[Yebes40],
    wavelengths=["cm"],
    d_refs="Cernicharo et al. 2021 A&A 647, L3",
    l_refs="McCarthy et al. 2020 J Phys Chem A 124, 5170; Lee & McCarthy 2019 J Phys Chem Lett 10, 2408",
    notes="",
    Acon=25961,
    Bcon=2616,
    Ccon=2413,
    mua=0.630,
    mub=0.011,
    census_version='2021.0.0',
)
HCCCHCHCN = Molecule(
    name="cyanovinylacetylene",
    formula="HCCCHCHCN",
    year=2021,
    label="HCCCHCHCN",
    astromol_name="HCCCHCHCN",
    sources=[TMC1],
    telescopes=[GBT],
    wavelengths=["cm"],
    d_refs="Lee et al. 2021 ApJL 908, L11",
    l_refs="Halter et al. 2001 J Am Chem Soc 123, 12353; Thorwirth et al. 2004 J Mol Spectrosc 225, 93; McCarthy et al. 2020 J Phys Chem A 124, 5170",
    notes="",
    Acon=7098,
    Bcon=2683,
    Ccon=1943,
    mua=4.2,
    mub=0.6,
    census_version='2021.0.0',
)
H2CCHC3N = Molecule(
    name="vinylcyanoacetylene",
    formula="H2CCHC3N",
    year=2021,
    label="H2CCHC3N",
    astromol_name="H2CCHC3N",
    sources=[TMC1],
    telescopes=[GBT],
    wavelengths=["cm"],
    d_refs="Lee et al. 2021 ApJL 908, L11",
    l_refs="Halter et al. 2001 J Am Chem Soc 123, 12353; Thorwirth et al. 2004 J Mol Spectrosc 225, 93; McCarthy et al. 2020 J Phys Chem A 124, 5170",
    notes="",
    Acon=39920,
    Bcon=1377,
    Ccon=1330,
    mua=5.3,
    mub=0.3,
    census_version='2021.0.0',
)

######################################################################
#                             Ten Atoms                              #
######################################################################

CH3COCH3 = Molecule(
    name="acetone",
    formula="CH3COCH3",
    year=1987,
    label="CH3COCH3",
    astromol_name="CH3COCH3",
    sources=[SgrB2],
    telescopes=[IRAM30, NRAO140, NRAOARO12],
    wavelengths=["cm", "mm"],
    d_refs="Combes et al. 1987 A&A 180, L13",
    l_refs="Vacherand et al. 1986 JMS 118, 355",
    notes="*Confirmed in 2002 ApJ 578, 245",
    Acon=10165,
    Bcon=8515,
    Ccon=4910,
    mub=2.9,
    census_version='2018.0.0',
)
HOCH2CH2OH = Molecule(
    name="ethylene glycol",
    formula="HOCH2CH2OH",
    year=2002,
    label="HOCH2CH2OH",
    astromol_name="HOCH2CH2OH",
    sources=[SgrB2],
    telescopes=[NRAOARO12],
    wavelengths=["mm"],
    d_refs="Hollis et al. 2002 ApJ 571, L59",
    l_refs="Christen et al. 1995 JMS 172, 57",
    notes="*aGg' conformer in 2017 A&A 598, A59",
    Acon=15361,
    Bcon=5588,
    Ccon=4614,
    mua=2.1,
    mub=0.9,
    census_version='2018.0.0',
)
CH3CH2CHO = Molecule(
    name="propanal",
    formula="CH3CH2CHO",
    year=2004,
    label="CH3CH2CHO",
    astromol_name="CH3CH2CHO",
    sources=[SgrB2],
    telescopes=[GBT],
    wavelengths=["cm"],
    d_refs="Hollis et al. 2004 ApJ 610, L21",
    l_refs="Butcher & Wilson 1964 JCP 40, 1671",
    notes=None,
    Acon=16712,
    Bcon=5969,
    Ccon=4648,
    mua=1.7,
    mub=1.9,
    census_version='2018.0.0',
)
CH3C5N = Molecule(
    name="methylcyanodiacetylene",
    formula="CH3C5N",
    year=2006,
    label="CH3C5N",
    astromol_name="CH3C5N",
    sources=[TMC1],
    telescopes=[GBT],
    wavelengths=["cm"],
    d_refs="Snyder et al. 2006 ApJ 647, 412",
    l_refs="Chen et al. 1998 JMS 192, 1",
    notes=None,
    Acon=158099,
    Bcon=778,
    Ccon=778,
    mua=5.4,
    census_version='2018.0.0',
)
CH3CHCH2O = Molecule(
    name="propylene oxide",
    formula="CH3CHCH2O",
    year=2016,
    label="CH3CHCH2O",
    astromol_name="CH3CHCH2O",
    sources=[SgrB2],
    telescopes=[GBT],
    wavelengths=["cm"],
    cyclic=True,
    d_refs="McGuire & Carroll et al. 2016 Science 352, 1449",
    l_refs="McGuire & Carroll et al. 2016 Science 352, 1449",
    notes="*First chiral molecule",
    Acon=18024,
    Bcon=6682,
    Ccon=5951,
    mua=1.0,
    mub=1.7,
    muc=0.6,
    census_version='2018.0.0',
)
CH3OCH2OH = Molecule(
    name="methoxymethanol",
    formula="CH3OCH2OH",
    year=2017,
    label="CH3OCH2OH",
    astromol_name="CH3OCH2OH",
    sources=[NGC6334],
    telescopes=[ALMA],
    wavelengths=["mm"],
    d_refs="McGuire et al. 2017 ApJ 851, L46",
    l_refs="Motiyenko et al. 2018 PCCP 20, 5509",
    notes=None,
    Acon=17238,
    Bcon=5568,
    Ccon=4813,
    mua=0.2,
    mub=0.1,
    muc=0.1,
    census_version='2018.0.0',
)
C6H4 = Molecule(
    name="ortho-benzyne",
    formula="C6H4",
    year=2021,
    label="C6H4",
    astromol_name="C6H4",
    sources=[TMC1],
    telescopes=[Yebes40, IRAM30],
    wavelengths=["cm"],
    d_refs="Cernicharo et al. 2021 A&A 652, L9",
    d_ref_bib_ids=["2021A&A...652L...9C"],
    l_refs="Brown et al. 1986 JACS 108, 1296; Kukolich et al. 2003 JCP 119, 4353; Robertson et al. 2003 JMS 217, 123",
    l_ref_bib_ids=["Brown:1986lp", "2003JChPh.119.4353K", "2003JMoSp.217..123R"],
    cyclic=True,
    Acon=6990,
    Bcon=5707,
    Ccon=3140,
    mub=1.38,
    census_version='2021.1.0',
    change_log = {
    				'2021.1.0' : 'Initial Entry',
    }
)

C2H5NCO = Molecule(
    name="ethyl isocyanate",
    formula="C2H5NCO",
    year=2021,
    label="C2H5NCO",
    astromol_name="C2H5NCO",
    sources=[G0693],
    telescopes=[IRAM30, Yebes40],
    wavelengths=["cm"],
    d_refs="Rodriguez-Almeida et al. 2021 A&A 654, L1",
    d_ref_bib_ids = ["2021A&A...654L...1R"],
    l_refs="Sakaizumi et al. 1976 Bull. Chem. Soc. Japan 49, 2908; Heineking et al. 1994 Mol. Phys. 81, 1177; Kolesnikova et al. 2018 A&A 616, A173",
    l_ref_bib_ids = ["Sakaizumi:1976uu", "Heineking:1994op", "2018A&A...616A.173K"],
    notes="",
    Acon=339670,
    Bcon=38086,
    Ccon=34124,
    mua=3.83,
    census_version='2021.3.0',
    change_log = {
    				'2021.3.0' : 'Initial Entry',
    }
)
HC7NHp = Molecule(
    name="",
    formula="HC7NH+",
    table_formula="HC7NH+",
    year=2022,
    label="HC7NH+",
    astromol_name="HC7NHp",
    sources=[TMC1],
    telescopes=[Yebes40],
    wavelengths=["cm"],
    d_refs="Cabezas et al. 2022 A&A 659, L8",
    d_ref_bib_ids = ["2022A&A...659L...8C"],
    l_refs="Cabezas et al. 2022 A&A 659, L8",
    l_ref_bib_ids = ["2022A&A...659L...8C"],
    notes = "Identification based solely on quantum chemical calculations.",
    Bcon=554,
    mua=6.4,
    census_version='2021.4.0',
    change_log = {
    				'2021.4.0' : 'Initial entry',
    			}
)
CH3CHCHCN = Molecule(
    name="crotonitrile",
    formula="CH3CHCHCN",
    table_formula="CH3CHCHCN",
    year=2022,
    label="CH3CHCHCN",
    astromol_name="CH3CHCHCN",
    sources=[TMC1],
    telescopes=[Yebes40],
    wavelengths=["cm"],
    d_refs="Cernicharo et al. 2022 A&A 663, L5",
    d_ref_bib_ids = ["Cernicharo:2022:L5"],
    l_refs="Lesarri et al. 1995 JMS 172, 520; McCarthy et al. 2020 J Phys Chem A 124, 5170",
    l_ref_bib_ids = ["Lesarri:1995:520","McCarthy:2020:5170"],
    notes = "Both the trans- and cis-conformers are reported. Parameters here are for the A-state of the trans version, which is lower in energy.  Dipole moments obtained from Beaudet 1963 J Chem Phys 38, 2548 and Suzuki & Kozima 1970 JMS 33, 407",
    Acon=38054,
    Bcon=2297,
    Ccon=2195,
    mua=4.35,
    census_version='2021.5.0',
    change_log = {
    				'2021.5.0' : 'Initial entry',
    			}
)
CH2CCH3CN = Molecule(
    name="methacrylonitrile",
    formula="CH2CCH3CN",
    table_formula="CH2CCH3CN",
    year=2022,
    label="CH2CCH3CN",
    astromol_name="CH2CCH3CN",
    sources=[TMC1],
    telescopes=[Yebes40],
    wavelengths=["cm"],
    d_refs="Cernicharo et al. 2022 A&A 663, L5",
    d_ref_bib_ids = ["Cernicharo:2022:L5"],
    l_refs="Lopez et al. 1990 JMS 141, 317; Lesarri et al. 1995 JMS 172, 520",
    l_ref_bib_ids = ["Lopez:1990:317","Lesarri:1995:520"],
    notes = "",
    Acon=9291,
    Bcon=4166,
    Ccon=2925,
    mua=3.940,
    mub=0.26,
    census_version='2021.5.0',
    change_log = {
    				'2021.5.0' : 'Initial entry',
    			}
)
CH2CHCH2CN = Molecule(
    name="allyl cyanide",
    formula="CH2CHCH2CN",
    table_formula="CH2CHCH2CN",
    year=2022,
    label="CH2CHCH2CN",
    astromol_name="CH2CHCH2CN",
    sources=[TMC1],
    telescopes=[Yebes40],
    wavelengths=["cm"],
    d_refs="Cernicharo et al. 2022 A&A 663, L5",
    d_ref_bib_ids = ["Cernicharo:2022:L5"],
    l_refs="Demaison et al. 1991 JMS 146, 455; McCarthy et al. 2020 J Phys Chem A 124, 5170",
    l_ref_bib_ids = ["Demaison:1991:455","McCarthy:2020:5170"],
    notes = "Both the gauche and cis forms were found; parameters here are for the cis, which is the lower energy conformer.  Dipole moments from Sastry et al. 1968 Can. J. Phys. 46, 959",
    Acon=11323,
    Bcon=3739,
    Ccon=2858,
    mua=3.26,
    mub=2.16,
    census_version='2021.5.0',
    change_log = {
    				'2021.5.0' : 'Initial entry',
    			}
)

######################################################################
#                           Eleven Atoms                             #
######################################################################

HC9N = Molecule(
    name="cyanotetraacetylene",
    formula="HC9N",
    year=1978,
    label="HC9N",
    astromol_name="HC9N",
    sources=[TMC1],
    telescopes=[Algonquin46, NRAO140],
    wavelengths=["cm"],
    d_refs="Broten et al. 1978 ApJ 223, L105",
    l_refs="Iida et al. 1991 ApJ 371, L45",
    notes=None,
    Bcon=291,
    mua=5.2,
    census_version='2018.0.0',
)
CH3C6H = Molecule(
    name="methyltriacetylene",
    formula="CH3C6H",
    year=2006,
    label="CH3C6H",
    astromol_name="CH3C6H",
    sources=[TMC1],
    telescopes=[GBT],
    wavelengths=["cm"],
    d_refs="Remijan et al. 2006 ApJ 643, L37",
    l_refs="Alexander et al. 1978 JMS 70, 84",
    notes=None,
    Acon=159140,
    Bcon=778,
    Ccon=778,
    mua=1.5,
    census_version='2018.0.0',
)
C2H5OCHO = Molecule(
    name="ethyl formate",
    formula="C2H5OCHO",
    year=2009,
    label="C2H5OCHO",
    astromol_name="C2H5OCHO",
    sources=[SgrB2],
    telescopes=[IRAM30],
    wavelengths=["mm"],
    d_refs="Belloche et al. 2009 A&A 499, 215",
    l_refs="Medvedev et al. 2009 ApJS 181, 433",
    notes=None,
    Acon=17747,
    Bcon=2905,
    Ccon=2579,
    mua=1.9,
    mub=0.7,
    muc=0.0,
    census_version='2018.0.0',
)
CH3COOCH3 = Molecule(
    name="methyl acetate",
    formula="CH3COOCH3",
    year=2013,
    label="CH3COOCH3",
    astromol_name="CH3COOCH3",
    sources=[Orion],
    telescopes=[IRAM30],
    wavelengths=["mm"],
    d_refs="Tercero et al. 2013 ApJ 770, L13",
    l_refs="Tudorie et al. 2011 JMS 269, 211",
    notes=None,
    Acon=10247,
    Bcon=4170,
    Ccon=3077,
    mua=0.0,
    mub=1.6,
    census_version='2018.0.0',
)
CH3COCH2OH = Molecule(
    name="hydroxyacetone",
    formula="CH3COCH2OH",
    year=2020,
    label="CH3COCH2OH",
    astromol_name="CH3COCH2OH",
    sources=[IRAS16293],
    telescopes=[ALMA],
    wavelengths=["mm"],
    d_refs="Zhou et al. 2020 Res. Astron. & Astrophys. 20, 125",
    l_refs="Kattija-Ari & Harmony et al. 1980 International Journal of Quantum Chemistry 18, 443; Apponi et al. 2006 ApJ 652, 1787; Braakman et al. 2010 JMS 264, 43",
    notes=None,
    Acon=10074,
    Bcon=3817,
    Ccon=2867,
    mua=2.22,
    mub=2.17,
    census_version='2021.0.0',
)
C5H6 = Molecule(
    name="cyclopentadiene",
    formula="C5H6",
    year=2021,
    label="C5H6",
    astromol_name="C5H6",
    sources=[TMC1],
    telescopes=[Yebes40],
    wavelengths=["cm"],
    d_refs="Cernicharo et al. 2021 A&A 649, L15",
    l_refs="Laurie 1956 J Chem Phys 24, 635; Scharpen & Laurie 1965 J Chem Phys 43, 2765; Benson & Flygare 1970 J Am Chem Soc 92, 7523; Bogey et al. 1988 J Mol Spectrosc 132, 277",
    cyclic=True,
    Acon=8426,
    Bcon=8226,
    Ccon=4271,
    mub=0.416,
    census_version='2021.0.0',
)
NH2CH2CH2OH = Molecule(
    name="ethanolamine",
    formula="NH2CH2CH2OH",
    year=2021,
    label="NH2CH2CH2OH",
    astromol_name="NH2CH2CH2OH",
    sources=[G0693],
    telescopes=[Yebes40, IRAM30],
    wavelengths=["cm","mm"],
    d_refs="Rivilla et al. 2021 PNAS 118, e2101314118",
    d_ref_bib_ids=["2021PNAS..11801314R"],
    l_refs="Penn & Curl 1971 JCP 53, 651; Kaushik & Woods 1982 Z Physik Chemie Neue Folge 132, 117; Widicus Weaver et al. 2003 JMS 217, 278",
    l_ref_bib_ids=["1971JChPh..55..651P", "Kaushik:1982ld", "2003JMoSp.217..278W"],
    Acon=14509,
    Bcon=5546,
    Ccon=4570,
    mua=2.65,
    mub=0.89,
    muc=0.42,
    census_version='2021.1.0',
)
CH2CCHC4H = Molecule(
    name="allenyl diacetylene",
    formula="CH2CCHC4H",
    table_formula="CH2CCHC4H",
    year=2022,
    label="CH2CCHC4H",
    astromol_name="CH2CCHC4H",
    sources=[TMC1],
    telescopes=[Yebes40],
    wavelengths=["cm"],
    d_refs="Fuentetaja et al. 2022 A&A 663, L3",
    d_ref_bib_ids = ["Fuentetaja:2022:L3"],
    l_refs="McCarthy et al. 2020 J Phys Chem A 124, 5170",
    l_ref_bib_ids = ["McCarthy:2020:5170"],
    notes = "Dipole moment from Lee & McCarthy 2019 J Phys Chem Lett 10, 2408",
    Acon=15656,
    Bcon=948,
    Ccon=898,
    mua=1.20,
    mub=0.20,
    census_version='2021.5.0',
    change_log = {
    				'2021.5.0' : 'Initial entry',
    			}
)

######################################################################
#                           Twelve Atoms                             #
######################################################################

C6H6 = Molecule(
    name="benzene",
    formula="C6H6",
    year=2001,
    label="C6H6",
    astromol_name="C6H6",
    sources=[CRL618],
    telescopes=[ISO],
    wavelengths=["IR"],
    cyclic=True,
    d_refs="Cernicharo et al. 2001 ApJ 546, L123",
    l_refs="Lindenmayer et al. 1988 JMS 128 172",
    notes=None,
    exgal=True,
    exgal_d_refs="Bernard-Salas et al. 2006 ApJ 652, L29",
    exgal_d_bib_ids=["2006ApJ...652L..29B"],
    exgal_sources="SMP LMC 11",
    mua=0.0,
    census_version='2018.0.0',
)
nC3H7CN = Molecule(
    name="n-propyl cyanide",
    formula="C3H7CN",
    table_formula="n-C3H7CN",
    year=2009,
    label="n-C3H7CN",
    astromol_name="nC3H7CN",
    sources=[SgrB2],
    telescopes=[IRAM30],
    wavelengths=["mm"],
    d_refs="Belloche et al. 2009 A&A 499, 215",
    l_refs="Belloche et al. 2009 A&A 499, 215",
    notes=None,
    Acon=23668,
    Bcon=2268,
    Ccon=2153,
    mua=4.0,
    mub=1.0,
    muc=0.0,
    census_version='2018.0.0',
)
iC3H7CN = Molecule(
    name="isopropyl cyanide",
    formula="C3H7CN",
    table_formula="i-C3H7CN",
    year=2014,
    label="i-C3H7CN",
    astromol_name="iC3H7CN",
    sources=[SgrB2],
    telescopes=[ALMA],
    wavelengths=["mm"],
    d_refs="Belloche et al. 2014 Science 345, 1584",
    l_refs="Muller et al. 2011 JMS 267, 100",
    notes=None,
    Acon=7941,
    Bcon=3968,
    Ccon=2901,
    mua=4.0,
    mub=0.6,
    census_version='2018.0.0',
)
C5H5CN1 = Molecule(
    name="1-cyano-1,3-cyclopentadiene",
    formula="C5H5CN",
    table_formula="1-C5H5CN",
    year=2021,
    label="1-C5H5CN",
    astromol_name="C5H5CN1",
    sources=[TMC1],
    telescopes=[GBT],
    wavelengths=["cm"],
    cyclic=True,
    d_refs="McCarthy et al. 2021 Nature Astronomy 5, 176",
    l_refs="McCarthy et al. 2021 Nature Astronomy 5, 176; Ford & Seitzman 1978 JMS 69, 326; Sakaizumi et al. 1987 Bull. Chem. Soc. Jap. 60, 3903",
    l_ref_bib_ids=["1978JMoSp..69..326F","1987BCSJ..60..3903"],
    notes="First molecule with a 5-membered ring detected in the ISM.",
    Acon=8353,
    Bcon=1904,
    Ccon=1565,
    mua=4.15,
    mub=0.27,
    census_version='2021.0.0',
)
C5H5CN2 = Molecule(
    name="2-cyano-1,3-cyclopentadiene",
    formula="C5H5CN",
    table_formula="2-C5H5CN",
    year=2021,
    label="2-C5H5CN",
    astromol_name="C5H5CN2",
    sources=[TMC1],
    telescopes=[GBT],
    wavelengths=["cm"],
    cyclic=True,
    d_refs="Lee et al. 2021 ApJL 910, L2",
    l_refs="Lee et al. 2021 ApJL 910, L2; McCarthy et al. 2021 Nature Astronomy 5, 176; Ford & Seitzman 1978 JMS 69, 326; Sakaizumi et al. 1987 Bull. Chem. Soc. Jap. 60, 3903",
    notes="",
    Acon=8236,
    Bcon=1902,
    Ccon=1560,
    mua=4.36,
    mub=0.77,
    census_version='2021.0.0',
)
nCH3CH2CH2OH = Molecule(
    name="n-propanol",
    formula="CH3CH2CH2OH",
    table_formula="n-CH3CH2CH2OH",
    year=2022,
    label="nCH3CH2CH2OH",
    astromol_name="nCH3CH2CH2OH",
    sources=[G0693],
    telescopes=[IRAM30,Yebes40],
    wavelengths=["cm","mm"],
    d_refs="Jimenez-Serra et al. 2022 A&A 663, A181",
    d_ref_bib_ids = ["Jimenez-Serra:2022:A181"],
    l_refs="Kisiel et al. 2010 PCCP 12, 8329; Dreizier & Scappini 1981 Z Naturforsch A 36, 1187",
    l_ref_bib_ids = ["Kisiel:2020:8329","Dreizier:1981:1187"],
    notes = "",
    Acon=14330,
    Bcon=5119,
    Ccon=4324,
    mua=0.4914,
    mub=0.9705,
    muc=0.9042,
    census_version='2021.5.0',
    change_log = {
    				'2021.5.0' : 'Initial entry',
    			}
)
iCH3CH2CH2OH = Molecule(
    name="i-propanol",
    formula="CH3CH2CH2OH",
    table_formula="i-CH3CH2CH2OH",
    year=2022,
    label="iCH3CH2CH2OH",
    astromol_name="iCH3CH2CH2OH",
    sources=[SgrB2],
    telescopes=[ALMA],
    wavelengths=["mm"],
    d_refs="Belloche et al. 2022 A&A 662, A110",
    d_ref_bib_ids = ["Belloche:2022:A110"],
    l_refs="Maeda et al. 2006 ApJS 166, 650",
    l_ref_bib_ids = ["Maeda:2006:650"],
    notes = "Both gauche and anti-conformers are detected, but not cataloged separately here.  Parameters given are for the gauche, which is lower in energy.",
    Acon=8640,
    Bcon=8063,
    Ccon=4768,
    mua=1.114,
    mub=0.737,
    muc=0.813,
    census_version='2021.5.0',
    change_log = {
    				'2021.5.0' : 'Initial entry',
    			}
)
######################################################################
#                         Thirteen Atoms                             #
######################################################################

cC6H5CN = Molecule(
    name="benzonitrile",
    formula="C6H5CN",
    year=2018,
    label="C6H5CN",
    astromol_name="cC6H5CN",
    sources=[TMC1],
    telescopes=[GBT],
    wavelengths=["cm"],
    cyclic=True,
    d_refs="McGuire et al. 2018 Science 359, 202",
    l_refs="Wohlfart et al. 2008 JMS 247, 119",
    notes=None,
    Acon=5655,
    Bcon=1547,
    Ccon=1214,
    mua=4.5,
    census_version='2018.0.0',
)
HC11N = Molecule(
    name="cyanopentaacetylene",
    formula="HC11N",
    year=2021,
    label="HC11N",
    astromol_name="HC11N",
    sources=[TMC1],
    telescopes=[GBT],
    wavelengths=["cm"],
    d_refs="Loomis et al. 2021 Nature Astronomy 5, 188",
    l_refs="Travers et al. 1996 ApJL 469, L65",
    notes=None,
    Bcon=169,
    mua=5.47,
    census_version='2021.0.0',
)
cC5H4CCH2 = Molecule(
    name="fulveneallene",
    formula="C5H4CCH2",
    table_formula="c-C5H4CCH2",
    year=2022,
    label="cC5H4CCH2",
    astromol_name="cC5H4CCH2",
    cyclic=True,
    sources=[TMC1],
    telescopes=[Yebes40],
    wavelengths=["cm"],
    d_refs="Cernicharo et al. 2022 A&A 663, L9",
    d_ref_bib_ids = ["Cernicharo:2022:L9"],
    l_refs="Sakaizumi et al. 1993 JMS 159, 112; McCarthy et al. 2020 J Phys Chem A 124, 5170",
    l_ref_bib_ids = ["Sakaizumi:1993:112", "McCarthy:2020:5170"],
    notes = "",
    Acon=8181,
    Bcon=1887,
    Ccon=1549,
    mua=0.69,
    census_version='2021.5.0',
    change_log = {
    				'2021.5.0' : 'Initial entry',
    			}
)
######################################################################
#                                PAHs                                #
######################################################################

CNN1 = Molecule(
    name="1-cyanonaphthalene",
    formula="C10H7CN",
    table_formula="1-C10H7CN",
    year=2021,
    label="CNN1",
    astromol_name="CNN1",
    sources=[TMC1],
    telescopes=[GBT],
    wavelengths=["cm"],
    cyclic=True,
    pah=True,
    d_refs="McGuire et al. 2021 Science 371, 1265",
    l_refs="McNaughton et al. 2018 MNRAS 476, 5268",
    notes="First individually detected PAH molecule in the ISM, alongside 2-cyanonaphthalene.",
    Acon=1479,
    Bcon=957,
    Ccon=581,
    mua=3.56,
    mub=2.96,
    census_version='2021.0.0',
)
CNN2 = Molecule(
    name="2-cyanonaphthalene",
    formula="C10H7CN",
    table_formula="2-C10H7CN",
    year=2021,
    label="CNN2",
    astromol_name="CNN2",
    sources=[TMC1],
    telescopes=[GBT],
    wavelengths=["cm"],
    cyclic=True,
    pah=True,
    d_refs="McGuire et al. 2021 Science 371, 1265",
    l_refs="McNaughton et al. 2018 MNRAS 476, 5268",
    notes="First individually detected PAH molecule in the ISM, alongside 1-cyanonaphthalene.",
    Acon=2707,
    Bcon=606,
    Ccon=495,
    mua=5.09,
    mub=0.98,
    census_version='2021.0.0',
)
C9H8 = Molecule(
    name="indene",
    formula="C9H8",
    year=2021,
    label="C9H8",
    astromol_name="C9H8",
    sources=[TMC1],
    telescopes=[GBT, Yebes40],
    wavelengths=["cm"],
    cyclic=True,
    pah=True,
    d_refs="Burkhardt et al. 2021 ApJL 913, L18; Cernicharo et al. 2021 A&AL 649, 15",
    l_refs="Burkhardt et al. 2021 ApJL 913, L18; Li et al. 1979 J Mol Struct 51, 171",
    notes="",
    Acon=3775,
    Bcon=1581,
    Ccon=1122,
    mua=0.59,
    mub=0.43,
    census_version='2021.0.0',
)

######################################################################
#                            Fullerenes                             #
######################################################################

C60 = Molecule(
    name="buckminsterfullerene",
    formula="C60",
    year=2010,
    label="C60",
    astromol_name="C60",
    sources=[TC1, NGC7023],
    telescopes=[Spitzer],
    wavelengths=["IR"],
    fullerene=True,
    cyclic=True,
    d_refs="Cami et al. 2010 Science 329, 1180",
    l_refs="Nemes et al. 1994 CPL 218, 295",
    notes="*See also Sellgren et al. 2010 ApJ 722, L54 and Werner 2004b, Sellgren 2007 therein",
    mua=0.0,
    census_version='2018.0.0',
)
C60p = Molecule(
    name="buckminsterfullerene cation",
    formula="C60+",
    year=2013,
    label="C60+",
    astromol_name="C60p",
    sources=[NGC7023],
    telescopes=[Spitzer],
    wavelengths=["IR"],
    fullerene=True,
    cyclic=True,
    d_refs="Berné et al. 2013 A&A 550, L4",
    l_refs="Kern et al. 2013 JPCA 117, 8251",
    notes="*See also Campbell et al. 2015 Nature 523, 322",
    mua=0.0,
    census_version='2018.0.0',
)
C70 = Molecule(
    name="rugbyballene",
    formula="C70",
    year=2010,
    label="C70",
    astromol_name="C70",
    sources=[TC1],
    telescopes=[Spitzer],
    wavelengths=["IR"],
    fullerene=True,
    cyclic=True,
    d_refs="Cami et al. 2010 Science 329, 1180",
    l_refs="Nemes et al. 1994 CPL 218, 295",
    notes=None,
    mua=0.0,
    census_version='2018.0.0',
)

#############################################################
# 							Full List						#
#############################################################

all_molecules = [
    # two atoms
    CH,
    CN,
    CHp,
    OH,
    CO,
    H2,
    SiO,
    CS,
    SO,
    SiS,
    NS,
    C2,
    NO,
    HCl,
    NaCl,
    AlCl,
    KCl,
    AlF,
    PN,
    SiC,
    CP,
    NH,
    SiN,
    SOp,
    COp,
    HF,
    N2,
    CFp,
    PO,
    O2,
    AlO,
    CNm,
    OHp,
    SHp,
    HClp,
    SH,
    TiO,
    ArHp,
    NSp,
    HeHp,
    VO,
    POp,
    SiP,
    # three atoms
    H2O,
    HCOp,
    HCN,
    OCS,
    HNC,
    H2S,
    N2Hp,
    C2H,
    SO2,
    HCO,
    HNO,
    HCSp,
    HOCp,
    SiC2,
    C2S,
    C3,
    CO2,
    CH2,
    C2O,
    MgNC,
    NH2,
    NaCN,
    N2O,
    MgCN,
    H3p,
    SiCN,
    AlNC,
    SiNC,
    HCP,
    CCP,
    AlOH,
    H2Op,
    H2Clp,
    KCN,
    FeCN,
    HO2,
    TiO2,
    CCN,
    SiCSi,
    S2H,
    HCS,
    HSC,
    NCO,
    CaNC,
    NCS,
    MgC2,
    # four atoms
    NH3,
    H2CO,
    HNCO,
    H2CS,
    C2H2,
    C3N,
    HNCS,
    HOCOp,
    C3O,
    lC3H,
    HCNHp,
    H3Op,
    C3S,
    cC3H,
    HC2N,
    H2CN,
    SiC3,
    CH3,
    C3Nm,
    PH3,
    HCNO,
    HOCN,
    HSCN,
    HOOH,
    lC3Hp,
    HMgNC,
    HCCO,
    CNCN,
    HONO,
    MgCCH,
    HCCS,
    HNCN,
    H2NC,
    HCCSp,
    # five atoms
    HC3N,
    HCOOH,
    CH2NH,
    NH2CN,
    H2CCO,
    C4H,
    SiH4,
    cC3H2,
    CH2CN,
    C5,
    SiC4,
    H2CCC,
    CH4,
    HCCNC,
    HNCCC,
    H2COHp,
    C4Hm,
    CNCHO,
    HNCNH,
    CH3O,
    NH3Dp,
    H2NCOp,
    NCCNHp,
    CH3Cl,
    MgC3N,
    HC3Op,
    NH2OH,
    HC3Sp,
    H2CCS,
    C4S,
    CHOSH,
    HCSCN,
    HC3O,
    # six atoms
    CH3OH,
    CH3CN,
    NH2CHO,
    CH3SH,
    C2H4,
    C5H,
    CH3NC,
    HC2CHO,
    H2C4,
    C5S,
    HC3NHp,
    C5N,
    HC4H,
    HC4N,
    cH2C3O,
    CH2CNH,
    C5Nm,
    HNCHCN,
    SiH3CN,
    MgC4H,
    CH3COp,
    H2CCCS,
    CH2CCH,
    HCSCCH,
    C5O,
    C5Hp,
    cC5H,
    HC4S,
    # seven atoms
    CH3CHO,
    CH3CCH,
    CH3NH2,
    CH2CHCN,
    HC5N,
    C6H,
    cC2H4O,
    CH2CHOH,
    C6Hm,
    CH3NCO,
    HC5O,
    HOCH2CN,
    HC4NC,
    HC3HNH,
    C3HCCH,
    MgC5N,
    CH2C3N,
    NC4NHp,
    # eight atoms
    HCOOCH3,
    CH3C3N,
    C7H,
    CH3COOH,
    H2C6,
    CH2OHCHO,
    HC6H,
    CH2CHCHO,
    CH2CCHCN,
    NH2CH2CN,
    CH3CHNH,
    CH3SiH3,
    NH2CONH2,
    HCCCH2CN,
    CH2CHCCH,
    MgC6H,
    C2H3NH2,
    HOCHCHOH,
    HCCCHCCC,
    # nine atoms
    CH3OCH3,
    CH3CH2OH,
    CH3CH2CN,
    HC7N,
    CH3C4H,
    C8H,
    CH3CONH2,
    C8Hm,
    CH2CHCH3,
    CH3CH2SH,
    HC7O,
    CH3NHCHO,
    H2CCCHCCH,
    HCCCHCHCN,
    H2CCHC3N,
    # ten atoms
    CH3COCH3,
    HOCH2CH2OH,
    CH3CH2CHO,
    CH3C5N,
    CH3CHCH2O,
    CH3OCH2OH,
    C6H4,
    C2H5NCO,
    HC7NHp,
    CH3CHCHCN,
    CH2CCH3CN,
    CH2CHCH2CN,
    # eleven atoms
    HC9N,
    CH3C6H,
    C2H5OCHO,
    CH3COOCH3,
    CH3COCH2OH,
    C5H6,
    NH2CH2CH2OH,
    CH2CCHC4H,
    # twelve atoms
    C6H6,
    nC3H7CN,
    iC3H7CN,
    C5H5CN1,
    C5H5CN2,
    nCH3CH2CH2OH,
    iCH3CH2CH2OH,
    # thirteen atoms
    cC6H5CN,
    HC11N,
    cC5H4CCH2,
    # PAHS
    CNN1,
    CNN2,
    C9H8,
    # fullerenes
    C60,
    C60p,
    C70,
]
