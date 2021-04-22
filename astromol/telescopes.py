'''
An exceptional resource for these is 'Observatories and Telescopes of Modern Times' 
by David Leverington.
'''

class Telescope(object):

	def __init__(self,name=None,shortname=None,type=None,wavelength=None,latitude=None,longitude=None,diameter=None,built=None,decommissioned=None,notes=None):
	
		self.name = name
		self.shortname = shortname
		self.type = type
		self.wavelength = wavelength
		self.latitude = latitude
		self.longitude = longitude
		self.diameter = diameter
		self.built = built
		self.decommissioned = decommissioned
		self.notes = notes
		self.ndetects = None
		self.mol_list = None
	
		return
		
	def update_stats(self,mol_list):
	
		'''
		Takes a list of molecules to loop over and add to the list of molecules detected with this telescope, as well as the detects counter.
		'''
		
		my_mols = []	
		for mol in mol_list:		
			if self in mol.telescopes:			
				my_mols.append(mol)
				
		self.mol_list = my_mols				
		self.ndetects = len(my_mols)
					
		return	


GBT = Telescope(name='Green Bank Telescope',
				shortname='GBT',
				type='Single Dish',
				wavelength=['cm','mm'],
				latitude=38.433056,
				longitude=-79.839722,
				diameter=100,
				built=2004)
				
IRAM30 = Telescope(name='IRAM 30-m',
					shortnam='IRAM',
					type='Single Dish',
					wavelength=['mm','sub-mm'],
					latitude=37.066161,
					longitude=-3.392719,
					diameter=30,
					built=1984)
					
Spitzer = Telescope(name='Spizter',
					shortname='Spitzer',
					type='Space',
					wavelength=['IR'],
					diameter=0.85,
					built=2003,
					decommissioned=2020)
					
Hubble = Telescope(name='Hubble Space Telescope',
					shortname='Hubble',
					type='Space',
					wavelength=['IR','Vis', 'UV'],
					diameter=2.4,
					built=1990)
					
ALMA = Telescope(name='Atacama Large Millimeter/sub-millimeter Array',
					shortname='ALMA',
					type='Interferometer',
					wavelength=['mm','sub-mm'],
					latitude=-23.0193,
					longitude=-67.7532,
					built=2011)
					
ISO = Telescope(name='Infrared Space Observatory',
				shortname='ISO',
				type='Space',
				wavelength=['IR'],
				diameter=0.6,
				built=1995,
				decommissioned=1998)
				
NRAO140 = Telescope(name='NRAO 140-ft',
					shortname='NRAO 140-ft',
					type='Single Dish',
					wavelength=['cm'],
					latitude=38.433056,
					longitude=-79.839722,
					diameter=43,
					built=1965,
					decommissioned=2008,
					notes='Technically started operations again in 2014, but not for PI science.')
					
Algonquin46 = Telescope(name='Algonquin 46-m Telescope',
						shortname='Algonquin 46-m',
						type='Single Dish',
						wavelength=['cm','mm'],
						latitude=45.955503,
						longitude=-78.073042,
						diameter=46,
						built=1966,
						decommissioned=1987,
						notes='Technically in operation much longer, but seems to have ceased PI science in 1987.')
						
NRAOARO12 = Telescope(name='NRAO/ARO 12-m Telescope',
						shortname='NRAO/ARO 12-m',
						type='Single Dish',
						wavelength=['mm'],
						latitude=31.9533,
						longitude=-111.615,
						diameter=12,
						built=1984,
						notes='Originally the NRAO 36-ft (11 m) telescope until 1984, when the dish was replaced (these are counted as separate facilities).  The observatory was handed over to the ARO in 2000 and renamed.  In 2013, the antenna was replaced with a 12-m ALMA prototype antenna.')
						
NRAO36 = Telescope(name='NRAO 36-ft Telescope',
					shortname='NRAO 36-ft',
					type='Single Dish',
					wavelength=['mm'],
					latitude=31.9533,
					longitude=-111.615,
					diameter=11,
					built=1967,
					decommissioned=1984,
					notes='Became the NRAO/ARO 12-m telescope in 1984.')
					
Nobeyama45 = Telescope(name='Nobeyama 45-m Telescope',
						shortname='Nobeyama',
						type='Single Dish',
						wavelength=['cm','mm'],
						latitude=35.9417,
						longitude=138.4758,
						diameter=45,
						built=1982)
						
Effelsberg100 = Telescope(name='Effelsberg 100-m Telescope',
							shortname='Effelsberg',
							type='Single Dish',
							wavelength=['cm'],
							latitude=50.5247,
							longitude=-6.8828,
							diameter=100,
							built=1972)
							
Haystack37 = Telescope(name='Haystack 37-m Telescope',
						shortname='Haystack',
						type='Single Dish',
						wavelength=['cm','mm'],
						latitude=42.6233,
						longitude=-71.4882,
						diameter=37,
						built=1964)
						
PdBI = Telescope(name='Plateu de Bure Interferometer',
					shortname='PdBI',
					type='Interferometer',
					wavelength=['mm'],
					latitude=44.63389,
					longitude=5.90792,
					built=1988,
					decommissioned=2016)
					
NOEMA = Telescope(name='Northern Extended Millimeter Array',
					shortname='NOEMA',
					type='Interferometer',
					wavelength=['mm'],
					latitude=44.63389,
					longitude=5.90792,
					built=2016)
					
BIMA = Telescope(name='Berkeley-Illinois-Maryland Array',
					shortname='BIMA',
					type='Interferometer',
					wavelength=['mm'],
					latitude=40.8178,
					longitude=-121.473,
					built=1986,
					decommissioned=2005,
					notes='Became part of CARMA.')
					
OVRO = Telescope(name='Caltech Owens Valley Radio Observatory Millimeter Array',
					shortname='OVRO',
					type='Interferometer',
					wavelength=['mm'],
					latitude=37.2339,
					longitude=-118.282,
					built=1984,
					decommissioned=2005,
					notes='Became part of CARMA.')
					
Yebes40 = Telescope(name='Yebes RT40-m Telescope',
						shortname='Yebes',
						type='Single Dish',
						wavelength=['cm','mm'],
						latitude=40525208,
						longitude=-3.088725,
						built=2007)
						
NRL85 = Telescope(name='Maryland Point Observatory Naval Research Lab 85-foot Telescope',
					shortname='NRL 85-ft',
					type='Single Dish',
					wavelength=['cm'],
					diameter=26,
					latitude=38.3741667,
					longitude=-77.230833,
					built=1965,
					decommissioned=1994,
					notes='Primarily used for VLBI for much of its later years.')
					
ATCA = Telescope(name='Australia Telescope Compact Array',
					shortname='ATCA',
					type='Interferometer',
					wavelength=['cm'],
					latitude=-30.312778,
					longitude=149.550278,
					built=1988)

Parkes64 = Telescope(name='Parkes 64-m Telescope',
						shortname='Parkes',
						type='Single Dish',
						wavelength=['cm'],
						diameter=64,
						latitude=-32.99778,
						longitude=148.26292,
						built=1961)

SMT10 = Telescope(name='ARO 10-m Submillimeter Telescope',
					shortname='SMT',
					type='Single Dish',
					wavelength=['mm','sub-mm'],
					diameter=10,
					latitude=32.701658,
					longitude=-109.871391,
					built=1993)

SEST15 = Telescope(name='Swedish-ESO 15-m Submillimetre Telescope',
					shortname='SEST',
					type='Single Dish',
					wavelength=['mm','sub-mm'],
					diameter=15,
					latitude=-29.26,
					longitude=-70.73,
					built=1987,
					decommissioned=2003)

Goldstone70 = Telescope(name='Goldstone 72-m (DSS-14; "Mars")',
						shortname='Goldstone',
						type='Single Dish',
						wavelength=['cm'],
						diameter=70,
						latitude=35.426667,
						longitude=-116.89,
						built=1966,
						notes='Originally a 64-m dish; become 70-m in 1988. Conceivably still PI Science Capable?')

Mitaka6 = Telescope(name='Tokyo Astronomical Observatory Mitaka 6-m',
					shortname='Mitaka 6-m',
					type='Single Dish',
					wavelength=['mm'],
					diameter=6,
					latitude=35.675217,
					longitude=139.538083,
					built=1970,
					decommissioned=2018,
					notes='Moved around quite a bit within Japan until returning (and retiring) in 2018.')

McMath = Telescope(name='McMath-Pierce Solar Telescope',
					shortname='McMath Solar Telescope',
					type='Optical',
					wavelength=['IR','Vis','UV'],
					diameter=1.6,
					latitude=31.9584,
					longitude=-111.595,
					built=1962)

Bell7m = Telescope(name='AT&T Bell Laboratories 7-m Telescope',
					shortname='Bell 7-m',
					type='Single Dish',
					wavelength=['cm'],
					diameter=7,
					built=1976,
					decommissioned=1992)

IRTF = Telescope(name='NASA Infrared Telescope Facility',
					shortname='IRTF',
					type='Optical',
					wavelength=['IR'],
					diameter=3,
					latitude=19.8263,
					longitude=-155.473,
					built=1974)

KPNO4m = Telescope(name='Mayall 4-m Telescope',
					shortname='KPNO 4-m',
					type='Optical',
					wavelength=['IR'],
					diameter=4,
					latitude=31.9583,
					longitude=-111.5967,
					built=1973)

Onsala20m = Telescope(name='Onsala 20-m Telescope',
						shortname='Onsala 20-m',
						type='Single Dish',
						wavelength=['cm','mm'],
						diameter=20,
						latitude=57.393056,
						longitude=11.917778,
						built=1976)

FCRAO14m = Telescope(name='Five College Radio Observatory 14-m Telescope',
						shortname='FCRAO 14-m',
						type='Single Dish',
						wavelength=['cm','mm'],
						latitude=42.391925,
						longitude=-72.344097,
						built=1976,
						decommissioned=2005)

APEX = Telescope(name='Atacama Pathfinder Experiment',
					shortname='APEX',
					type='Single Dish',
					wavelength=['mm', 'sub-mm'],
					diameter=12,
					latitude=-23.0058,
					longitude=-67.7592,
					built=2005)

CSO = Telescope(name='Caltech Submillimeter Observatory',
				shortname='CSO',
				type='Single Dish',
				wavelength=['mm', 'sub-mm'],
				diameter=10.4,
				latitude=19.8225,
				longitude=-155.70694,
				built=1986,
				decommissioned=2015)

MWO4m = Telescope(name='University of Texas Millimeter Wave Observatory 4.9-m Telescope',
					shortname='MWO 4.9-m',
					type='Single Dish',
					wavelength=['mm'],
					latitude=30.3866,
					longitude=-97.7269,
					built=1971,
					decommissioned=1988)

HatCreek = Telescope(name='Hat Creek Station 20-ft Telescope',
						shortname='Hat Creek 20-ft',
						type='Single Dish',
						wavelength=['cm','mm'],
						latitude=40.8178,
						longitude=-121.473,
						built=1965,
						decommissioned=1983,
						notes='Best build date found was "mid 1960s", so 1965 is an estimate.  It appears to have either been subsummed into BIMA or decommissioned when BIMA came online.  The decommissioning date is an estimate.')

SMA = Telescope(name='Submillimeter Array',
				shortname='SMA',
				type='Interferometer',
				wavelength=['mm'],
				latitude=19.8225,
				longitude=-155.70694,
				built=2003)

Herschel = Telescope(name='Herschel Space Telescope',
						shortname='Herschel',
						type='Space',
						wavelength=['sub-mm','IR'],
						built=2009,
						decommissioned=2013)

UKIRT = Telescope(name='United Kingdom Infrared Telescope',
					shortname='UKIRT',
					type='Optical',
					wavelength=['IR'],
					diameter=3.8,
					latitude=19.8225,
					longitude=-155.70694,
					built=1979)

SOFIA = Telescope(name='Stratospheric Observatory for Infrared Astronomy',
					shortname='SOFIA',
					type='Airborne',
					wavelength=['sub-mm','IR'],
					diameter=2.5,
					built=2010)

Odin = Telescope(name='Odin',
					shortname='Odin',
					type='Space',
					wavelength=['sub-mm'],
					diameter=1.1,
					built=2001)

FUSE = Telescope(name='Far Ultraviolet Spectroscopic Explorer',
					shortname='FUSE',
					type='Space',
					wavelength=['UV'],
					built=1999,
					decommissioned=2007)

Kuiper = Telescope(name='Kuiper Airborne Observatory',
					shortname='KAO',
					type='Airborne',
					wavelength=['sub-mm','IR'],
					built=1974,
					decommissioned=1995)

MtHopkins = Telescope(name='Tillinghast 60 inch',
						shortname='Mt. Hopkins 60-in',
						type='Optical',
						diameter=1.5,
						wavelength=['IR'],
						built=1969,
						latitude=31.6811,
						longitude=-110.878)

Aerobee = Telescope(name='Aerobee-150 Rocket',
					shortname='Aerobee-150 Rocket',
					type='Airborne',
					wavelength=['UV'],
					built=1970,
					decommissioned=1970,
					notes='They literally put a spectrometer on a rocket and shot it into the sky, then used the same spectrometer to measure H2 in the laboratory.')

Millstone = Telescope(name='Lincoln Laboratory Millstone Hill Observatory 84-ft',
						shortname='Millstone Hill 84-ft',
						type='Single Dish',
						wavelength=['cm'],
						built=1956,
						decommissioned=1978,
						diameter=26,
						latitude=42.6233,
						longitude=-71.4882,
						notes='Originally the Ballistic Missile Early Warning System radar antenna.  Decommissioning date is a best guess based on the installation of larger telescopes to the site at that time.')

MtWilson = Telescope(name='Mount Wilson 100-in',
						shortname='Mt. Wilson',
						type='Optical',
						wavelength=['UV','VIS'],
						diameter=2.54,
						built=1917,
						decommissioned=1989,
						notes='Now known as the Hooker Telescope.')

IRAS = Telescope(name='Infrared Astronomical Satellite',
					shortname='IRAS',
					type='Space',
					wavelength=['IR'],
					diameter=0.60,
					built=1983,
					decommissioned=1983)



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
	NOEMA,
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
