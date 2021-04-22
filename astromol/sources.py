class Source(object):

	def __init__(self,name=None,type=None,ra=None,dec=None,detects=0,mols=None,simbad_url=None):
	
		self.name = name
		self.type = type
		self.ra = ra
		self.dec = dec
		self.detects = detects
		self.mols = mols
		self.simbad_url = simbad_url
					
		return	
		
	def update_stats(self,list):	
		for x in list:		
			if self in x.sources:			
				if self.mols is None:
					self.mols = [x]					
				else:				
					self.mols.append(x)					
				self.detects += 1				
		return
		
'''
Source coordinates are generalized for simplicity, and individual detections may have 
been made (read: absolutely have been made) toward various pointing positions within 
these sources.  

Similarly, source types are highly generalized, to allow for some aggregate analysis.  

Finally, some detections have been made along the line of sight to these sources.  In 
these cases,  the source is catagorized as a 'LOS Cloud', regardless of the actual source
type, as it was only used as an absorbing background.

Diffuse, translucent, and dense clouds are all classified this way for simplicity.

Please check the individual papers for details.
'''

AFGL890LOS = Source(name="AFGL 890 LOS",
					type='LOS Cloud',
					ra='06:10:48.0',
					dec='-06:12:00',
					simbad_url='http://simbad.u-strasbg.fr/simbad/sim-id?Ident=AFGL+890')

AFGL961LOS = Source(name="AFGL 961 LOS",
					type='LOS Cloud',
					ra='06:34:37.741',
					dec='+04:12:44.20',
					simbad_url='http://simbad.u-strasbg.fr/simbad/sim-id?Ident=AFGL961')

AFGL989LOS = Source(name="AFGL 989 LOS",
					type='LOS Cloud',
					ra='06:41:10.06',
					dec='+09:29:35.8',
					simbad_url='http://simbad.u-strasbg.fr/simbad/sim-id?Ident=AFGL989')

B1b = Source(name="B1-b",
			type='Dark Cloud',
			ra='03:33:20.8',
			dec='+31:07:34',
			simbad_url='http://simbad.u-strasbg.fr/simbad/sim-id?Ident=%5BHKM99%5D+B1-b')

CRL2688 = Source(name="CRL 2688",
					type='Carbon Star',
					ra='21:02:18.27',
					dec='+36:41:37.0',
					simbad_url='http://simbad.u-strasbg.fr/simbad/sim-id?Ident=CRL+2688')

CRL618 = Source(name="CRL 618",
				type='Carbon Star',
				ra='04:42:53.62',
				dec='+36:06:53.40',
				simbad_url='http://simbad.u-strasbg.fr/simbad/sim-id?Ident=CRL+618')

CasALOS = Source(name="Cas A LOS",
					type='LOS Cloud',
					ra='23:23:24.00',
					dec='+58:48:54.0',
					simbad_url='http://simbad.u-strasbg.fr/simbad/sim-id?Ident=Cas+A')

CrabNebula = Source(name="Crab Nebula",
					type='SNR',
					ra='05:34:31.94',
					dec='+22:00:52.2',
					simbad_url='http://simbad.u-strasbg.fr/simbad/sim-id?Ident=Crab+Nebula')

CygnusOB212LOS = Source(name="Cygnus OB2 - 12",
						type='LOS Cloud',
						ra='20:32:40.96',
						dec='+41:04:13.2',
						simbad_url='http://simbad.u-strasbg.fr/simbad/sim-id?Ident=%4011680932&Name=Schulte%2012')

DR21 = Source(name="DR 21",
				type='SFR',
				ra='20:39:01.6',
				dec='+42:19:38',
				simbad_url='http://simbad.u-strasbg.fr/simbad/sim-basic?Ident=DR+21')

DR21LOS = Source(name="DR 21 LOS",
					type='LOS Cloud',
					ra='20:39:01.6',
					dec='+42:19:38',
					simbad_url='http://simbad.u-strasbg.fr/simbad/sim-basic?Ident=DR+21')

DR21OH = Source(name="DR 21(OH)",
				type='SFR',
				ra='20:39:01.01',
				dec='+42:22:50.22',
				simbad_url='http://simbad.u-strasbg.fr/simbad/sim-basic?Ident=DR+21%28OH%29')

G0693 = Source(name="G+0.693-0.027",
				type='Shock',
				ra='17:47:21.86',
				dec='-28:22:43.00')

G327306LOS = Source(name="G327.3-0.6 LOS",
					type='LOS Cloud',
					ra='15:53:05.0',
					dec='-54:35:24',
					simbad_url='http://simbad.u-strasbg.fr/simbad/sim-id?Ident=G327.3-0.6')

GL2136LOS = Source(name="GL2136 LOS",
					type='LOS Cloud',
					ra='18:27:18.43',
					dec='-25:04:02.84',
					simbad_url='http://simbad.u-strasbg.fr/simbad/sim-id?Ident=%402520798&Name=GJ%20%202136%20B')

GalacticCenter = Source(name="Galactic Center",
						type='SFR')

HD124314LOS = Source(name="HD 124314 LOS",
						type='LOS Cloud',
						ra='14:15:01.61',
						dec='-61:42:24.38',
						simbad_url='http://simbad.u-strasbg.fr/simbad/sim-id?Ident=HD+124314')

HD27778LOS = Source(name="HD 27778 LOS",
					type='LOS Cloud',
					ra='04:23:59.78',
					dec='24:18:03.53',
					simbad_url='http://simbad.u-strasbg.fr/simbad/sim-id?Ident=HD+27778')

HorseheadPDR = Source(name="Horsehead PDR",
						type='PDR',
						ra='05:40:53.936',
						dec='-02:28:00',
						simbad_url='http://simbad.u-strasbg.fr/simbad/sim-id?Ident=%40828287&Name=NAME%20Horsehead%20Nebula')

IC443G = Source(name="IC 443G",
				type='SNR',
				ra='06:16:43.4',
				dec='+22:32:24',
				simbad_url='http://simbad.u-strasbg.fr/simbad/sim-id?Ident=IC+443G')

IRAS16293 = Source(name="IRAS 16293",
					type='Protostar',
					ra='16:32:22.56',
					dec='-24:28:31.8',
					simbad_url='http://simbad.u-strasbg.fr/simbad/sim-id?Ident=IRAS+16293-2422')

IRC10216 = Source(name="IRC+10216",
					type='Carbon Star',
					ra='09:47:57.406',
					dec='+13:16:43.56',
					simbad_url='http://simbad.u-strasbg.fr/simbad/sim-id?Ident=IRC%2B10216')

K350 = Source(name="K3-50",
				type='HII',
				ra='20:04:45.59',
				dec='33:32:42.0',
				simbad_url='http://simbad.u-strasbg.fr/simbad/sim-id?Ident=%402904320&Name=NAME%20K%203-50A')

L134 = Source(name="L134",
				type='Dark Cloud',
				ra='15:53:36.3',
				dec='-04:35:26.0',
				simbad_url='http://simbad.u-strasbg.fr/simbad/sim-id?Ident=%402622026&Name=LDN%20%20134')

L1527 = Source(name="L1527",
				type='Dark Cloud',
				ra='04:39:53.0',
				dec='+25:45:00',
				simbad_url='http://simbad.u-strasbg.fr/simbad/sim-id?Ident=L1527')

L1544 = Source(name="L1544",
				type='Dark Cloud',
				ra='05:04:16.6',
				dec='+25:10:48',
				simbad_url='http://simbad.u-strasbg.fr/simbad/sim-id?Ident=L1544')

L183 = Source(name="L183",
				type='Dark Cloud',
				ra='15:54:12.2',
				dec='-02:49:42',
				simbad_url='http://simbad.u-strasbg.fr/simbad/sim-basic?Ident=L183')

L483 = Source(name="L483",
				type='Dark Cloud',
				ra='18:17:35.0',
				dec='-04:39:48',
				simbad_url='http://simbad.u-strasbg.fr/simbad/sim-id?Ident=L483')

LOSCloud = Source(name="LOS Cloud",
					type='LOS Cloud')

Lupus1A = Source(name="Lupus-1A",
					type='Dark Cloud',
					ra='15:42:52.4',
					dec='-34:07:53.5',
					simbad_url='http://simbad.u-strasbg.fr/simbad/sim-id?Ident=%405549180&Name=NAME%20Lupus-1A')

M17LOS = Source(name="M17 LOS",
				type='LOS Cloud',
				ra='18:20:47',
				dec='-16:10:18',
				simbad_url='http://simbad.u-strasbg.fr/simbad/sim-id?Ident=M17')

M17SW = Source(name="M17SW",
				type='PDR',
				ra='18:20:23.1',
				dec='-16:11:43',
				simbad_url='http://simbad.u-strasbg.fr/simbad/sim-id?Ident=M17+SW')

M3LOS = Source(name="M3 LOS",
				type='LOS Cloud',
				ra='13:42:11.62',
				dec='+28:22:38.2',
				simbad_url='http://simbad.u-strasbg.fr/simbad/sim-id?Ident=M3')

NGC2024 = Source(name="NGC 2024",
					type='PDR',
					ra='05:41:43',
					dec='-01:50:30',
					simbad_url='http://simbad.u-strasbg.fr/simbad/sim-id?Ident=NGC+2024')

NGC2024LOS = Source(name="NGC 2024 LOS",
					type='LOS Cloud',
					ra='05:41:43',
					dec='-01:50:30',
					simbad_url='http://simbad.u-strasbg.fr/simbad/sim-id?Ident=NGC+2024')

NGC2264 = Source(name="NGC 2264",
					type='YSO',
					ra='06:40:58',
					dec='+09:53:42',
					simbad_url='http://simbad.u-strasbg.fr/simbad/sim-id?Ident=NGC+2264')

NGC6334 = Source(name="NGC 6334",
					type='SFR',
					ra='17:20:53.3',
					dec='-35:46:59',
					simbad_url='http://simbad.u-strasbg.fr/simbad/sim-id?Ident=%402361705&Name=NAME%20NGC%206334-I')

NGC6334LOS = Source(name="NGC 6334 LOS",
					type='LOS Cloud',
					ra='17:20:53.3',
					dec='-35:46:59',
					simbad_url='http://simbad.u-strasbg.fr/simbad/sim-id?Ident=%402361705&Name=NAME%20NGC%206334-I')

NGC7023 = Source(name="NGC 7023",
					type='PDR',
					ra='21:01:36.9',
					dec='+68:09:48',
					simbad_url='http://simbad.u-strasbg.fr/simbad/sim-id?Ident=NGC+7023')

NGC7027 = Source(name="NGC 7027",
					type='PN',
					ra='21:07:01.8',
					dec='+42:14:10',
					simbad_url='http://simbad.u-strasbg.fr/simbad/sim-id?Ident=NGC+7023')

NGC7538 = Source(name="NGC 7538",
					type='YSO',
					ra='23:13:37.2',
					dec='61:30:00',
					simbad_url='http://simbad.u-strasbg.fr/simbad/sim-id?Ident=NGC+7538')

NGC7538LOS = Source(name="NGC 7538 LOS",
					type='LOS Cloud',
					ra='23:13:37.2',
					dec='61:30:00',
					simbad_url='http://simbad.u-strasbg.fr/simbad/sim-id?Ident=NGC+7538')

Orion = Source(name="Orion",
					type='SFR',
					ra='05:35:14.16',
					dec='-05:22:21.5',
					simbad_url='http://simbad.u-strasbg.fr/simbad/sim-id?Ident=Orion+KL')

OrionBar = Source(name="Orion Bar",
					type='PDR',
					ra='05:35:22.30',
					dec='-05:24:33.0',
					simbad_url='http://simbad.u-strasbg.fr/simbad/sim-id?Ident=Orion+Bar')

rhoOphA = Source(name="rho Ophiuchi A",
					type='SFR',
					ra='16:26:27.20',
					dec='-24:24:04',
					simbad_url='http://simbad.u-strasbg.fr/simbad/sim-id?Ident=%402488602&Name=NAME%20rho%20Oph%20A%20SM%201')

SgrA = Source(name="Sgr A",
					type='Sgr A',
					ra='17:45:40.0',
					dec='-29:00:28.2',
					simbad_url='http://simbad.u-strasbg.fr/simbad/sim-id?Ident=Sgr+A')

SgrALOS = Source(name="Sgr A LOS",
					type='LOS Cloud',
					ra='17:45:40.0',
					dec='-29:00:28.2',
					simbad_url='http://simbad.u-strasbg.fr/simbad/sim-id?Ident=Sgr+A')

SgrB2 = Source(name="Sgr B2",
					type='SFR',
					ra='17:47:20.4',
					dec='-28:23:07',
					simbad_url='http://simbad.u-strasbg.fr/simbad/sim-id?Ident=Sgr+B2')

SgrB2LOS = Source(name="Sgr B2 LOS",
					type='LOS Cloud',
					ra='17:47:20.4',
					dec='-28:23:07',
					simbad_url='http://simbad.u-strasbg.fr/simbad/sim-id?Ident=Sgr+B2')

TC1 = Source(name="TC 1",
					type='PN',
					ra='17:45:35.29',
					dec='-46:05:23.7',
					simbad_url='http://simbad.u-strasbg.fr/simbad/sim-id?Ident=PN%20Tc%201%20')

TMC1 = Source(name="TMC-1",
					type='Dark Cloud',
					ra='04:41:45.9',
					dec='+25:41:27',
					simbad_url='http://simbad.u-strasbg.fr/simbad/sim-id?Ident=TMC-1')

VYCaMaj = Source(name="VY Ca Maj",
					type='Oxygen Star',
					ra='07:22:58.3',
					dec='-25:46:03.2',
					simbad_url='http://simbad.u-strasbg.fr/simbad/sim-id?Ident=VY+Canis+Majoris')

W3 = Source(name="W3",
					type='LOS Cloud',
					ra='02:27:04.10',
					dec='+61:52:27.1',
					simbad_url='http://simbad.u-strasbg.fr/simbad/sim-id?Ident=W3')

W3OH = Source(name="W3(OH)",
					type='SFR',
					ra='02:27:04.1',
					dec='+61:52:52',
					simbad_url='http://simbad.u-strasbg.fr/simbad/sim-id?Ident=W3+%28OH%29&NbIdent=1')

W31LOS = Source(name="W31 LOS",
					type='LOS Cloud',
					ra='18:10:28.6',
					dec='-19:55:51',
					simbad_url='http://simbad.u-strasbg.fr/simbad/sim-id?Ident=W31')

W33LOS = Source(name="W33 LOS",
					type='LOS Cloud',
					ra='18:14:14.0',
					dec='-17:55:50',
					simbad_url='http://simbad.u-strasbg.fr/simbad/sim-id?Ident=W33')

W43LOS = Source(name="W43 LOS",
					type='LOS Cloud',
					ra='18:47:32.4',
					dec='-01:56:31',
					simbad_url='http://simbad.u-strasbg.fr/simbad/sim-id?Ident=W43')

W44LOS = Source(name="W44 LOS",
					type='LOS Cloud',
					ra='18:56:10.65',
					dec='+01:13:21.30',
					simbad_url='http://simbad.u-strasbg.fr/simbad/sim-id?Ident=W43')

W49 = Source(name="W49",
					type='SFR',
					ra='19:10:19.6',
					dec='+09:07:42',
					simbad_url='http://simbad.u-strasbg.fr/simbad/sim-id?Ident=W49')

W49LOS = Source(name="W49 LOS",
					type='LOS Cloud',
					ra='19:10:19.6',
					dec='+09:07:42',
					simbad_url='http://simbad.u-strasbg.fr/simbad/sim-id?Ident=W49')

W51 = Source(name="W51",
					type='SFR',
					ra='19:23:50',
					dec='+14:06:0',
					simbad_url='http://simbad.u-strasbg.fr/simbad/sim-id?Ident=W51')

W51LOS = Source(name="W51 LOS",
					type='LOS Cloud',
					ra='19:23:50',
					dec='+14:06:0',
					simbad_url='http://simbad.u-strasbg.fr/simbad/sim-id?Ident=W51')

XiPerLOS = Source(name="Xi Per LOS",
					type='LOS Cloud',
					ra='03:58:57.9',
					dec='+35:47:27.74',
					simbad_url='http://simbad.u-strasbg.fr/simbad/sim-id?Ident=Xi+Per')

rhoOphA = Source(name="rho Oph A",
					type='SFR',
					ra='16:25:35.14',
					dec='-23:26:49.9',
					simbad_url='http://simbad.u-strasbg.fr/simbad/sim-id?Ident=Rho+Oph+A')

source_list = [
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
		