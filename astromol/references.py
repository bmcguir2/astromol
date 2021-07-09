import bibtexparser as btp
from pathlib import Path
module_path = Path(__file__)

####################################################################################################
#                                      Extragalactic Molecules                                     # 
####################################################################################################         

#load in the bibtex database
with open(module_path.parents[0].joinpath('bibtex/exgal_refs.bib'), 'r') as input:
    exgal_bib_database = btp.bparser.BibTexParser(common_strings=True).parse_file(input) 
    #split the authors into a list 
    for entry in exgal_bib_database.entries:
        entry = btp.customization.author(entry)


#get rid of {} in authors list and re-assign to an actual dictionary with the ADS keys as the dictionary keys.
exgal_bibs = {}
for entry in exgal_bib_database.entries:
    entry['author'] = [x.replace("}","").replace("{","") for x in entry['author']]
    exgal_bibs[entry['ID']] = entry

####################################################################################################
#                                  Protoplanetary Disk Molecules                                   # 
####################################################################################################         

#load in the bibtex database
with open(module_path.parents[0].joinpath('bibtex/ppd_refs.bib'), 'r') as input:
    ppd_bib_database = btp.bparser.BibTexParser(common_strings=True).parse_file(input) 
    #split the authors into a list 
    for entry in ppd_bib_database.entries:
        entry = btp.customization.author(entry)


#get rid of {} in authors list and re-assign to an actual dictionary with the ADS keys as the dictionary keys.
ppd_bibs = {}
for entry in ppd_bib_database.entries:
    entry['author'] = [x.replace("}","").replace("{","") for x in entry['author']]
    ppd_bibs[entry['ID']] = entry    


####################################################################################################
#                                  Exoplanet Atmosphere Molecules                                  # 
####################################################################################################         

#load in the bibtex database
with open(module_path.parents[0].joinpath('bibtex/exo_refs.bib'), 'r') as input:
    exo_bib_database = btp.bparser.BibTexParser(common_strings=True).parse_file(input) 
    #split the authors into a list 
    for entry in exo_bib_database.entries:
        entry = btp.customization.author(entry)


#get rid of {} in authors list and re-assign to an actual dictionary with the ADS keys as the dictionary keys.
exo_bibs = {}
for entry in exo_bib_database.entries:
    entry['author'] = [x.replace("}","").replace("{","") for x in entry['author']]
    exo_bibs[entry['ID']] = entry        


####################################################################################################
#                                    Interstellar Ice Molecules                                    # 
####################################################################################################         

#load in the bibtex database
with open(module_path.parents[0].joinpath('bibtex/ice_refs.bib'), 'r') as input:
    ice_bib_database = btp.bparser.BibTexParser(common_strings=True).parse_file(input) 
    #split the authors into a list 
    for entry in ice_bib_database.entries:
        entry = btp.customization.author(entry)


#get rid of {} in authors list and re-assign to an actual dictionary with the ADS keys as the dictionary keys.
ice_bibs = {}
for entry in ice_bib_database.entries:
    entry['author'] = [x.replace("}","").replace("{","") for x in entry['author']]
    ice_bibs[entry['ID']] = entry            