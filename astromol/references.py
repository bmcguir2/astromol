import bibtexparser as btp

#load in the bibtex database
with open('astromol.bib', 'r') as input:
    bib_database = btp.bparser.BibTexParser(common_strings=True).parse_file(input) 

#re-assign to an actual dictionary with the ADS keys as the dictionary keys.
bibs = {}
for entry in bib_database.entries:
    bibs[entry['ID']] = entry