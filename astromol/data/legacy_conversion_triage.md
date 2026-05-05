# Legacy Conversion Triage

Generated from `scripts/convert_legacy_molecules.py`.

## Summary

- Converted molecules: 325
- Converted core detections: 325
- Reference mappings applied: 867
- Remaining issues: 170

## Remaining Issues by Kind

- `extra_context_detection_omitted`: 139
- `formula_normalized`: 2
- `missing_name_filled_from_formula`: 4
- `nested_isotopologues_omitted`: 15
- `unresolved_free_text_reference`: 10

## Reference Alias Mappings

- line 1264 `CN` [lab] Thomas & Dalby 1968 Can. J. Phys. 46, 2815 -> Thomson:1968:2815
- line 1348 `CO` [lab] Cord et al. 1968 Microwave Spectral Tables V5 -> Cord:1968:
- line 2471 `N2Hp` [observation] Thaddues & Turner 1975 ApJ 201, L25 -> Thaddeus:1975:L25
- line 3811 `SiC3` [lab] McCarthy et al. JCP 110, 1064 -> McCarthy:1999:10645
- line 4014 `HCCO` [lab] Oshima & Endo 1993 JMS 159, 458 -> Ohshima:1993:458
- line 4318 `HCOOH` [lab] Zukerman et al. 1971 ApJ 163, L41 -> Zuckerman:1971:L41
- line 4318 `HCOOH` [observation] Zukerman et al. 1971 ApJ 163, L41 -> Zuckerman:1971:L41
- line 4456 `SiH4` [lab] Goldhaber and Betz 1977 ApJ 279, L55 -> Goldhaber:1984:L55
- line 4456 `SiH4` [observation] Goldhaber and Betz 1977 ApJ 279, L55 -> Goldhaber:1984:L55
- line 4805 `NCCNHp` [lab] Gottlieb et al. 200 JCP 113, 1910 -> Gottlieb:2000:1910
- line 4821 `CH3Cl` [observation] Fayolle et al. 2017 Nature Astron. 1, 702 -> Fayolle:2017:703
- line 5286 `CH3CN` [lab] Cord et al. 1968 Microwave Spectral Tables V5 -> Cord:1968:
- line 5286 `CH3CN` [lab] Kessler et al. Phys Rev 79, 54 -> Kessler:1950:54
- line 5494 `HC3NHp` [lab] Lee & Amano 1987 ApJ 323 -> Lee:1987:L145
- line 5648 `HNCHCN` [lab] Zaleski et al. 2013 ApJ 765, L9 -> Zaleski:2013:L10
- line 5648 `HNCHCN` [observation] Zaleski et al. 2013 ApJ 765, L9 -> Zaleski:2013:L10
- line 5771 `HCSCCH` [lab] Brown:1982ur -> Brown:1982:1747
- line 5771 `HCSCCH` [lab] Crabtree:2016fj -> Crabtree:2016:124201
- line 6259 `CH2CHOH` [lab] Kaushik 1977 CPL 49, 90 -> Kaushik:1977:89
- line 6777 `HC6H` [lab] Haas etal. 1994 JMS 167, 176 -> Haas:1994:176
- line 6820 `CH2CCHCN` [lab] Bouche et al. 1973 J Mol Struct 18, 211 -> Bouchy:1973:211
- line 6868 `CH3CHNH` [lab] Loomis et al. 2013 ApJL 765, L10 -> Loomis:2013:L9
- line 6868 `CH3CHNH` [observation] Loomis et al. 2013 ApJL 765, L10 -> Loomis:2013:L9
- line 7185 `CH3OCH3` [lab] Kasai & Myers JCP 30, 1096 -> Kasai:1959:1096
- line 7222 `CH3CH2OH` [lab] Takano et al. 1986 JMS 26, 157 -> Takano:1968:157
- line 7222 `CH3CH2OH` [observation] Zukerman et al. 1975 ApJ 196, L99 -> Zuckerman:1975:L99
- line 7801 `C6H4` [lab] Brown:1986lp -> Brown:1986:1296
- line 7824 `C2H5NCO` [lab] Sakaizumi:1976uu -> Sakaizumi:1976:2908
- line 7824 `C2H5NCO` [lab] Heineking:1994op -> Heineking:1994:1177
- line 8127 `C5H6` [lab] Benson & Flygare 1970 J Am Chem Soc 92, 7523 -> Flygare:1970:7523
- line 8145 `NH2CH2CH2OH` [lab] Kaushik:1982ld -> Kaushik:1982:117
- line 8733 `C9H8` [observation] Cernicharo et al. 2021 A&AL 649, 15 -> Cernicharo:2021:L15a

## Ambiguous References

- None

## Unresolved Reference Keys

- None

## Unresolved Free-Text References

- 1x Additional work used in Belloche et al. 2019 A&A 628, A10 to be reported in Medvedev et al. in prep as of 9/16/2019.
  - line 6915 `NH2CONH2` formula `NH2CONH2` name `urea` [lab]
- 1x Dixon 1959 Can J. Phys. 37, 1171 and Klaus et al. 1997 A&A 322, L1
  - line 1786 `NH` formula `NH` name `imidogen radical` [lab]
- 1x Guélin et al. 1998 A&A 355, L1
  - line 5510 `C5N` formula `C5N` name `cyanobutadiynyl radical` [observation]
- 1x Jevons 1932 Phys Soc. pp 177-179
  - line 1237 `CH` formula `CH` name `methylidyne` [lab]
- 1x Kattija-Ari & Harmony et al. 1980 International Journal of Quantum Chemistry 18, 443
  - line 8103 `CH3COCH2OH` formula `CH3COCH2OH` name `hydroxyacetone` [lab]
- 1x Kruger et al. 2010 Ang. Chem. 23, 1644
  - line 4622 `HCCNC` formula `HCCNC` name `isocyanoacetylene` [lab]
- 1x Miller et al. 1962 JMS 8, 153
  - line 4381 `NH2CN` formula `NH2CN` name `cyanamide` [lab]
- 1x Nakimi et al. 1998 JMS 191, 176
  - line 2056 `TiO` formula `TiO` name `titanium monoxide` [lab]
- 1x Shinegari 1967 J Phys Soc Jpn 23, 404
  - line 3412 `H2CO` formula `H2CO` name `formaldehyde` [lab]
- 1x Steenbeckeliers 1968 Ann. Soc. Sci. Brux 82, 331
  - line 2532 `SO2` formula `SO2` name `sulfur dioxide` [lab]

## Deferred Categories

- `extra_context_detection_omitted`: 139
- `nested_isotopologues_omitted`: 15
- `missing_name_filled_from_formula`: 4
- `formula_normalized`: 2
