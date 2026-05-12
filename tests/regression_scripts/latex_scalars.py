from pathlib import Path
from tempfile import TemporaryDirectory

from astromol.census import CensusView
from astromol.database import Database
from astromol.latex import scalar_fragments, write_scalar_fragments


db = Database()
view_2021 = CensusView.for_census(db, "2021")

expected = {
    "ndetects.tex": r"240\endinput",
    "nelems.tex": r"19\endinput",
    "nppds.tex": r"25\endinput",
    "nppdisos.tex": r"15\endinput",
    "nexgal.tex": r"73\endinput",
    "nexgalpercent.tex": r"30\endinput",
    "nexos.tex": r"9\endinput",
    "nices.tex": r"9\endinput",
    "radiopercent.tex": r"90\endinput",
    "nscopes.tex": r"46\endinput",
    "unsatpercent.tex": r"93\endinput",
    "satlist.tex": (
        r"\ce{CH3Cl}, \ce{CH4}, \ce{CH3OH}, \ce{CH3SH}, "
        r"\ce{CH3NH2}, \ce{CH3CH2OH}, \ce{CH3CH2SH}, "
        r"\ce{CH3OCH3}, \ce{CH3OCH2OH}, and \ce{HOCH2CH2OH}\endinput"
    ),
    "nsats.tex": r"10\endinput",
    "satpercent.tex": r"7\endinput",
    "sfr_rad_percent.tex": r"8\endinput",
    "dark_rad_percent.tex": r"25\endinput",
    "rate_since_1968.tex": r"3.9\endinput",
    "rate_since_2005.tex": r"6.0\endinput",
}

assert scalar_fragments(view_2021) == expected

with TemporaryDirectory() as tmp:
    output_dir = Path(tmp)
    assert write_scalar_fragments(view_2021, output_dir) == expected
    for filename, content in expected.items():
        assert (output_dir / filename).read_text() == content

print("LaTeX scalar verification passed")
