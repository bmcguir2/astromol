from pathlib import Path
from tempfile import TemporaryDirectory

from astromol.census import CensusView
from astromol.database import Database
from astromol.latex import (
    rate_by_atoms_fits,
    rate_by_atoms_table_fragments,
    write_rate_by_atoms_table,
)


db = Database()
view_2021 = CensusView.for_census(db, "2021")
view_2026 = CensusView.for_census(db, "2026")

fits_2021 = rate_by_atoms_fits(view_2021)
values_2021 = {
    fit.label: (round(fit.slope, 2), round(fit.r_squared, 2), fit.onset_year)
    for fit in fits_2021
}
assert values_2021 == {
    "2": (0.68, 0.98, 1968),
    "3": (0.78, 0.98, 1968),
    "4": (0.53, 0.98, 1968),
    "5": (0.48, 0.95, 1971),
    "6": (0.40, 0.98, 1970),
    "7": (0.17, 0.88, 1973),
    "8": (0.31, 0.92, 1975),
    "9": (0.19, 0.89, 1974),
    "10": (0.21, 0.86, 2001),
    "11": (0.23, 0.90, 2004),
    "12": (0.16, 0.84, 2001),
    "13+": (0.30, 0.60, 2018),
    "Fullerenes": (0.09, 0.57, 2010),
}
assert "PAHs" not in values_2021

legacy_displayed_r_values_2021 = {
    fit.label: (round(fit.slope, 2), round(fit.r_value, 2), fit.onset_year)
    for fit in fits_2021
}
assert legacy_displayed_r_values_2021["2"] == (0.68, 0.99, 1968)
assert legacy_displayed_r_values_2021["12"] == (0.16, 0.92, 2001)
assert legacy_displayed_r_values_2021["13+"] == (0.30, 0.77, 2018)
assert legacy_displayed_r_values_2021["Fullerenes"] == (0.09, 0.75, 2010)

fragments_2021 = rate_by_atoms_table_fragments(view_2021)
assert set(fragments_2021) == {"rates_by_atoms_table.tex"}
content_2021 = fragments_2021["rates_by_atoms_table.tex"]
assert "least-squares linear fits using NumPy" in content_2021
assert "scipy.stats.linregress" not in content_2021
assert "$R^2$" in content_2021
assert "13+\t&\t0.30\t&\t0.60\t&\t2018" + r"\\" in content_2021
assert "PAHs" not in content_2021
assert content_2021.count(r"\\") == 14

with TemporaryDirectory() as tmp:
    output_dir = Path(tmp)
    assert write_rate_by_atoms_table(view_2021, output_dir) == fragments_2021
    assert (output_dir / "rates_by_atoms_table.tex").read_text() == content_2021

fits_2026 = rate_by_atoms_fits(view_2026)
values_2026 = {
    fit.label: (round(fit.slope, 2), round(fit.r_squared, 2), fit.onset_year)
    for fit in fits_2026
}
assert values_2026["13+"] == (0.93, 0.90, 2018)
assert values_2026["PAHs"] == (1.43, 0.82, 2021)
assert values_2026["Fullerenes"] == (0.05, 0.44, 2010)

content_2026 = rate_by_atoms_table_fragments(view_2026)["rates_by_atoms_table.tex"]
assert "13+\t&\t0.93\t&\t0.90\t&\t2018" + r"\\" in content_2026
assert "PAHs\t&\t1.43\t&\t0.82\t&\t2021" + r"\\" in content_2026

print("LaTeX rate-by-atoms table verification passed")
