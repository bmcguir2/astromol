from pathlib import Path
from tempfile import TemporaryDirectory

from astromol.census import CensusView
from astromol.database import Database
from astromol.latex import (
    facility_table_entries,
    facility_table_fragments,
    write_facility_table,
)


db = Database()
view_2021 = CensusView.for_census(db, "2021")
view_2026 = CensusView.for_census(db, "2026")

entries_2021 = facility_table_entries(view_2021)
assert len(entries_2021) == 46
assert sum(count for _, count in entries_2021) == 304
assert entries_2021[:5] == [
    ("IRAM 30-m", 64),
    ("NRAO 36-ft", 33),
    ("GBT 100-m", 28),
    ("NRAO/ARO 12-m", 27),
    ("Yebes 40-m", 19),
]
assert ("VLT", 0) not in entries_2021

fragments_2021 = facility_table_fragments(view_2021)
assert set(fragments_2021) == {"facilities_table.tex"}
content_2021 = fragments_2021["facilities_table.tex"]
assert r"\label{detects_by_scope}" in content_2021
assert "Facility\t&\t\\# \t&\tFacility\t&\t\\# \\\\" in content_2021
assert "IRAM 30-m\t&\t64\t&\tHat Creek 20-ft\t&\t2\t\\\\" in content_2021
assert "NRAO 36-ft\t&\t33\t&\tIRTF\t&\t2\t\\\\" in content_2021
assert content_2021.count(r"\\") == 24

with TemporaryDirectory() as tmp:
    output_dir = Path(tmp)
    assert write_facility_table(view_2021, output_dir) == fragments_2021
    assert (output_dir / "facilities_table.tex").read_text() == content_2021

entries_2026 = facility_table_entries(view_2026)
assert len(entries_2026) == 47
assert sum(count for _, count in entries_2026) == 415
assert entries_2026[:5] == [
    ("IRAM 30-m", 87),
    ("Yebes 40-m", 86),
    ("GBT 100-m", 36),
    ("NRAO 36-ft", 33),
    ("NRAO/ARO 12-m", 29),
]

content_2026 = facility_table_fragments(view_2026)["facilities_table.tex"]
assert "IRAM 30-m\t&\t87\t&\tIRTF\t&\t2\t\\\\" in content_2026
assert "Yebes 40-m\t&\t86\t&\tMWO 4.9-m\t&\t2\t\\\\" in content_2026
assert content_2026.count(r"\\") == 25

print("LaTeX facility table verification passed")
