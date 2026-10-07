from pathlib import Path
from tempfile import TemporaryDirectory

from baseline import load_production_baseline

from astromol.census import CensusView
from astromol.database import Database
from astromol.latex import (
    facility_table_entries,
    facility_table_fragments,
    write_facility_table,
)


db = Database()
counts_current = load_production_baseline()["regression_counts"][
    "latex_facility_table_current"
]
view_2021 = CensusView.for_census(db, "2021")
view_current = CensusView.current(db, end_year=load_production_baseline()["analysis_end_year"])

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

entries_current = facility_table_entries(view_current)
assert len(entries_current) == counts_current["entries"]
assert sum(count for _, count in entries_current) == counts_current["credited_detections"]
assert entries_current[:5] == [
    tuple(entry) for entry in counts_current["top_entries"]
]

content_current = facility_table_fragments(view_current)["facilities_table.tex"]
top_facility, top_facility_count = counts_current["top_entries"][0]
assert f"{top_facility}\t&\t{top_facility_count}\t&" in content_current
split_at = (len(entries_current) + 1) // 2
assert content_current.count(r"\\") == split_at + 1

print("LaTeX facility table verification passed")
