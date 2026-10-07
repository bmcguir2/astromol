from pathlib import Path
from tempfile import TemporaryDirectory

from baseline import load_production_baseline

from astromol.census import CensusView
from astromol.database import Database
from astromol.latex import (
    source_table_entries,
    source_table_fragments,
    write_source_table,
)


db = Database()
counts_current = load_production_baseline()["regression_counts"][
    "latex_source_table_current"
]
view_2021 = CensusView.for_census(db, "2021")
view_current = CensusView.current(db, end_year=load_production_baseline()["analysis_end_year"])

entries_2021 = source_table_entries(view_2021)
labels_2021 = {label for label, _ in entries_2021}
assert len(entries_2021) == 38
assert sum(count for _, count in entries_2021) == 332
assert entries_2021[:5] == [
    ("Sgr B2", 68),
    ("TMC-1", 57),
    ("IRC+10216", 55),
    ("Diffuse Cloud", 42),
    ("Orion", 24),
]
assert "Diffuse Cloud" in labels_2021
assert "DiffuseCloud" not in labels_2021
assert "Sgr B2 LOS" not in labels_2021

fragments_2021 = source_table_fragments(view_2021)
assert set(fragments_2021) == {"source_table.tex"}
content_2021 = fragments_2021["source_table.tex"]
assert r"\label{detects_by_source}" in content_2021
assert "Source\t&\t\\# \t&\tSource\t&\t\\# \\\\" in content_2021
assert "Sgr B2\t&\t68\t&\tL1527\t&\t2\t\\\\" in content_2021
assert "TMC-1\t&\t57\t&\tL1544\t&\t2\t\\\\" in content_2021
assert "Diffuse Cloud\t&\t42\t&\tNGC 7023\t&\t2\t\\\\" in content_2021
assert content_2021.count(r"\\") == 20

with TemporaryDirectory() as tmp:
    output_dir = Path(tmp)
    assert write_source_table(view_2021, output_dir) == fragments_2021
    assert (output_dir / "source_table.tex").read_text() == content_2021

entries_current = source_table_entries(view_current)
labels_current = {label for label, _ in entries_current}
assert len(entries_current) == counts_current["entries"]
assert sum(count for _, count in entries_current) == counts_current["credited_detections"]
assert entries_current[:5] == [
    tuple(entry) for entry in counts_current["top_entries"]
]
assert "Diffuse Cloud" in labels_current
assert "DiffuseCloud" not in labels_current
assert "Sgr B2 LOS" not in labels_current

content_current = source_table_fragments(view_current)["source_table.tex"]
top_source, top_source_count = counts_current["top_entries"][0]
assert f"{top_source}\t&\t{top_source_count}\t&" in content_current
split_at = (len(entries_current) + 1) // 2
left_entries = entries_current[:split_at]
right_entries = entries_current[split_at:]
for (left_label, left_count), (right_label, right_count) in zip(
    left_entries[:5],
    right_entries[:5],
):
    assert (
        f"{left_label}\t&\t{left_count}\t&\t"
        f"{right_label}\t&\t{right_count}\t\\\\"
    ) in content_current
assert content_current.count(r"\\") == split_at + 1

print("LaTeX source table verification passed")
