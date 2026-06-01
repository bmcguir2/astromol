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
counts_2026 = load_production_baseline()["regression_counts"][
    "latex_source_table_2026"
]
view_2021 = CensusView.for_census(db, "2021")
view_2026 = CensusView.for_census(db, "2026")

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

entries_2026 = source_table_entries(view_2026)
labels_2026 = {label for label, _ in entries_2026}
assert len(entries_2026) == counts_2026["entries"]
assert sum(count for _, count in entries_2026) == counts_2026["credited_detections"]
assert entries_2026[:5] == [
    tuple(entry) for entry in counts_2026["top_entries"]
]
assert "Diffuse Cloud" in labels_2026
assert "DiffuseCloud" not in labels_2026
assert "Sgr B2 LOS" not in labels_2026

content_2026 = source_table_fragments(view_2026)["source_table.tex"]
top_source, top_source_count = counts_2026["top_entries"][0]
assert f"{top_source}\t&\t{top_source_count}\t&" in content_2026
split_at = (len(entries_2026) + 1) // 2
left_entries = entries_2026[:split_at]
right_entries = entries_2026[split_at:]
for (left_label, left_count), (right_label, right_count) in zip(
    left_entries[:5],
    right_entries[:5],
):
    assert (
        f"{left_label}\t&\t{left_count}\t&\t"
        f"{right_label}\t&\t{right_count}\t\\\\"
    ) in content_2026
assert content_2026.count(r"\\") == split_at + 1

print("LaTeX source table verification passed")
