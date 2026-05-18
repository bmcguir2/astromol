from __future__ import annotations

import matplotlib.image as mpimg

from astromol.census import CensusView
from astromol.database import Database
from astromol.figures.style import FIGURE_DPI, FIGURE_SIZE
from astromol.registry import (
    FIGURE_OUTPUTS,
    SLIDE_OUTPUTS,
    TABLE_OUTPUTS,
    OutputContext,
    TableOutputSpec,
    figure_output_names,
    slide_output_names,
    table_output_names,
)


def test_registry_names_are_unique() -> None:
    """Registry names are stable public identifiers for generated outputs."""
    for names in (figure_output_names(), table_output_names(), slide_output_names()):
        assert len(names) == len(set(names))


def test_standard_registry_inventory() -> None:
    """The registry exposes the current standard output inventory."""
    assert len(FIGURE_OUTPUTS) == 22
    assert len(TABLE_OUTPUTS) == 9
    assert len(SLIDE_OUTPUTS) == 2
    assert "cumulative_detections" in figure_output_names()
    assert "ism_tables" in table_output_names()
    assert "ppd_detection_slide" in slide_output_names()


def test_registry_figure_spec_writes_output(tmp_path) -> None:
    """A registry figure spec can build data and write a file."""
    db = Database()
    context = OutputContext(
        view=CensusView.for_census(db, "2021"),
        view_choice="2021",
        baseline_view=CensusView.for_census(db, "2021"),
    )
    spec = next(item for item in FIGURE_OUTPUTS if item.name == "cumulative_detections")
    output_path = tmp_path / "cumulative_detections.png"

    written_path = spec.write(context, output_path)

    assert written_path == output_path
    assert output_path.exists()
    assert output_path.stat().st_size > 0


def test_registry_figure_png_uses_publication_dpi(tmp_path) -> None:
    """Registry PNG figures are written at the shared publication DPI."""
    db = Database()
    context = OutputContext(
        view=CensusView.for_census(db, "2021"),
        view_choice="2021",
        baseline_view=CensusView.for_census(db, "2021"),
    )
    spec = next(item for item in FIGURE_OUTPUTS if item.name == "cumulative_detections")
    output_path = tmp_path / "cumulative_detections.png"

    spec.write(context, output_path)

    image = mpimg.imread(output_path)
    expected_width = int(FIGURE_SIZE[0] * FIGURE_DPI)
    expected_height = int(FIGURE_SIZE[1] * FIGURE_DPI)
    assert image.shape[:2] == (expected_height, expected_width)


def test_slide_registry_filenames_use_view_choice() -> None:
    """Slide filenames are derived from registry metadata."""
    db = Database()
    context = OutputContext(view=CensusView.for_census(db, "2026"), view_choice="2026")

    filenames = {spec.name: spec.filename(context) for spec in SLIDE_OUTPUTS}
    report_filenames = {spec.name: spec.report_filename(context) for spec in SLIDE_OUTPUTS}

    assert filenames["ism_molecule_slide"] == "astro_molecules_2026.pptx"
    assert filenames["ppd_detection_slide"] == "ppd_molecules_2026.pptx"
    assert report_filenames["ism_molecule_slide"] == "astro_molecules_2026_layout.md"
    assert report_filenames["ppd_detection_slide"] == "ppd_molecules_2026_layout.md"


def test_table_registry_returns_generated_paths_on_rerun(tmp_path) -> None:
    """Table specs report generated files even when overwriting existing output."""

    def write_demo_table(_view, output_dir):
        output_path = output_dir / "demo_table.tex"
        output_path.write_text("demo\n", encoding="utf-8")
        return {"demo_table.tex": "demo\n"}

    spec = TableOutputSpec(
        "demo_table",
        "Demo table",
        "Synthetic table for registry path tests.",
        write_demo_table,
    )
    db = Database()
    view = CensusView.for_census(db, "2021")

    assert spec.write(view, tmp_path) == [tmp_path / "demo_table.tex"]
    assert spec.write(view, tmp_path) == [tmp_path / "demo_table.tex"]
