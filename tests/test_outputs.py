from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from astromol import outputs


@dataclass
class _FakeDb:
    molecules: tuple[int, ...] = (1, 2, 3)
    detections: tuple[int, ...] = (1, 2, 3, 4)
    sources: tuple[int, ...] = (1, 2)
    telescopes: tuple[int, ...] = (1,)


@dataclass
class _FakeView:
    name: str


class _FakeFigureSpec:
    description = "Synthetic registry figure for output-bundle tests."

    def __init__(self, name: str, label: str) -> None:
        self.name = name
        self.label = label
        self.calls: list[tuple[str, str | Path]] = []

    def write(self, context, output_path: str | Path) -> Path:
        self.calls.append((context.view_choice, output_path))
        output_path = Path(output_path)
        output_path.write_text(f"figure:{context.view_choice}\n", encoding="utf-8")
        return output_path


class _FakeTableSpec:
    name = "demo_table"
    label = "Demo table"
    description = "Synthetic registry table for output-bundle tests."

    def write(self, view, output_dir: str | Path) -> list[Path]:
        output_path = Path(output_dir) / f"{view.name}_demo_table.tex"
        output_path.write_text("table\n", encoding="utf-8")
        return [output_path]


class _FakeSlideSpec:
    description = "Synthetic registry slide for output-bundle tests."
    report_label = "Demo slide layout report"
    report_description = "Synthetic layout report for output-bundle tests."

    def __init__(self, name: str, label: str, stem: str) -> None:
        self.name = name
        self.label = label
        self.stem = stem
        self.slide_calls: list[str] = []
        self.report_calls: list[str] = []

    def write_slide(self, context, output_dir: str | Path) -> Path:
        self.slide_calls.append(context.view_choice)
        output_path = Path(output_dir) / f"{self.stem}_{context.view_choice}.pptx"
        output_path.write_text("slide\n", encoding="utf-8")
        return output_path

    def write_report(self, context, output_dir: str | Path) -> Path:
        self.report_calls.append(context.view_choice)
        output_path = Path(output_dir) / f"{self.stem}_{context.view_choice}_layout.md"
        output_path.write_text("report\n", encoding="utf-8")
        return output_path


def test_generate_standard_outputs_uses_registry_specs(monkeypatch, tmp_path: Path) -> None:
    """Bundle generation should iterate the registry specs and write the index."""
    selected_view = _FakeView("selected")
    baseline_view = _FakeView("baseline")
    cumulative_figure_spec = _FakeFigureSpec("cumulative_detections", "Cumulative ISM/CSM detections")
    figure_spec = _FakeFigureSpec("demo_figure", "Demo figure")
    table_spec = _FakeTableSpec()
    ism_slide_spec = _FakeSlideSpec("ism_molecule_slide", "Demo ISM slide", "astro_molecules")
    ppd_slide_spec = _FakeSlideSpec("ppd_detection_slide", "Demo PPD slide", "ppd_molecules")

    monkeypatch.setattr(outputs, "Database", lambda: _FakeDb())
    monkeypatch.setattr(outputs, "view_from_choice", lambda db, choice: selected_view)
    monkeypatch.setattr(outputs.CensusView, "for_census", lambda db, choice: baseline_view)
    monkeypatch.setattr(outputs, "FIGURE_OUTPUTS", (cumulative_figure_spec, figure_spec))
    monkeypatch.setattr(outputs, "TABLE_OUTPUTS", (table_spec,))
    monkeypatch.setattr(outputs, "SLIDE_OUTPUTS", (ism_slide_spec, ppd_slide_spec))

    products = outputs.generate_standard_outputs(
        output_dir=tmp_path,
        view_choice="2026",
        formats=("png", "pdf"),
    )

    figure_paths = [
        tmp_path / "figures" / suffix / f"cumulative_detections.{suffix}"
        for suffix in ("png", "pdf")
    ]
    extra_figure_paths = [
        tmp_path / "figures" / suffix / f"demo_figure.{suffix}"
        for suffix in ("png", "pdf")
    ]
    table_path = tmp_path / "tables" / "selected_demo_table.tex"
    slide_paths = [
        tmp_path / "slides" / "astro_molecules_2026.pptx",
        tmp_path / "slides" / "ppd_molecules_2026.pptx",
    ]
    report_paths = [
        tmp_path / "slides" / "astro_molecules_2026_layout.md",
        tmp_path / "slides" / "ppd_molecules_2026_layout.md",
    ]
    zip_path = tmp_path / "astromol_latest_outputs.zip"
    index_path = tmp_path / "index.html"

    for path in (*figure_paths, *extra_figure_paths, table_path, *slide_paths, *report_paths, zip_path, index_path):
        assert path.exists()

    assert cumulative_figure_spec.calls == [("2026", figure_paths[0]), ("2026", figure_paths[1])]
    assert figure_spec.calls == [("2026", extra_figure_paths[0]), ("2026", extra_figure_paths[1])]
    assert ism_slide_spec.slide_calls == ["2026"]
    assert ism_slide_spec.report_calls == ["2026"]
    assert ppd_slide_spec.slide_calls == ["2026"]
    assert ppd_slide_spec.report_calls == ["2026"]

    product_paths = {product.path for product in products}
    assert product_paths == {
        figure_paths[0],
        figure_paths[1],
        extra_figure_paths[0],
        extra_figure_paths[1],
        table_path,
        slide_paths[0],
        slide_paths[1],
        report_paths[0],
        report_paths[1],
        zip_path,
    }

    index_html = index_path.read_text(encoding="utf-8")
    assert "Primary Downloads" in index_html
    assert "Presentation Decks" in index_html
    assert "Figures" in index_html
    assert "Cumulative ISM/CSM detections" in index_html
    assert "PNG Preview" in index_html
    assert "PDF File" in index_html
    assert "Additional Files" not in index_html
    assert "Layout Reports" not in index_html
    assert "LaTeX Tables" not in index_html
    assert "Demo table" not in index_html
    assert "tables/selected_demo_table.tex" not in index_html
    assert "The latest set of standard downloads generated from the selected census view." in index_html
    assert "registry-driven" not in index_html
    assert "Demo figure" in index_html
    assert "Demo ISM slide" in index_html
    assert "Demo PPD slide" in index_html
    assert "Demo figure" in index_html
    assert "Complete output bundle" in index_html
    assert '<span class="file">' not in index_html
    assert 'href="figures/png/demo_figure.png"' in index_html
    assert 'href="figures/pdf/demo_figure.pdf"' in index_html
    assert 'target="_blank"' in index_html
    assert 'rel="noopener noreferrer"' in index_html
    assert "figures/png/demo_figure.png" in index_html
    assert "figures/pdf/demo_figure.pdf" in index_html
    assert 'href="slides/astro_molecules_2026.pptx" target="_blank"' not in index_html
    assert "slides/astro_molecules_2026.pptx" in index_html
    assert "slides/ppd_molecules_2026.pptx" in index_html
    assert "slides/astro_molecules_2026_layout.md" not in index_html
    assert "slides/ppd_molecules_2026_layout.md" not in index_html
