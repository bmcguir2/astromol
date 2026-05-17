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
    name = "demo_figure"
    label = "Demo figure"
    description = "Synthetic registry figure for output-bundle tests."

    def __init__(self) -> None:
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
    name = "demo_slide"
    label = "Demo slide"
    description = "Synthetic registry slide for output-bundle tests."

    def __init__(self) -> None:
        self.slide_calls: list[str] = []
        self.report_calls: list[str] = []

    def write_slide(self, context, output_dir: str | Path) -> Path:
        self.slide_calls.append(context.view_choice)
        output_path = Path(output_dir) / f"demo_{context.view_choice}.pptx"
        output_path.write_text("slide\n", encoding="utf-8")
        return output_path

    def write_report(self, context, output_dir: str | Path) -> Path:
        self.report_calls.append(context.view_choice)
        output_path = Path(output_dir) / f"demo_{context.view_choice}_layout.md"
        output_path.write_text("report\n", encoding="utf-8")
        return output_path

    report_label = "Demo slide layout report"
    report_description = "Synthetic layout report for output-bundle tests."


def test_generate_standard_outputs_uses_registry_specs(monkeypatch, tmp_path: Path) -> None:
    """Bundle generation should iterate the registry specs and write the index."""
    selected_view = _FakeView("selected")
    baseline_view = _FakeView("baseline")
    figure_spec = _FakeFigureSpec()
    table_spec = _FakeTableSpec()
    slide_spec = _FakeSlideSpec()

    monkeypatch.setattr(outputs, "Database", lambda: _FakeDb())
    monkeypatch.setattr(outputs, "view_from_choice", lambda db, choice: selected_view)
    monkeypatch.setattr(outputs.CensusView, "for_census", lambda db, choice: baseline_view)
    monkeypatch.setattr(outputs, "FIGURE_OUTPUTS", (figure_spec,))
    monkeypatch.setattr(outputs, "TABLE_OUTPUTS", (table_spec,))
    monkeypatch.setattr(outputs, "SLIDE_OUTPUTS", (slide_spec,))

    products = outputs.generate_standard_outputs(
        output_dir=tmp_path,
        view_choice="2026",
        formats=("png", "pdf"),
    )

    figure_paths = [tmp_path / "figures" / suffix / f"demo_figure.{suffix}" for suffix in ("png", "pdf")]
    table_path = tmp_path / "tables" / "selected_demo_table.tex"
    slide_path = tmp_path / "slides" / "demo_2026.pptx"
    report_path = tmp_path / "slides" / "demo_2026_layout.md"
    zip_path = tmp_path / "astromol_latest_outputs.zip"
    index_path = tmp_path / "index.html"

    for path in (*figure_paths, table_path, slide_path, report_path, zip_path, index_path):
        assert path.exists()

    assert figure_spec.calls == [("2026", figure_paths[0]), ("2026", figure_paths[1])]
    assert slide_spec.slide_calls == ["2026"]
    assert slide_spec.report_calls == ["2026"]

    product_paths = {product.path for product in products}
    assert product_paths == {figure_paths[0], figure_paths[1], table_path, slide_path, report_path, zip_path}

    index_html = index_path.read_text(encoding="utf-8")
    assert "Primary Downloads" in index_html
    assert "Full Inventory" in index_html
    assert "Slide Decks" in index_html
    assert "PNG Figures" in index_html
    assert "LaTeX Tables" in index_html
    assert "Demo figure" in index_html
    assert "Demo slide" in index_html
    assert "Demo table" in index_html
    assert "Complete output bundle" in index_html
    assert "figures/png/demo_figure.png" in index_html
    assert "slides/demo_2026.pptx" in index_html
    assert "tables/selected_demo_table.tex" in index_html
