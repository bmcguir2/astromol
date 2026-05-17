from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from types import SimpleNamespace

from astromol import cli


@dataclass
class _FakeFigureSpec:
    name: str = "demo_figure"
    label: str = "Demo figure"
    description: str = "Synthetic figure for CLI tests."

    def write(self, context, output_path: str | Path) -> Path:
        output_path = Path(output_path)
        output_path.write_text(f"figure:{context.view_choice}\n", encoding="utf-8")
        return output_path


@dataclass
class _FakeTableSpec:
    name: str = "demo_table"
    label: str = "Demo table"
    description: str = "Synthetic table for CLI tests."

    def write(self, view, output_dir: str | Path) -> list[Path]:
        output_path = Path(output_dir) / f"{view.name}_{self.name}.tex"
        output_path.write_text("table\n", encoding="utf-8")
        return [output_path]


@dataclass
class _FakeSlideSpec:
    name: str = "demo_slide"
    label: str = "Demo slide"
    description: str = "Synthetic slide for CLI tests."
    report_label: str = "Demo slide report"
    report_description: str = "Synthetic report for CLI tests."

    def write_slide(self, context, output_dir: str | Path) -> Path:
        output_path = Path(output_dir) / f"{self.name}_{context.view_choice}.pptx"
        output_path.write_text("slide\n", encoding="utf-8")
        return output_path

    def write_report(self, context, output_dir: str | Path) -> Path:
        output_path = Path(output_dir) / f"{self.name}_{context.view_choice}_layout.md"
        output_path.write_text("report\n", encoding="utf-8")
        return output_path


def _fake_context(view_choice: str):
    return SimpleNamespace(
        view=SimpleNamespace(name="selected"),
        view_choice=view_choice,
        baseline_view=SimpleNamespace(name="baseline"),
    )


def test_cli_list_prints_registry_names(monkeypatch, capsys) -> None:
    """The public CLI lists names from the registry specs."""
    monkeypatch.setattr(cli, "FIGURE_OUTPUTS", (_FakeFigureSpec(),))
    monkeypatch.setattr(cli, "TABLE_OUTPUTS", (_FakeTableSpec(),))
    monkeypatch.setattr(cli, "SLIDE_OUTPUTS", (_FakeSlideSpec(),))

    assert cli.main(["list"]) == 0

    output = capsys.readouterr().out
    assert "Figures:" in output
    assert "demo_figure" in output
    assert "Tables:" in output
    assert "demo_table" in output
    assert "Slides:" in output
    assert "demo_slide" in output


def test_cli_generates_one_figure(monkeypatch, tmp_path: Path, capsys) -> None:
    """The figure subcommand dispatches one registry figure by name."""
    output_path = tmp_path / "figure.pdf"
    monkeypatch.setattr(cli, "FIGURE_OUTPUTS", (_FakeFigureSpec(),))
    monkeypatch.setattr(cli, "_build_context", _fake_context)

    assert cli.main(["figure", "demo_figure", "--view", "2026", "--output", str(output_path)]) == 0

    assert output_path.read_text(encoding="utf-8") == "figure:2026\n"
    assert str(output_path) in capsys.readouterr().out


def test_cli_generates_one_table_group(monkeypatch, tmp_path: Path, capsys) -> None:
    """The table subcommand dispatches one registry table group by name."""
    monkeypatch.setattr(cli, "TABLE_OUTPUTS", (_FakeTableSpec(),))
    monkeypatch.setattr(cli, "_build_context", _fake_context)

    assert cli.main(["table", "demo_table", "--view", "2026", "--output-dir", str(tmp_path)]) == 0

    output_path = tmp_path / "selected_demo_table.tex"
    assert output_path.exists()
    assert str(output_path) in capsys.readouterr().out


def test_cli_generates_one_slide_and_report(monkeypatch, tmp_path: Path, capsys) -> None:
    """The slide subcommand can include the optional layout report."""
    monkeypatch.setattr(cli, "SLIDE_OUTPUTS", (_FakeSlideSpec(),))
    monkeypatch.setattr(cli, "_build_context", _fake_context)

    assert cli.main(["slide", "demo_slide", "--view", "2026", "--output-dir", str(tmp_path), "--report"]) == 0

    slide_path = tmp_path / "demo_slide_2026.pptx"
    report_path = tmp_path / "demo_slide_2026_layout.md"
    assert slide_path.exists()
    assert report_path.exists()
    output = capsys.readouterr().out
    assert str(slide_path) in output
    assert str(report_path) in output


def test_cli_outputs_delegates_to_bundle_generator(monkeypatch, tmp_path: Path, capsys) -> None:
    """The outputs subcommand preserves the existing full-bundle generator."""
    calls = []

    def fake_generate_standard_outputs(*args, **kwargs):
        calls.append((args, kwargs))
        return [SimpleNamespace(path=tmp_path / "one.txt"), SimpleNamespace(path=tmp_path / "two.txt")]

    monkeypatch.setattr(cli, "generate_standard_outputs", fake_generate_standard_outputs)

    assert (
        cli.main(
            [
                "outputs",
                "--output-dir",
                str(tmp_path),
                "--view",
                "2026",
                "--formats",
                "png",
                "--skip-slides",
            ]
        )
        == 0
    )

    assert calls == [
        (
            (tmp_path,),
            {
                "view_choice": "2026",
                "formats": ("png",),
                "clean": True,
                "include_tables": True,
                "include_slides": False,
            },
        )
    ]
    assert f"Generated 2 product(s) in {tmp_path}" in capsys.readouterr().out
