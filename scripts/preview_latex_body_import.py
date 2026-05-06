"""Preview import of molecule prose from the 2021 census TeX source.

This script intentionally writes staging artifacts only. It does not modify
the production molecules.json file.
"""

from __future__ import annotations

from collections import Counter, defaultdict
import json
from pathlib import Path
import re

try:
    import bibtexparser
except ImportError as exc:  # pragma: no cover - mirrors runtime dependency
    raise SystemExit("This script requires bibtexparser.") from exc


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "astromol" / "data"
TEX_PATH = DATA / "2021_census_arxiv.tex"
MOLECULES_PATH = DATA / "molecules.json"
REFERENCES_PATH = DATA / "references.bib"
PREVIEW_PATH = DATA / "molecules.latex_body.preview.json"
IMPORT_PATH = DATA / "latex_body_import.preview.json"
REPORT_PATH = DATA / "latex_body_import_report.md"

MOLECULE_SECTIONS = {
    "Known Interstellar Molecules",
    "Tentative Detections",
    "Disputed Detections",
}

FORMULA_ALIASES = {
    "(NH2)2CO": "NH2CONH2",
    "(CH3)2CO": "CH3COCH3",
    "CH3CH2OCHO": "C2H5OCHO",
    "HO(CH2)2OH": "HOCH2CH2OH",
    "c-SiC2": "SiC2",
}

CITATION_ALIASES = {
    "Agundez:2018tm": "Agundez:2018:L22",
    "Brunken:2007fx": "Brunken:2007:L43",
    "Buhl:1973tp": "Buhl:1973:187",
    "Campbell:2016hl": "Campbell:2016:17",
    "Cord:1968kw": "Cord:1968:",
    "Decin:2018ju": "Decin:2018:113",
    "Godfrey:1997wv": "Godfrey:1997:405",
    "Gottlieb:1973vc": "Gottlieb:1973:181",
    "Gusten:2019cj": "Gusten:2019:357",
    "Hoeft:1965fh": "Hoeft:1965:1327",
    "Hollis:1991ua": "Hollis:1991:407",
    "Jevons:1932id": "Jevons:1932:177",
    "Kaminski:2013gk": "Kaminski:2013:A113",
    "Katija:1980pd": "Kattija-Ari:2009:443",
    "Kroto:1984hu": "Kroto:1984:993",
    "Krueger:2010gi": "Guarnieri:1992:39",
    "Kukolich:1965ey": "Kukolich:1965:1322",
    "Lee:2021pd": "Lee:2021:L11",
    "Lee:2021ud": "Lee:2021:L2a",
    "Morino:2000ff": "Morino:2000:367",
    "RodriguezAlmeida:2021ht": "Rodriguez-Almeida:2021:L11",
    "Torring:1968pd": "Torring:1968:777",
}


def unwrap_command(text: str, command: str) -> str:
    """Remove wrappers like \\update{...}, preserving their contents."""
    token = "\\" + command + "{"
    out = []
    i = 0
    while i < len(text):
        if text.startswith(token, i):
            i += len(token)
            depth = 1
            inner = []
            while i < len(text) and depth:
                ch = text[i]
                if ch == "{":
                    depth += 1
                    inner.append(ch)
                elif ch == "}":
                    depth -= 1
                    if depth:
                        inner.append(ch)
                else:
                    inner.append(ch)
                i += 1
            out.append(unwrap_command("".join(inner), command))
        else:
            out.append(text[i])
            i += 1
    return "".join(out)


def unwrap_common(text: str) -> str:
    for command in ("update", "emph", "textbf", "textsc", "mathrm"):
        text = unwrap_command(text, command)
    return text


def tex_to_plain(text: str) -> str:
    """Normalize just enough TeX for matching/reporting."""
    text = unwrap_common(text)
    replacements = {
        r"{\'u}": "u",
        r"{\"u}": "u",
        r"{\'i}": "i",
        r"{\i}": "i",
        r"{\ss}": "ss",
        r"\&": "&",
        r"~": " ",
    }
    for old, new in replacements.items():
        text = text.replace(old, new)
    text = re.sub(r"\\[a-zA-Z]+\{([^{}]*)\}", r"\1", text)
    text = text.replace("{", "").replace("}", "")
    text = re.sub(r"\s+", " ", text)
    return text.strip()


def normalize_formula(text: str) -> str:
    text = unwrap_common(text)
    text = re.sub(r"\\ce\{([^{}]*)\}", r"\1", text)
    text = text.replace("$", "")
    text = text.replace("\\emph{", "").replace("}", "")
    text = text.replace("{", "").replace("}", "")
    text = text.replace("\\", "")
    text = text.replace(" ", "")
    text = text.replace("^+", "+").replace("^-", "-")
    text = text.replace("^", "")
    return FORMULA_ALIASES.get(text, text)


def split_title_formula_name(title: str) -> tuple[str, str | None]:
    """Split a subsection title into formula-ish and common-name parts."""
    title = unwrap_common(title).strip()
    depth = 0
    in_math = False
    split_at = None
    for idx, ch in enumerate(title):
        if ch == "$":
            in_math = not in_math
        elif ch == "{":
            depth += 1
        elif ch == "}":
            depth = max(0, depth - 1)
        elif ch == "(" and depth == 0 and not in_math:
            split_at = idx
            break

    if split_at is None:
        return normalize_formula(title), None

    formula = normalize_formula(title[:split_at].strip())
    name_part = title[split_at:].strip()
    name = None
    if name_part.startswith("(") and name_part.endswith(")"):
        name = tex_to_plain(name_part[1:-1])
    return formula, name


def command_arg_from_line(line: str, command: str) -> str | None:
    token = "\\" + command
    pos = line.find(token)
    if pos == -1:
        return None
    brace = line.find("{", pos)
    if brace == -1:
        return None
    depth = 0
    out = []
    for ch in line[brace:]:
        if ch == "{":
            if depth:
                out.append(ch)
            depth += 1
        elif ch == "}":
            depth -= 1
            if depth == 0:
                return "".join(out)
            out.append(ch)
        elif depth:
            out.append(ch)
    return None


def section_title(line: str) -> str | None:
    arg = command_arg_from_line(line.strip(), "section")
    if arg is None:
        return None
    return tex_to_plain(arg)


def subsection_title(line: str) -> str | None:
    arg = command_arg_from_line(line.strip(), "subsection")
    if arg is None:
        return None
    return arg


def strip_standalone_labels(lines: list[str]) -> tuple[list[str], str | None]:
    label = None
    kept = []
    for line in lines:
        stripped = line.strip()
        match = re.fullmatch(r"\\label\{([^{}]+)\}", stripped)
        if match:
            if label is None:
                label = match.group(1)
            continue
        kept.append(line)
    return kept, label


def trim_blank_lines(lines: list[str]) -> list[str]:
    while lines and not lines[0].strip():
        lines.pop(0)
    while lines and not lines[-1].strip():
        lines.pop()
    return lines


def extract_tex_subsections(tex: str) -> list[dict]:
    body_tex = tex.split(r"\begin{thebibliography}", 1)[0]
    lines = body_tex.splitlines()
    headings = []
    current_section = None

    for idx, line in enumerate(lines):
        stripped = line.strip()
        if stripped.startswith(r"\section"):
            current_section = section_title(stripped)
            headings.append(("section", idx, current_section, stripped))
        elif stripped.startswith(r"\subsection"):
            title = subsection_title(stripped)
            headings.append(("subsection", idx, current_section, title))

    subsections = []
    for pos, (kind, idx, current_section, title) in enumerate(headings):
        if kind != "subsection" or current_section not in MOLECULE_SECTIONS:
            continue
        next_idx = len(lines)
        for next_kind, next_line, _, _ in headings[pos + 1 :]:
            if next_kind in {"section", "subsection"}:
                next_idx = next_line
                break
        raw_body_lines = lines[idx + 1 : next_idx]
        body_lines, label = strip_standalone_labels(raw_body_lines)
        body_lines = trim_blank_lines(body_lines)
        body = "\n".join(body_lines)
        body = unwrap_command(body, "update").strip()
        title_plain = tex_to_plain(title)
        formula, name = split_title_formula_name(title)
        subsections.append(
            {
                "tex_line": idx + 1,
                "section": current_section,
                "tex_label": label,
                "tex_title": title,
                "tex_title_plain": title_plain,
                "tex_formula": formula,
                "tex_name": name,
                "latex_body": body,
            }
        )
    return subsections


def unique_index(pairs: list[tuple[str, str]]) -> dict[str, str]:
    buckets = defaultdict(list)
    for key, value in pairs:
        if key:
            buckets[key].append(value)
    return {key: values[0] for key, values in buckets.items() if len(set(values)) == 1}


def molecule_indexes(molecules: list[dict]) -> dict[str, dict[str, str]]:
    by_suffix = {}
    table_pairs = []
    formula_pairs = []
    name_pairs = []
    for molecule in molecules:
        label = molecule["label"]
        suffix = label.removeprefix("mol:")
        by_suffix[suffix] = label
        table_pairs.append((molecule.get("table_formula"), label))
        formula_pairs.append((molecule.get("formula"), label))
        name_pairs.append(((molecule.get("name") or "").lower(), label))
    return {
        "suffix": by_suffix,
        "table_formula": unique_index(table_pairs),
        "formula": unique_index(formula_pairs),
        "name": unique_index(name_pairs),
    }


def map_subsection_to_molecule(entry: dict, indexes: dict[str, dict[str, str]]) -> dict:
    tex_label = entry.get("tex_label")
    tex_formula = entry.get("tex_formula")
    tex_name = (entry.get("tex_name") or "").lower()

    candidates = []
    if tex_label and tex_label in indexes["suffix"]:
        candidates.append(("tex_label", indexes["suffix"][tex_label]))
    if tex_formula in indexes["suffix"]:
        candidates.append(("formula_as_label", indexes["suffix"][tex_formula]))
    if tex_formula in indexes["table_formula"]:
        candidates.append(("table_formula", indexes["table_formula"][tex_formula]))
    if tex_formula in indexes["formula"]:
        candidates.append(("formula", indexes["formula"][tex_formula]))
    if tex_name in indexes["name"]:
        candidates.append(("name", indexes["name"][tex_name]))

    labels = list(dict.fromkeys(label for _, label in candidates))
    if len(labels) == 1:
        methods = [method for method, label in candidates if label == labels[0]]
        confidence = "high" if methods[0] != "name" else "medium"
        return {
            "status": "mapped",
            "confidence": confidence,
            "label": labels[0],
            "mapping_methods": methods,
        }
    if len(labels) > 1:
        return {
            "status": "ambiguous",
            "confidence": "low",
            "label": None,
            "mapping_methods": [method for method, _ in candidates],
            "candidate_labels": labels,
        }
    return {
        "status": "unmapped",
        "confidence": "none",
        "label": None,
        "mapping_methods": [],
    }


def norm_bibcode(value: str | None) -> str | None:
    if not value:
        return None
    return value.replace(r"\&", "&").replace("\\", "").strip()


def first_page(pages: str | None) -> str | None:
    if not pages:
        return None
    page = re.split(r"--|-|,", str(pages), 1)[0].strip()
    return page or None


def normalize_doi(value: str | None) -> str | None:
    if not value:
        return None
    return value.lower().replace("https://doi.org/", "").replace("http://doi.org/", "").strip()


def current_reference_indexes() -> dict[str, dict[str, str]]:
    with REFERENCES_PATH.open() as handle:
        parsed = bibtexparser.load(handle)

    ads = {}
    doi = {}
    tuple_buckets = defaultdict(list)
    for entry in parsed.entries:
        key = entry.get("ID")
        annotation = entry.get("annotation", "")
        match = re.search(r"ADS Bibcode:\s*([^}]+)", annotation)
        if match:
            ads[norm_bibcode(match.group(1))] = key
        if entry.get("doi"):
            doi[normalize_doi(entry.get("doi"))] = key
        tup = (
            str(entry.get("year", "")).strip(),
            str(entry.get("volume", "")).strip(),
            first_page(entry.get("pages") or entry.get("page")),
        )
        if all(tup):
            tuple_buckets[tup].append(key)

    tuples = {
        key: values[0]
        for key, values in tuple_buckets.items()
        if len(set(values)) == 1
    }
    direct = {entry.get("ID"): entry.get("ID") for entry in parsed.entries}
    return {"direct": direct, "ads": ads, "doi": doi, "tuple": tuples}


def parse_old_bibliography(tex: str) -> dict[str, dict]:
    if r"\begin{thebibliography}" not in tex:
        return {}
    bib = tex.split(r"\begin{thebibliography}", 1)[1]
    bib = bib.split(r"\end{thebibliography}", 1)[0]
    lines = bib.splitlines()
    starts = [idx for idx, line in enumerate(lines) if line.strip().startswith(r"\bibitem")]
    entries = {}

    for pos, start in enumerate(starts):
        end = starts[pos + 1] if pos + 1 < len(starts) else len(lines)
        block_lines = lines[start:end]
        header_lines = []
        body_start = 0
        key = None
        for offset, line in enumerate(block_lines):
            header_lines.append(line.strip())
            header = " ".join(header_lines)
            match = re.search(r"\]\{([^{}]+)\}", header) or re.search(r"\\bibitem\{([^{}]+)\}", header)
            if match:
                key = match.group(1)
                body_start = offset + 1
                break
        if not key:
            continue
        body = " ".join(line.strip() for line in block_lines[body_start:]).strip()
        doi_match = re.search(r"\\dodoi\{([^{}]+)\}", body)
        year_match = re.search(r"(\d{4})(?:\{\\natexlab\{[a-z]\}\})?,\s*([^,]+),\s*([^,]+),\s*([^,\s]+)", body)
        year = volume = page = None
        if year_match:
            year = year_match.group(1)
            volume = tex_to_plain(year_match.group(3))
            page = tex_to_plain(year_match.group(4))
        entries[key] = {
            "old_key": key,
            "body": body,
            "doi": normalize_doi(doi_match.group(1)) if doi_match else None,
            "year": year,
            "volume": volume,
            "first_page": first_page(page),
        }
    return entries


def build_citation_map(old_refs: dict[str, dict], current_refs: dict[str, dict[str, str]]) -> dict[str, dict]:
    citation_map = {}
    for old_key, old in old_refs.items():
        new_key = None
        method = None
        normalized_old_key = norm_bibcode(old_key)
        if old_key in CITATION_ALIASES:
            new_key = CITATION_ALIASES[old_key]
            method = "manual_alias"
        elif old_key in current_refs["direct"]:
            new_key = current_refs["direct"][old_key]
            method = "exact"
        elif normalized_old_key in current_refs["ads"]:
            new_key = current_refs["ads"][normalized_old_key]
            method = "ads_bibcode"
        elif old.get("doi") and old["doi"] in current_refs["doi"]:
            new_key = current_refs["doi"][old["doi"]]
            method = "doi"
        else:
            tup = (old.get("year"), old.get("volume"), old.get("first_page"))
            if all(tup) and tup in current_refs["tuple"]:
                new_key = current_refs["tuple"][tup]
                method = "year_volume_page"

        citation_map[old_key] = {
            "new_key": new_key,
            "method": method,
            "old_reference": old,
        }
    return citation_map


CITE_RE = re.compile(r"\\cite(?P<kind>t|p|alt)?\{(?P<keys>[^{}]+)\}")


def rewrite_citations(text: str, citation_map: dict[str, dict]) -> tuple[str, list[dict], list[str]]:
    replacements = []
    unresolved = []

    def repl(match: re.Match) -> str:
        kind = match.group("kind") or ""
        old_keys = [key.strip() for key in match.group("keys").split(",") if key.strip()]
        new_keys = []
        for old_key in old_keys:
            mapped = citation_map.get(old_key, {})
            new_key = mapped.get("new_key")
            if new_key:
                new_keys.append(new_key)
                if new_key != old_key:
                    replacements.append(
                        {
                            "old": old_key,
                            "new": new_key,
                            "method": mapped.get("method"),
                        }
                    )
            else:
                new_keys.append(old_key)
                unresolved.append(old_key)
        return "\\cite" + kind + "{" + ",".join(new_keys) + "}"

    rewritten = CITE_RE.sub(repl, text)
    return rewritten, replacements, sorted(set(unresolved))


def prefix_marker(entry: dict) -> str | None:
    formula = entry.get("tex_formula") or ""
    title = entry.get("tex_title_plain") or ""
    for prefix in ("n-", "i-", "c-", "l-", "E-", "Z-", "1-", "2-", "5-", "4-"):
        if formula.startswith(prefix) or title.startswith(prefix):
            return prefix
    return None


def make_report(import_rows: list[dict], citation_map: dict[str, dict]) -> str:
    mapped_rows = [row for row in import_rows if row["status"] == "mapped"]
    problem_rows = [row for row in import_rows if row["status"] != "mapped"]
    citation_counts = Counter(
        row["status"] for row in citation_map.values() for row["status"] in []
    )
    mapped_cites = [old for old, row in citation_map.items() if row.get("new_key")]
    unresolved_cites = [old for old, row in citation_map.items() if not row.get("new_key")]
    used_unresolved = Counter(
        key
        for row in import_rows
        for key in row.get("unresolved_citations", [])
    )

    lines = []
    lines.append("# LaTeX Body Import Preview")
    lines.append("")
    lines.append("Generated by `scripts/preview_latex_body_import.py`.")
    lines.append("")
    lines.append("## Summary")
    lines.append("")
    lines.append(f"- Molecule subsections extracted: {len(import_rows)}")
    lines.append(f"- Mapped to modern molecules: {len(mapped_rows)}")
    lines.append(f"- Unmapped or ambiguous molecule sections: {len(problem_rows)}")
    lines.append(f"- Old bibliography entries parsed: {len(citation_map)}")
    lines.append(f"- Old citation keys mapped to current keys: {len(mapped_cites)}")
    lines.append(f"- Old citation keys not mapped: {len(unresolved_cites)}")
    lines.append("")
    lines.append("Preview artifacts:")
    lines.append(f"- `{PREVIEW_PATH.relative_to(ROOT)}`: full molecule JSON preview with imported `latex_body` values.")
    lines.append(f"- `{IMPORT_PATH.relative_to(ROOT)}`: per-subsection mapping, citation replacements, and body preview data.")
    lines.append("")

    lines.append("## Molecule Mapping Issues")
    lines.append("")
    if not problem_rows:
        lines.append("No molecule mapping issues found.")
    else:
        for row in problem_rows:
            lines.append(
                f"- line {row['tex_line']} `{row.get('tex_label')}` "
                f"`{row['tex_title_plain']}` -> {row['status']}"
            )
            if row.get("candidate_labels"):
                lines.append(f"  - candidates: {', '.join(row['candidate_labels'])}")
    lines.append("")

    lines.append("## Structural Prefix Review")
    lines.append("")
    prefix_rows = [row for row in import_rows if row.get("prefix_marker")]
    if not prefix_rows:
        lines.append("No structural prefixes found.")
    else:
        for row in prefix_rows:
            lines.append(
                f"- line {row['tex_line']} `{row['tex_title_plain']}` "
                f"uses `{row['prefix_marker']}` -> `{row.get('label')}` "
                f"({row['status']}, {row['confidence']})"
            )
            if row.get("modern_table_formula"):
                lines.append(f"  - modern table_formula: `{row['modern_table_formula']}`")
    lines.append("")

    lines.append("## Custom Label And Title Mappings")
    lines.append("")
    custom_rows = [
        row for row in import_rows
        if row["status"] == "mapped"
        and (
            row.get("tex_label") != row.get("label", "").removeprefix("mol:")
            or row.get("tex_formula") != row.get("modern_table_formula")
        )
    ]
    if not custom_rows:
        lines.append("No custom label/title mappings found.")
    else:
        for row in custom_rows:
            lines.append(
                f"- line {row['tex_line']} tex label `{row.get('tex_label')}` / "
                f"title formula `{row.get('tex_formula')}` -> `{row.get('label')}`"
            )
            if row.get("modern_table_formula"):
                lines.append(f"  - modern table_formula: `{row['modern_table_formula']}`")
    lines.append("")

    lines.append("## Mapped Non-Known Sections")
    lines.append("")
    non_known_rows = [
        row for row in import_rows
        if row["status"] == "mapped" and row["section"] != "Known Interstellar Molecules"
    ]
    if not non_known_rows:
        lines.append("No non-known sections mapped to modern molecules.")
    else:
        lines.append(
            "These imported bodies came from old tentative/disputed sections and may need prose updates."
        )
        lines.append("")
        for row in non_known_rows:
            lines.append(
                f"- line {row['tex_line']} `{row['tex_title_plain']}` "
                f"from `{row['section']}` -> `{row.get('label')}`"
            )
    lines.append("")

    lines.append("## Citation Replacements Used In Imported Bodies")
    lines.append("")
    replacement_counts = Counter(
        (replacement["old"], replacement["new"], replacement["method"])
        for row in import_rows
        for replacement in row.get("citation_replacements", [])
    )
    if not replacement_counts:
        lines.append("No citation replacements were made.")
    else:
        for (old, new, method), count in replacement_counts.most_common():
            lines.append(f"- {count}x `{old}` -> `{new}` ({method})")
    lines.append("")

    lines.append("## Unresolved Citations Used In Imported Bodies")
    lines.append("")
    if not used_unresolved:
        lines.append("No unresolved citations are used in mapped molecule bodies.")
    else:
        for key, count in used_unresolved.most_common():
            lines.append(f"- {count}x `{key}`")
    lines.append("")

    lines.append("## Mapped Body Previews")
    lines.append("")
    for row in mapped_rows:
        preview = " ".join(row["latex_body"].split())
        if len(preview) > 360:
            preview = preview[:357] + "..."
        lines.append(
            f"### line {row['tex_line']} `{row['tex_title_plain']}` -> `{row['label']}`"
        )
        lines.append("")
        lines.append(f"- confidence: `{row['confidence']}`")
        lines.append(f"- methods: `{', '.join(row['mapping_methods'])}`")
        if row.get("unresolved_citations"):
            lines.append(f"- unresolved citations: `{', '.join(row['unresolved_citations'])}`")
        lines.append("")
        lines.append(preview)
        lines.append("")

    return "\n".join(lines).rstrip() + "\n"


def main() -> None:
    tex = TEX_PATH.read_text()
    molecules = json.loads(MOLECULES_PATH.read_text())
    molecules_by_label = {molecule["label"]: molecule for molecule in molecules}

    subsections = extract_tex_subsections(tex)
    indexes = molecule_indexes(molecules)
    current_refs = current_reference_indexes()
    old_refs = parse_old_bibliography(tex)
    citation_map = build_citation_map(old_refs, current_refs)

    import_rows = []
    preview_molecules = json.loads(json.dumps(molecules))
    preview_by_label = {molecule["label"]: molecule for molecule in preview_molecules}

    for subsection in subsections:
        mapping = map_subsection_to_molecule(subsection, indexes)
        rewritten_body, replacements, unresolved = rewrite_citations(
            subsection["latex_body"], citation_map
        )
        row = {
            **subsection,
            **mapping,
            "latex_body": rewritten_body,
            "citation_replacements": replacements,
            "unresolved_citations": unresolved,
            "prefix_marker": prefix_marker(subsection),
        }
        if mapping["label"] in molecules_by_label:
            row["modern_table_formula"] = molecules_by_label[mapping["label"]].get("table_formula")
            row["existing_latex_body"] = molecules_by_label[mapping["label"]].get("latex_body")
            preview_by_label[mapping["label"]]["latex_body"] = rewritten_body
        import_rows.append(row)

    PREVIEW_PATH.write_text(json.dumps(preview_molecules, indent=2) + "\n")
    IMPORT_PATH.write_text(json.dumps(import_rows, indent=2) + "\n")
    REPORT_PATH.write_text(make_report(import_rows, citation_map))

    print(f"Wrote {PREVIEW_PATH.relative_to(ROOT)}")
    print(f"Wrote {IMPORT_PATH.relative_to(ROOT)}")
    print(f"Wrote {REPORT_PATH.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
