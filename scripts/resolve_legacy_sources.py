"""Stage candidate source records for legacy extra-context detections.

This script does not modify ``sources.json`` or ``detections.preview.json``.
It reads source gaps from ``legacy_conversion_issues.json``, resolves unique
legacy source strings where possible, and writes review artifacts:

* ``sources_additional.preview.json``: new candidate Source-like records.
* ``sources_additional_report.md``: human-readable review report.
* ``sources_additional_report.json``: machine-readable summary used by triage.
"""

from __future__ import annotations

import argparse
import json
import re
import time
import urllib.error
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from collections import Counter, defaultdict
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = REPO_ROOT / "astromol" / "data"
ISSUES_PATH = DATA_DIR / "legacy_conversion_issues.json"
SOURCES_PATH = DATA_DIR / "sources.json"

DEFAULT_ADDITIONAL_OUT = DATA_DIR / "sources_additional.preview.json"
DEFAULT_REPORT_MD_OUT = DATA_DIR / "sources_additional_report.md"
DEFAULT_REPORT_JSON_OUT = DATA_DIR / "sources_additional_report.json"

SESAME_ENDPOINT = "https://cds.unistra.fr/cgi-bin/nph-sesame/-oxp"

CONTEXT_SOURCE_TYPES = {
    "exgal": "External Galaxy",
    "exo": "Exoplanet",
    "ppd": "YSO",
}

QUERY_OVERRIDES = {
    "Cloverleaf LOS": ["Cloverleaf Quasar", "H 1413+117"],
    "HD 209458b": ["HD 209458 b", "HD 209458"],
    "PKS 1830-211 LOS": ["PKS 1830-211"],
    "QSO Mrk 231": ["Mrk 231"],
}

SOURCE_DERIVATIONS = {
    "Cloverleaf LOS": {
        "query": "Cloverleaf Quasar",
        "proposed_nick": "CloverleafLOS",
        "proposed_name": "Cloverleaf LOS",
        "reason": (
            "User-approved line-of-sight source: preserve the LOS qualifier "
            "while using the resolved Cloverleaf quasar coordinates."
        ),
    },
    "PKS 1830-211 LOS": {
        "source_nick": "PKS1830",
        "proposed_nick": "PKS1830LOS",
        "proposed_name": "PKS 1830-211 LOS",
        "reason": (
            "User-approved line-of-sight source: preserve the LOS qualifier "
            "while copying coordinates and link metadata from PKS1830."
        ),
    },
}

SOURCE_APPROVALS = {
    "Arp 220": {
        "reason": "User-approved high-confidence Sesame source match.",
    },
    "GG Tau": {
        "reason": "User-approved PPD source staging as source type YSO.",
        "note": "PPD source staged as YSO by user approval.",
    },
    "HD 209458b": {
        "reason": "User-approved high-confidence Sesame source match.",
    },
    "IC 342": {
        "reason": "User-approved high-confidence Sesame source match.",
    },
    "IRAS 08572+3915": {
        "reason": "User-approved high-confidence Sesame source match.",
    },
    "LMC": {
        "reason": "User-approved high-confidence Sesame source match.",
    },
    "M33": {
        "reason": "User-approved high-confidence Sesame source match.",
    },
    "M82": {
        "reason": "User-approved high-confidence Sesame source match.",
    },
    "Maffei 2": {
        "reason": "User-approved high-confidence Sesame source match.",
    },
    "NGC 1068": {
        "reason": "User-approved high-confidence Sesame source match.",
    },
    "NGC 253": {
        "reason": "User-approved high-confidence Sesame source match.",
    },
    "NGC 4418": {
        "reason": "User-approved high-confidence Sesame source match.",
    },
    "NGC 4945": {
        "reason": "User-approved high-confidence Sesame source match.",
    },
    "NGC 5128": {
        "reason": "User-approved high-confidence Sesame source match.",
    },
    "NGC 6946": {
        "reason": "User-approved high-confidence Sesame source match.",
    },
    "QSO Mrk 231": {
        "reason": "User-approved QSO-prefix source staging as Mrk 231.",
    },
    "SMC": {
        "reason": "User-approved high-confidence Sesame source match.",
    },
    "SMP LMC 11": {
        "reason": (
            "User-approved extragalactic source staging despite resolver "
            "matching an internal LMC substructure."
        ),
    },
}

SOURCE_TYPE_OVERRIDES = {
    "SMP LMC 11": "External Galaxy",
}


def main():
    args = parse_args()
    issues = load_json(args.issues)
    sources = load_json(args.sources)
    source_index = build_source_index(sources)
    usages = collect_source_usages(issues)

    records = []
    additional_sources = []
    used_nicks = {source["nick"] for source in sources}

    for legacy_source in sorted(usages):
        record = resolve_source_record(
            legacy_source,
            usages[legacy_source],
            source_index,
            used_nicks,
            skip_network=args.skip_network,
            timeout=args.timeout,
            sleep=args.sleep,
        )
        records.append(record)
        if record["action"] == "add_source":
            additional_sources.append(record_to_source_entry(record))
            used_nicks.add(record["proposed_nick"])
        if not args.skip_network and args.sleep:
            time.sleep(args.sleep)

    report = build_report(records, usages)
    write_json(args.additional_out, additional_sources)
    write_json(args.report_json_out, report)
    args.report_md_out.write_text(write_markdown_report(report), encoding="utf-8")

    print(
        "Reviewed "
        f"{report['summary']['unique_source_count']} unique sources across "
        f"{report['summary']['detection_count']} detections and "
        f"{report['summary']['source_assignment_count']} source assignments."
    )
    print(f"Wrote {len(additional_sources)} candidate sources to {args.additional_out}")
    print(f"Wrote source report to {args.report_md_out}")
    print(f"Wrote machine-readable source report to {args.report_json_out}")


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--issues", type=Path, default=ISSUES_PATH)
    parser.add_argument("--sources", type=Path, default=SOURCES_PATH)
    parser.add_argument("--additional-out", type=Path, default=DEFAULT_ADDITIONAL_OUT)
    parser.add_argument("--report-md-out", type=Path, default=DEFAULT_REPORT_MD_OUT)
    parser.add_argument("--report-json-out", type=Path, default=DEFAULT_REPORT_JSON_OUT)
    parser.add_argument(
        "--skip-network",
        action="store_true",
        help="Do not query Sesame; only report existing local source matches.",
    )
    parser.add_argument("--timeout", type=float, default=20.0)
    parser.add_argument(
        "--sleep",
        type=float,
        default=0.2,
        help="Seconds to sleep between Sesame requests.",
    )
    return parser.parse_args()


def load_json(path):
    return json.loads(path.read_text(encoding="utf-8"))


def write_json(path, data):
    path.write_text(json.dumps(data, indent=2, ensure_ascii=True) + "\n")


def build_source_index(sources):
    index = {}
    for source in sources:
        for value in (
            source.get("nick"),
            source.get("name"),
            source.get("latex_name"),
        ):
            if value:
                index.setdefault(normalize_name(value), source)
    return index


def collect_source_usages(issue_report):
    usages = defaultdict(list)
    for issue in issue_report["issues"]:
        if issue["kind"] != "extra_context_detection_metadata_missing":
            continue

        details = issue["details"]
        source_field = details.get("source_field")
        if not source_field or details.get("sources"):
            continue

        for source in split_sources(source_field):
            usages[source].append(
                {
                    "legacy_var": issue["legacy_var"],
                    "line": issue["line"],
                    "formula": issue["formula"],
                    "name": issue.get("name"),
                    "molecule_label": details.get("molecule_label"),
                    "context": details.get("context"),
                    "raw_source_field": source_field,
                    "observation_refs": details.get("observation_refs", []),
                }
            )
    return usages


def split_sources(source_field):
    if isinstance(source_field, list):
        values = source_field
    else:
        values = str(source_field).split(",")
    return [value.strip() for value in values if value and value.strip()]


def resolve_source_record(
    legacy_source,
    usages,
    source_index,
    used_nicks,
    skip_network=False,
    timeout=20.0,
    sleep=0.2,
):
    contexts = sorted({usage["context"] for usage in usages if usage.get("context")})
    proposed_type = proposed_source_type(contexts)

    exact_existing = source_index.get(normalize_name(legacy_source))
    if exact_existing is not None:
        return base_record(
            legacy_source,
            usages,
            contexts,
            action="map_existing",
            status="already_available",
            confidence="high",
            matched_existing=exact_existing,
            proposed_type=proposed_type,
            reason="Legacy source already matches an existing source entry.",
        )

    derived = derived_source_record(
        legacy_source, usages, contexts, source_index, used_nicks, proposed_type
    )
    if derived is not None:
        return derived

    existing_query = existing_match_query(legacy_source)
    if existing_query is not None:
        loose_existing = source_index.get(normalize_name(existing_query))
        if loose_existing is not None:
            return base_record(
                legacy_source,
                usages,
                contexts,
                action="map_existing",
                status="needs_approval",
                confidence="medium",
                matched_existing=loose_existing,
                proposed_type=proposed_type,
                reason=(
                    "Legacy source matches an existing source after removing "
                    "a descriptive prefix/suffix."
                ),
            )

    resolution = None
    if not skip_network:
        for query in query_variants(legacy_source):
            resolution = query_sesame(query, contexts, timeout=timeout)
            if resolution is not None:
                resolution["query"] = query
                break
            time.sleep(sleep)

    if resolution is None:
        return base_record(
            legacy_source,
            usages,
            contexts,
            action="manual",
            status="manual_attention",
            confidence="manual",
            matched_existing=None,
            proposed_type=proposed_type,
            reason="No resolver match was found.",
        )

    confidence, reason = confidence_for_resolution(legacy_source, resolution, contexts)
    proposed_type = proposed_source_type(contexts, resolution)
    proposed_type = SOURCE_TYPE_OVERRIDES.get(legacy_source, proposed_type)
    proposed_name = display_source_name(legacy_source, resolution)
    proposed_nick = unique_nick(nick_from_name(proposed_name), used_nicks)
    status = "needs_approval"
    approval = SOURCE_APPROVALS.get(legacy_source)
    if approval is not None:
        status = "approved"
        reason = approval["reason"]

    return base_record(
        legacy_source,
        usages,
        contexts,
        action="add_source",
        status=status,
        confidence=confidence,
        matched_existing=None,
        proposed_type=proposed_type,
        proposed_name=proposed_name,
        proposed_nick=proposed_nick,
        resolver=resolution,
        reason=reason,
    )


def derived_source_record(
    legacy_source, usages, contexts, source_index, used_nicks, proposed_type
):
    derivation = SOURCE_DERIVATIONS.get(legacy_source)
    if derivation is None:
        return None

    proposed_nick = unique_nick(derivation["proposed_nick"], used_nicks)
    if "source_nick" in derivation:
        existing = source_index.get(normalize_name(derivation["source_nick"]))
        if existing is None:
            return None
        resolver = {
            "database": "sources.json",
            "query": derivation["source_nick"],
            "matched_name": existing.get("name"),
            "otype": None,
            "ra": existing.get("ra"),
            "dec": existing.get("dec"),
            "ra_deg": None,
            "dec_deg": None,
            "jpos": (
                f"{existing.get('ra')} {existing.get('dec')}"
                if existing.get("ra") and existing.get("dec")
                else None
            ),
            "url": None,
            "simbad_url": existing.get("simbad_url"),
            "derived_from_existing_nick": existing.get("nick"),
        }
    else:
        resolver = query_sesame(derivation["query"], contexts)
        if resolver is None:
            return None
        resolver["query"] = derivation["query"]

    return base_record(
        legacy_source,
        usages,
        contexts,
        action="add_source",
        status="approved",
        confidence="high",
        matched_existing=None,
        proposed_type=proposed_type,
        proposed_name=derivation["proposed_name"],
        proposed_nick=proposed_nick,
        resolver=resolver,
        reason=derivation["reason"],
    )


def base_record(
    legacy_source,
    usages,
    contexts,
    action,
    status,
    confidence,
    matched_existing=None,
    proposed_type=None,
    proposed_name=None,
    proposed_nick=None,
    resolver=None,
    reason=None,
):
    usage_rows = sorted(usages, key=lambda row: (row["line"], row["legacy_var"]))
    return {
        "legacy_source": legacy_source,
        "contexts": contexts,
        "usage_count": len(usage_rows),
        "used_by": usage_rows,
        "action": action,
        "status": status,
        "confidence": confidence,
        "matched_existing_nick": (
            matched_existing.get("nick") if matched_existing else None
        ),
        "matched_existing_name": (
            matched_existing.get("name") if matched_existing else None
        ),
        "proposed_name": proposed_name,
        "proposed_nick": proposed_nick,
        "proposed_type": proposed_type,
        "resolver": resolver,
        "reason": reason,
    }


def proposed_source_type(contexts, resolution=None):
    if resolution is not None:
        resolver_type = source_type_from_resolver(resolution)
        if resolver_type is not None:
            return resolver_type

    types = {CONTEXT_SOURCE_TYPES.get(context, "Other") for context in contexts}
    if len(types) == 1:
        return next(iter(types))
    return "Other"


def source_type_from_resolver(resolution):
    otype = resolution.get("otype")
    if otype == "Pl":
        return "Exoplanet"
    if otype == "PN":
        return "Planetary Nebula"
    if otype in {"Or*", "Y*O", "TT*"}:
        return "YSO"
    if otype in {"G", "GiG", "GPair", "Sy1", "Sy2", "H2G", "QSO"}:
        return "External Galaxy"
    return None


def existing_match_query(legacy_source):
    text = legacy_source.strip()
    if text.endswith(" LOS"):
        return text.removesuffix(" LOS").strip()
    if text.startswith("QSO "):
        return text.removeprefix("QSO ").strip()
    return None


def query_variants(legacy_source):
    seen = set()
    variants = [legacy_source]
    variants.extend(QUERY_OVERRIDES.get(legacy_source, []))

    loose = existing_match_query(legacy_source)
    if loose:
        variants.append(loose)

    for variant in variants:
        if variant and variant not in seen:
            seen.add(variant)
            yield variant


def query_sesame(query, contexts, timeout=20.0):
    databases = "NSV" if contexts == ["exgal"] else "SNV"
    url = (
        f"{SESAME_ENDPOINT}/{databases}?"
        f"{urllib.parse.quote_plus(query)}"
    )
    request = urllib.request.Request(
        url,
        headers={"User-Agent": "astromol-source-preview/0.1"},
    )
    try:
        with urllib.request.urlopen(request, timeout=timeout) as response:
            payload = response.read()
    except (urllib.error.URLError, TimeoutError):
        return None

    try:
        root = ET.fromstring(payload)
    except ET.ParseError:
        return None

    resolver = root.find(".//Resolver")
    if resolver is None:
        return None

    matched_name = text_or_none(resolver.findtext("oname")) or query
    jpos = text_or_none(resolver.findtext("jpos"))
    ra = dec = None
    if jpos:
        parts = jpos.split()
        if len(parts) >= 2:
            ra, dec = parts[0], parts[1]

    return {
        "database": resolver.attrib.get("name"),
        "query": query,
        "matched_name": matched_name,
        "otype": text_or_none(resolver.findtext("otype")),
        "ra": ra,
        "dec": dec,
        "ra_deg": text_or_none(resolver.findtext("jradeg")),
        "dec_deg": text_or_none(resolver.findtext("jdedeg")),
        "jpos": jpos,
        "url": url,
    }


def confidence_for_resolution(legacy_source, resolution, contexts):
    query = resolution.get("query") or legacy_source
    if legacy_source.endswith(" LOS"):
        return (
            "medium",
            "Resolved after preserving the legacy LOS name but querying without or around LOS.",
        )
    if legacy_source.startswith("QSO "):
        return (
            "medium",
            "Resolved after removing the legacy QSO prefix for name lookup.",
        )
    if legacy_source == "HD 209458b" and query == "HD 209458":
        return (
            "medium",
            "Resolver fell back to the host star coordinates for the exoplanet.",
        )
    if contexts == ["ppd"]:
        return (
            "medium",
            "PPD detections currently map staged sources to source type YSO.",
        )
    if resolution.get("otype") == "PN" and contexts == ["exgal"]:
        return (
            "medium",
            "Resolved object is a substructure in an extragalactic detection context.",
        )
    return "high", "Resolved directly with a single Sesame result."


def display_source_name(legacy_source, resolution):
    if legacy_source.endswith(" LOS"):
        return legacy_source
    if legacy_source == "HD 209458b":
        return "HD 209458 b"
    return legacy_source


def record_to_source_entry(record):
    resolver = record.get("resolver") or {}
    name = record["proposed_name"]
    query_name = resolver.get("matched_name") or resolver.get("query") or name
    return {
        "_break": "==========================================",
        "name": name,
        "nick": record["proposed_nick"],
        "type": record["proposed_type"],
        "ra": resolver.get("ra"),
        "dec": resolver.get("dec"),
        "simbad_url": resolver.get("simbad_url")
        or (simbad_url(query_name) if query_name else None),
        "latex_name": name,
        "note": source_note(record),
        "_legacy_source_string": record["legacy_source"],
        "_matched_name": resolver.get("matched_name"),
        "_resolver": resolver.get("database"),
        "_resolver_query": resolver.get("query"),
        "_resolver_otype": resolver.get("otype"),
        "_confidence": record["confidence"],
        "_needs_approval": record["status"] not in {"approved", "already_available"},
        "_used_by": record["used_by"],
    }


def source_note(record):
    if record["legacy_source"].endswith(" LOS"):
        if record["legacy_source"] == "PKS 1830-211 LOS":
            return (
                "Line-of-sight source staged separately from PKS1830; coordinates "
                "and link metadata copied from PKS1830 while preserving the LOS qualifier."
            )
        if record["status"] == "approved":
            return (
                "Line-of-sight source staged separately by user approval while "
                "preserving the LOS qualifier."
            )
        return (
            "Staged from a legacy LOS source string; review whether the source "
            "should remain line-of-sight-specific."
        )
    if record["legacy_source"] == "HD 209458b":
        return "Staged from exoplanet legacy source string."
    if record["contexts"] == ["ppd"]:
        approval = SOURCE_APPROVALS.get(record["legacy_source"])
        if approval and approval.get("note"):
            return approval["note"]
        return "Staged from a PPD legacy source string; review source type."
    return None


def simbad_url(query):
    return (
        "https://simbad.u-strasbg.fr/simbad/sim-id?Ident="
        f"{urllib.parse.quote_plus(query)}"
    )


def build_report(records, usages):
    counts = Counter(record["confidence"] for record in records)
    action_counts = Counter(record["action"] for record in records)
    status_counts = Counter(record["status"] for record in records)
    detections = {
        (
            usage["legacy_var"],
            usage["line"],
            usage["molecule_label"],
            usage["context"],
        )
        for values in usages.values()
        for usage in values
    }
    return {
        "summary": {
            "detection_count": len(detections),
            "source_assignment_count": sum(len(value) for value in usages.values()),
            "unique_source_count": len(usages),
            "new_source_candidates": action_counts.get("add_source", 0),
            "existing_mappings": action_counts.get("map_existing", 0),
            "manual_attention": status_counts.get("manual_attention", 0),
            "approved": status_counts.get("approved", 0),
            "by_confidence": dict(sorted(counts.items())),
            "by_action": dict(sorted(action_counts.items())),
            "by_status": dict(sorted(status_counts.items())),
        },
        "records": records,
    }


def write_markdown_report(report):
    lines = [
        "# Additional Source Resolution Report",
        "",
        "Generated from `scripts/resolve_legacy_sources.py`.",
        "",
        "This report is a staging artifact. It does not modify `sources.json` "
        "or `detections.preview.json`.",
        "",
        "## Summary",
        "",
    ]
    summary = report["summary"]
    lines.extend(
        [
            f"- Detections with legacy source text reviewed: {summary['detection_count']}",
            f"- Individual source assignments reviewed: {summary['source_assignment_count']}",
            f"- Unique legacy source names: {summary['unique_source_count']}",
            f"- New source candidates staged: {summary['new_source_candidates']}",
            f"- Existing source mapping candidates: {summary['existing_mappings']}",
            f"- Manual attention needed: {summary['manual_attention']}",
            f"- Confidence breakdown: {json.dumps(summary['by_confidence'], sort_keys=True)}",
            "",
        ]
    )

    for confidence in ("high", "medium", "manual"):
        records = [
            record
            for record in report["records"]
            if record["confidence"] == confidence
        ]
        lines.extend(["", f"## {confidence.title()} Confidence", ""])
        if not records:
            lines.append("- None")
            continue

        for record in records:
            lines.extend(markdown_record(record))

    return "\n".join(lines) + "\n"


def count_phrase(count, singular, plural=None):
    return f"{count} {singular if count == 1 else (plural or singular + 's')}"


def markdown_record(record):
    contexts = ", ".join(record["contexts"])
    lines = [
        (
            f"- `{record['legacy_source']}` [{contexts}] "
            f"({count_phrase(record['usage_count'], 'detection')})"
        )
    ]

    if record["action"] == "map_existing":
        lines.append(
            "  - Proposed action: map to existing source "
            f"`{record['matched_existing_nick']}` "
            f"({record['matched_existing_name']})"
        )
    elif record["action"] == "add_source":
        resolver = record.get("resolver") or {}
        lines.append(
            "  - Proposed action: add staged source "
            f"`{record['proposed_nick']}` "
            f"({record['proposed_name']}, {record['proposed_type']})"
        )
        lines.append(
            "  - Resolver: "
            f"{resolver.get('database')} query `{resolver.get('query')}` "
            f"matched `{resolver.get('matched_name')}` "
            f"otype `{resolver.get('otype')}` "
            f"at `{resolver.get('jpos')}`"
        )
    else:
        lines.append("  - Proposed action: manual source curation required")

    lines.append(f"  - Reason: {record['reason']}")
    lines.append("  - Used by:")
    for usage in record["used_by"]:
        lines.append(
            "    - "
            f"line {usage['line']} `{usage['legacy_var']}` "
            f"formula `{usage['formula']}` "
            f"label `{usage['molecule_label']}`"
        )
    return lines


def normalize_name(value):
    return re.sub(r"[^a-z0-9]+", "", str(value).lower())


def nick_from_name(name):
    text = name.replace("+", "p")
    text = re.sub(r"[^A-Za-z0-9]+", "", text)
    if not text:
        return "Source"
    if text[0].isdigit():
        return f"Source{text}"
    return text


def unique_nick(nick, used_nicks):
    if nick not in used_nicks:
        return nick
    index = 2
    while f"{nick}{index}" in used_nicks:
        index += 1
    return f"{nick}{index}"


def text_or_none(value):
    if value is None:
        return None
    text = str(value).strip()
    return text or None


if __name__ == "__main__":
    main()
