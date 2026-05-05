"""Fetch candidate missing references from NASA ADS into a temporary BibTeX file.

The script reads unresolved reference issues from
``astromol/data/legacy_conversion_issues.json`` and tries to resolve them
against ADS. It writes temporary import artifacts only; it never modifies
``references.bib``.

Authentication:
  - ``ADS_API_TOKEN`` environment variable, or
  - ``~/.ads/token``, or
  - ``--token-stdin`` to read one line from stdin.
"""

from __future__ import annotations

import argparse
import json
import os
import re
import sys
import time
import unicodedata
import urllib.error
import urllib.parse
import urllib.request
from collections import Counter
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = REPO_ROOT / "astromol" / "data"
ISSUES_PATH = DATA_DIR / "legacy_conversion_issues.json"
REFERENCES_PATH = DATA_DIR / "references.bib"

DEFAULT_BIB_OUT = DATA_DIR / "ads_missing_refs.bib"
DEFAULT_REPORT_OUT = DATA_DIR / "ads_missing_refs_report.md"
DEFAULT_JSON_OUT = DATA_DIR / "ads_missing_refs_report.json"

ADS_BASE_URL = "https://api.adsabs.harvard.edu/v1"


class ADSMissingReferenceFetcher:
    def __init__(self, token: str, sleep_seconds: float = 0.15):
        self.token = token
        self.sleep_seconds = sleep_seconds
        self.current_ref_ids = self._load_current_ref_ids()
        self.search_cache = {}
        self.export_cache = {}

    def run(self):
        items = self._load_items()
        resolved = []
        unresolved = []
        ambiguous = []

        for item in items:
            match = self._resolve_item(item)
            if match["status"] == "resolved":
                resolved.append(match)
            elif match["status"] == "ambiguous":
                ambiguous.append(match)
            else:
                unresolved.append(match)

        unique_bibcodes = sorted({match["bibcode"] for match in resolved})
        bibtex = self._export_bibtex(unique_bibcodes)
        bibtex, normalized_publications = self._zotero_friendly_bibtex(
            bibtex, resolved
        )
        bibtex, rewritten_citekeys, citekey_by_bibcode = (
            self._citekey_friendly_bibtex(bibtex)
        )
        report = self._build_report(
            items,
            resolved,
            ambiguous,
            unresolved,
            bibtex,
            normalized_publications,
            rewritten_citekeys,
            citekey_by_bibcode,
        )

        return bibtex, report

    def _load_current_ref_ids(self):
        if not REFERENCES_PATH.exists():
            return set()
        text = REFERENCES_PATH.read_text()
        return set(re.findall(r"@\w+\{([^,\s]+)", text))

    def _load_items(self):
        with open(ISSUES_PATH) as handle:
            report = json.load(handle)

        items = []
        seen = set()
        for issue in report["issues"]:
            item = self._item_from_issue(issue)
            if item is None:
                continue

            dedupe_key = (
                item["kind"],
                item.get("source_ref"),
                item.get("text"),
                item.get("year"),
                item.get("surname"),
                item.get("volume"),
                item.get("page"),
            )
            if dedupe_key in seen:
                continue
            seen.add(dedupe_key)
            items.append(item)

        return items

    def _item_from_issue(self, issue):
        kind = issue["kind"]
        details = issue["details"]

        if kind == "unresolved_reference":
            ref = details["bibcode"]
            if ref in self.current_ref_ids:
                return None

            if self._is_ads_bibcode(ref):
                return {
                    "kind": "ads_bibcode",
                    "source_ref": ref,
                    "issue": self._issue_context(issue),
                }

            parsed = self._parse_legacy_key(ref)
            if parsed is None:
                return {
                    "kind": "unparsed_reference_key",
                    "source_ref": ref,
                    "issue": self._issue_context(issue),
                }
            return {
                "kind": "parsed_reference_key",
                "source_ref": ref,
                "surname": parsed["surname"],
                "year": parsed["year"],
                "volume": None,
                "page": parsed["page"],
                "issue": self._issue_context(issue),
            }

        if kind in {"unresolved_free_text_reference", "free_text_reference_unparsed"}:
            text = details["text"]
            parsed = self._parse_free_text_reference(text)
            item = {
                "kind": "free_text_reference",
                "text": text,
                "issue": self._issue_context(issue),
            }
            if parsed is None:
                item["kind"] = "unparsed_free_text_reference"
            else:
                item.update(parsed)
            return item

        return None

    def _issue_context(self, issue):
        return {
            "legacy_var": issue.get("legacy_var"),
            "formula": issue.get("formula"),
            "line": issue.get("line"),
        }

    def _parse_legacy_key(self, ref):
        match = re.match(r"^([^:]+):(\d{4}):(.+)$", ref)
        if not match:
            return None
        page = self._normalize_page(match.group(3))
        if not page:
            return None
        return {
            "surname": self._normalize_lookup_text(match.group(1)),
            "year": match.group(2),
            "page": page,
        }

    def _parse_free_text_reference(self, text):
        year_match = re.search(r"\b(18|19|20)\d{2}\b", str(text))
        if not year_match:
            return None

        author_text = str(text)[: year_match.start()].strip(" ,")
        if not author_text:
            return None

        author_text = re.split(
            r"\s+(?:et\s+al\.?|and|&)(?:\s+|$)", author_text, maxsplit=1
        )[0]
        surname = self._normalize_lookup_text(author_text)
        volume, page = self._volume_page_from_reference_text(
            str(text)[year_match.end() :]
        )
        if not surname or not page:
            return None

        return {
            "surname": surname,
            "year": year_match.group(0),
            "volume": volume,
            "page": page,
        }

    def _volume_page_from_reference_text(self, text):
        text = text.strip()
        range_match = re.search(
            r"(?:,\s*|pp\.?\s+)([A-Za-z]?\d+[A-Za-z]?)(?:\s*[-–]\s*[A-Za-z]?\d+[A-Za-z]?)?\s*$",
            text,
            re.IGNORECASE,
        )
        if range_match:
            tokens = re.findall(r"\b\d+[A-Za-z]?\b", text[: range_match.start()])
            volume = self._normalize_volume(tokens[-1]) if tokens else None
            return volume, self._normalize_page(range_match.group(1))

        tokens = re.findall(r"\b[A-Za-z]?\d+[A-Za-z]?\b", text)
        if not tokens:
            return None, None
        volume = self._normalize_volume(tokens[-2]) if len(tokens) > 1 else None
        return volume, self._normalize_page(tokens[-1])

    def _resolve_item(self, item):
        if item["kind"] == "ads_bibcode":
            docs = self._search_ads(f'bibcode:"{item["source_ref"]}"')
            docs = [doc for doc in docs if doc.get("bibcode") == item["source_ref"]]
            if len(docs) == 1:
                return self._resolved(item, docs[0], "direct_bibcode")
            return self._unresolved(item, "ADS bibcode was not found exactly.")

        if item["kind"] in {
            "parsed_reference_key",
            "free_text_reference",
        }:
            docs = self._search_by_metadata(item)
            matching_docs = [doc for doc in docs if self._metadata_matches(item, doc)]
            if len(matching_docs) == 1:
                return self._resolved(item, matching_docs[0], "metadata_search")
            if len(matching_docs) > 1:
                return self._ambiguous(item, matching_docs, "multiple metadata matches")
            return self._unresolved(item, "No unique metadata-consistent ADS match.")

        return self._unresolved(item, "Reference could not be parsed into search metadata.")

    def _search_by_metadata(self, item):
        surname = item["surname"]
        year = item["year"]
        volume = item.get("volume")
        page = item["page"]

        queries = []
        if volume:
            queries.append(f'author:"{surname}" year:{year} volume:"{volume}" page:"{page}"')
            queries.append(f'author:"{surname}" year:{year} "{volume}" "{page}"')
        queries.append(f'author:"{surname}" year:{year} page:"{page}"')
        queries.append(f'author:"{surname}" year:{year} "{page}"')

        docs = []
        for query in queries:
            docs.extend(self._search_ads(query))
            if docs:
                break
        return self._unique_docs(docs)

    def _search_ads(self, query):
        if query in self.search_cache:
            return self.search_cache[query]

        params = {
            "q": query,
            "fl": "bibcode,title,author,year,pub,volume,page,doi",
            "rows": "10",
        }
        url = f"{ADS_BASE_URL}/search/query?{urllib.parse.urlencode(params)}"
        data = self._request_json("GET", url)
        docs = data.get("response", {}).get("docs", [])
        self.search_cache[query] = docs
        time.sleep(self.sleep_seconds)
        return docs

    def _export_bibtex(self, bibcodes):
        if not bibcodes:
            return ""

        chunks = []
        for start in range(0, len(bibcodes), 50):
            chunk = bibcodes[start : start + 50]
            chunks.append(self._export_bibtex_chunk(chunk))
        return "\n\n".join(part.strip() for part in chunks if part.strip()) + "\n"

    def _export_bibtex_chunk(self, bibcodes):
        key = tuple(bibcodes)
        if key in self.export_cache:
            return self.export_cache[key]

        url = f"{ADS_BASE_URL}/export/bibtex"
        payload = {
            "bibcode": bibcodes,
            "sort": "no sort",
            "authorlimit": 500,
        }
        data = self._request_json("POST", url, payload)
        export = data.get("export", "")
        self.export_cache[key] = export
        time.sleep(self.sleep_seconds)
        return export

    def _zotero_friendly_bibtex(self, bibtex, resolved):
        """Replace ADS journal macros with explicit publication strings."""
        pub_by_bibcode = {
            match["bibcode"]: match["pub"]
            for match in resolved
            if match.get("bibcode") and match.get("pub")
        }

        entries = re.split(r"\n(?=@)", bibtex.strip())
        normalized_entries = []
        normalized_count = 0

        for entry in entries:
            key_match = re.match(r"@\w+\{([^,\s]+),", entry)
            if key_match is None:
                normalized_entries.append(entry)
                continue

            bibcode = self._bibcode_from_entry(entry, key_match.group(1))
            publication = pub_by_bibcode.get(bibcode)
            if not publication:
                normalized_entries.append(entry)
                continue

            lines = entry.splitlines()
            changed = False
            for index, line in enumerate(lines):
                if re.match(r"\s*journal\s*=", line, flags=re.IGNORECASE):
                    lines[index] = (
                        f"      journal = {{{self._escape_bibtex_value(publication)}}},"
                    )
                    changed = True
                    break

            if not changed:
                insert_at = 1
                for index, line in enumerate(lines):
                    if re.match(r"\s*title\s*=", line, flags=re.IGNORECASE):
                        insert_at = index + 1
                        break
                lines.insert(
                    insert_at,
                    f"      journal = {{{self._escape_bibtex_value(publication)}}},",
                )

            normalized_count += 1
            normalized_entries.append("\n".join(lines))

        return "\n\n".join(normalized_entries).strip() + "\n", normalized_count

    def _citekey_friendly_bibtex(self, bibtex):
        """Rewrite ADS BibTeX keys to LastName:Year:FirstPage keys."""
        entries = self._split_bibtex_entries(bibtex)
        normalized_entries = []
        normalized_count = 0
        used_citekeys = set(self.current_ref_ids)
        citekey_by_bibcode = {}

        for entry in entries:
            key_match = re.match(r"@\w+\{([^,\s]+),", entry)
            if key_match is None:
                normalized_entries.append(entry)
                continue

            old_key = key_match.group(1)
            bibcode = self._bibcode_from_entry(entry, old_key)
            citekey_base = self._citekey_base_from_entry(entry)
            if not citekey_base:
                normalized_entries.append(entry)
                continue

            citekey = self._deduplicated_citekey(citekey_base, used_citekeys)
            used_citekeys.add(citekey)
            if bibcode:
                citekey_by_bibcode[bibcode] = citekey

            if old_key != citekey:
                entry = re.sub(
                    r"^(@\w+\{)[^,\s]+",
                    lambda match: f"{match.group(1)}{citekey}",
                    entry,
                    count=1,
                )
                normalized_count += 1
            normalized_entries.append(entry)

        return (
            "\n\n".join(normalized_entries).strip() + "\n",
            normalized_count,
            citekey_by_bibcode,
        )

    def _split_bibtex_entries(self, bibtex):
        text = bibtex.strip()
        if not text:
            return []
        return re.split(r"\n\s*(?=@)", text)

    def _bibcode_from_entry(self, entry, key):
        if self._is_ads_bibcode(key):
            return key

        adsurl = self._bibtex_field_value(entry, "adsurl")
        if adsurl:
            match = re.search(r"/abs/([^/\s}]+)", adsurl)
            if match:
                return match.group(1)
        return None

    def _citekey_base_from_entry(self, entry):
        author = self._bibtex_field_value(entry, "author")
        year = self._bibtex_field_value(entry, "year")
        pages = self._bibtex_field_value(entry, "pages")
        if not pages:
            pages = self._bibtex_field_value(entry, "eid")

        surname = self._first_author_surname_for_citekey(author)
        first_page = self._first_page_for_citekey(pages)
        year = re.sub(r"[^0-9]", "", str(year or ""))
        if not surname or not year or not first_page:
            return None
        return f"{surname}:{year}:{first_page}"

    def _first_author_surname_for_citekey(self, author_value):
        if not author_value:
            return None

        authors = self._split_bibtex_authors(author_value)
        if not authors:
            return None

        first_author = authors[0].strip()
        double_wrapped = first_author.startswith("{{") and first_author.endswith("}}")
        name = self._strip_wrapping_braces(first_author)
        comma_index = self._top_level_delimiter_index(name, ",")
        if comma_index >= 0:
            surname = name[:comma_index]
        elif double_wrapped:
            surname = name
        else:
            parts = self._latex_to_text(name).split()
            surname = parts[-1] if parts else name

        surname = self._latex_to_text(surname)
        surname = re.sub(r"\s+", "", surname)
        surname = re.sub(r"[^A-Za-z0-9-]", "", surname)
        surname = self._normalize_citekey_surname_case(surname)
        return surname or None

    def _normalize_citekey_surname_case(self, surname):
        if re.match(r"^Mc[a-z]", surname):
            return f"Mc{surname[2].upper()}{surname[3:]}"
        return surname

    def _split_bibtex_authors(self, author_value):
        authors = []
        start = 0
        depth = 0
        index = 0
        while index < len(author_value):
            char = author_value[index]
            if char == "{":
                depth += 1
            elif char == "}":
                depth = max(0, depth - 1)
            elif depth == 0 and author_value.startswith(" and ", index):
                authors.append(author_value[start:index])
                index += 5
                start = index
                continue
            index += 1

        authors.append(author_value[start:])
        return [author for author in authors if author.strip()]

    def _first_page_for_citekey(self, pages):
        if not pages:
            return None
        first_page = re.split(r"\s*[-–—]\s*", str(pages).strip(), maxsplit=1)[0]
        first_page = self._latex_to_text(first_page)
        first_page = re.sub(r"[^A-Za-z0-9]", "", first_page)
        return first_page or None

    def _deduplicated_citekey(self, base, used_citekeys):
        if base not in used_citekeys:
            return base

        suffix_index = 0
        while True:
            suffix = chr(ord("a") + suffix_index)
            candidate = f"{base}{suffix}"
            if candidate not in used_citekeys:
                return candidate
            suffix_index += 1

    def _bibtex_field_value(self, entry, field):
        match = re.search(rf"(?im)^\s*{re.escape(field)}\s*=\s*", entry)
        if match is None:
            return None

        index = match.end()
        while index < len(entry) and entry[index].isspace():
            index += 1
        if index >= len(entry):
            return None

        if entry[index] == "{":
            start = index + 1
            depth = 1
            index += 1
            while index < len(entry):
                if entry[index] == "{":
                    depth += 1
                elif entry[index] == "}":
                    depth -= 1
                    if depth == 0:
                        return entry[start:index].strip()
                index += 1
            return None

        if entry[index] == '"':
            start = index + 1
            index += 1
            while index < len(entry):
                if entry[index] == '"' and entry[index - 1] != "\\":
                    return entry[start:index].strip()
                index += 1
            return None

        start = index
        while index < len(entry) and entry[index] not in ",\n":
            index += 1
        return entry[start:index].strip()

    def _strip_wrapping_braces(self, text):
        text = text.strip()
        while text.startswith("{") and text.endswith("}"):
            depth = 0
            wraps_whole_text = False
            for index, char in enumerate(text):
                if char == "{":
                    depth += 1
                elif char == "}":
                    depth -= 1
                    if depth == 0:
                        wraps_whole_text = index == len(text) - 1
                        break
            if not wraps_whole_text:
                break
            text = text[1:-1].strip()
        return text

    def _top_level_delimiter_index(self, text, delimiter):
        depth = 0
        for index, char in enumerate(text):
            if char == "{":
                depth += 1
            elif char == "}":
                depth = max(0, depth - 1)
            elif char == delimiter and depth == 0:
                return index
        return -1

    def _latex_to_text(self, value):
        text = str(value or "")
        text = re.sub(r"\\['\"`^~=.]\\?([A-Za-z])", r"\1", text)
        text = re.sub(r"\\['\"`^~=.]\{\\?([A-Za-z])\}", r"\1", text)
        text = text.replace(r"\i", "i")
        text = text.replace(r"\o", "o")
        text = text.replace(r"\O", "O")
        text = text.replace(r"\ae", "ae")
        text = text.replace(r"\AE", "AE")
        text = text.replace(r"\aa", "a")
        text = text.replace(r"\AA", "A")
        text = re.sub(r"\\[A-Za-z]+", "", text)
        text = text.replace("{", "").replace("}", "")
        text = unicodedata.normalize("NFKD", text)
        return "".join(ch for ch in text if not unicodedata.combining(ch))

    def _escape_bibtex_value(self, value):
        return str(value).replace("\\", "\\textbackslash{}").replace("{", "").replace("}", "")

    def _request_json(self, method, url, payload=None):
        body = None
        headers = {
            "Authorization": f"Bearer {self.token}",
            "Content-Type": "application/json",
        }
        if payload is not None:
            body = json.dumps(payload).encode("utf-8")

        request = urllib.request.Request(url, data=body, headers=headers, method=method)
        try:
            with urllib.request.urlopen(request, timeout=30) as response:
                return json.loads(response.read().decode("utf-8"))
        except urllib.error.HTTPError as exc:
            message = exc.read().decode("utf-8", errors="replace")
            raise RuntimeError(f"ADS API request failed ({exc.code}): {message}") from exc

    def _metadata_matches(self, item, doc):
        if str(doc.get("year")) != str(item["year"]):
            return False

        first_author = self._first_author(doc.get("author") or [])
        if first_author and self._normalize_lookup_text(first_author) != item["surname"]:
            return False

        volume = item.get("volume")
        if volume and self._normalize_volume(doc.get("volume")) != volume:
            return False

        doc_pages = self._doc_pages(doc)
        if item["page"] not in doc_pages:
            return False

        return True

    def _doc_pages(self, doc):
        pages = doc.get("page") or []
        if isinstance(pages, str):
            pages = [pages]
        return {self._normalize_page(page) for page in pages if self._normalize_page(page)}

    def _first_author(self, authors):
        if not authors:
            return None
        author = authors[0]
        if "," in author:
            return author.split(",", 1)[0]
        return author.split()[-1] if author.split() else None

    def _resolved(self, item, doc, strategy):
        return {
            "status": "resolved",
            "strategy": strategy,
            "item": item,
            "bibcode": doc["bibcode"],
            "author": doc.get("author"),
            "title": self._first_value(doc.get("title")),
            "year": doc.get("year"),
            "pub": doc.get("pub"),
            "volume": doc.get("volume"),
            "page": doc.get("page"),
        }

    def _ambiguous(self, item, docs, reason):
        return {
            "status": "ambiguous",
            "reason": reason,
            "item": item,
            "candidates": [
                {
                    "bibcode": doc.get("bibcode"),
                    "title": self._first_value(doc.get("title")),
                    "year": doc.get("year"),
                    "pub": doc.get("pub"),
                    "volume": doc.get("volume"),
                    "page": doc.get("page"),
                }
                for doc in docs
            ],
        }

    def _unresolved(self, item, reason):
        return {
            "status": "unresolved",
            "reason": reason,
            "item": item,
        }

    def _build_report(
        self,
        items,
        resolved,
        ambiguous,
        unresolved,
        bibtex,
        normalized_publications,
        rewritten_citekeys,
        citekey_by_bibcode,
    ):
        by_kind = Counter(item["kind"] for item in items)
        resolved = self._resolved_with_citekeys(resolved, citekey_by_bibcode)
        return {
            "summary": {
                "items_considered": len(items),
                "items_by_kind": dict(sorted(by_kind.items())),
                "resolved": len(resolved),
                "ambiguous": len(ambiguous),
                "unresolved": len(unresolved),
                "unique_bibcodes": len({match["bibcode"] for match in resolved}),
                "bibtex_entries": len(re.findall(r"@\w+\{", bibtex)),
                "publication_fields_normalized": normalized_publications,
                "citekeys_normalized": len(citekey_by_bibcode),
                "citekeys_rewritten": rewritten_citekeys,
            },
            "resolved": resolved,
            "ambiguous": ambiguous,
            "unresolved": unresolved,
        }

    def _resolved_with_citekeys(self, resolved, citekey_by_bibcode):
        updated = []
        for match in resolved:
            match = dict(match)
            citekey = citekey_by_bibcode.get(match.get("bibcode"))
            if citekey:
                match["citekey"] = citekey
            updated.append(match)
        return updated

    def _unique_docs(self, docs):
        seen = set()
        unique = []
        for doc in docs:
            bibcode = doc.get("bibcode")
            if bibcode and bibcode not in seen:
                seen.add(bibcode)
                unique.append(doc)
        return unique

    def _normalize_page(self, value):
        text = re.sub(r"[^A-Za-z0-9]", "", str(value or "")).upper()
        return text or None

    def _normalize_volume(self, value):
        text = re.sub(r"[^A-Za-z0-9]", "", str(value or "")).upper()
        return text or None

    def _normalize_lookup_text(self, value):
        text = str(value or "")
        text = re.sub(r"\\['\"`^~=.]\\{?([A-Za-z])\\}?", r"\1", text)
        text = text.replace("{", "").replace("}", "")
        text = unicodedata.normalize("NFKD", text)
        text = "".join(ch for ch in text if not unicodedata.combining(ch))
        return re.sub(r"[^A-Za-z0-9]", "", text).lower()

    def _is_ads_bibcode(self, value):
        return bool(re.match(r"^\d{4}.{14}[A-Z]$", str(value)))

    def _first_value(self, value):
        if isinstance(value, list):
            return value[0] if value else None
        return value


def write_markdown_report(path, report, bib_out):
    summary = report["summary"]
    lines = [
        "# ADS Missing References Report",
        "",
        "This is a temporary import report. The active `references.bib` file was not modified.",
        "",
        f"Temporary BibTeX file: `{bib_out}`",
        "",
        "## Summary",
        "",
        f"- Items considered: {summary['items_considered']}",
        f"- Resolved: {summary['resolved']}",
        f"- Ambiguous: {summary['ambiguous']}",
        f"- Unresolved: {summary['unresolved']}",
        f"- Unique ADS bibcodes exported: {summary['unique_bibcodes']}",
        f"- BibTeX entries written: {summary['bibtex_entries']}",
        f"- Publication fields normalized: {summary.get('publication_fields_normalized', 0)}",
        f"- Citekeys normalized: {summary.get('citekeys_normalized', 0)}",
        f"- Citekeys rewritten this run: {summary.get('citekeys_rewritten', 0)}",
        "",
        "## Resolved",
        "",
    ]

    for match in report["resolved"]:
        item = match["item"]
        source = item.get("source_ref") or item.get("text")
        citekey = match.get("citekey") or match["bibcode"]
        lines.append(
            f"- `{source}` -> `{citekey}` (ADS `{match['bibcode']}`, {match.get('year')}, {match.get('pub')}, {match.get('volume')}, {match.get('page')})"
        )

    lines.extend(["", "## Ambiguous", ""])
    if report["ambiguous"]:
        for match in report["ambiguous"]:
            item = match["item"]
            source = item.get("source_ref") or item.get("text")
            lines.append(f"- `{source}`: {match['reason']}")
            for candidate in match["candidates"]:
                lines.append(
                    f"  - `{candidate['bibcode']}` ({candidate.get('year')}, {candidate.get('pub')}, {candidate.get('volume')}, {candidate.get('page')})"
                )
    else:
        lines.append("- None")

    lines.extend(["", "## Unresolved", ""])
    if report["unresolved"]:
        for match in report["unresolved"]:
            item = match["item"]
            source = item.get("source_ref") or item.get("text")
            lines.append(f"- `{source}`: {match['reason']}")
    else:
        lines.append("- None")

    path.write_text("\n".join(lines) + "\n")


def read_token(args):
    if args.token_stdin:
        token = sys.stdin.readline().strip()
    else:
        token = os.environ.get("ADS_API_TOKEN")
        if not token:
            token_path = Path.home() / ".ads" / "token"
            if token_path.exists():
                token = token_path.read_text().strip()

    if not token:
        raise SystemExit(
            "Missing ADS API token. Use ADS_API_TOKEN, ~/.ads/token, or --token-stdin."
        )
    return token


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--token-stdin", action="store_true")
    parser.add_argument(
        "--normalize-existing",
        action="store_true",
        help="Rewrite an existing BibTeX output using cached ADS metadata in the JSON report.",
    )
    parser.add_argument("--bib-out", type=Path, default=DEFAULT_BIB_OUT)
    parser.add_argument("--report-out", type=Path, default=DEFAULT_REPORT_OUT)
    parser.add_argument("--json-out", type=Path, default=DEFAULT_JSON_OUT)
    parser.add_argument("--sleep", type=float, default=0.15)
    return parser.parse_args()


def main():
    args = parse_args()

    if args.normalize_existing:
        bibtex = args.bib_out.read_text()
        report = json.loads(args.json_out.read_text())
        fetcher = ADSMissingReferenceFetcher(token="", sleep_seconds=args.sleep)
        bibtex, normalized_publications = fetcher._zotero_friendly_bibtex(
            bibtex, report.get("resolved", [])
        )
        bibtex, rewritten_citekeys, citekey_by_bibcode = (
            fetcher._citekey_friendly_bibtex(bibtex)
        )
        report["resolved"] = fetcher._resolved_with_citekeys(
            report.get("resolved", []), citekey_by_bibcode
        )
        report.setdefault("summary", {})[
            "publication_fields_normalized"
        ] = normalized_publications
        report["summary"]["citekeys_normalized"] = len(citekey_by_bibcode)
        report["summary"]["citekeys_rewritten"] = rewritten_citekeys
        report["summary"]["bibtex_entries"] = len(re.findall(r"@\w+\{", bibtex))

        args.bib_out.write_text(bibtex)
        args.json_out.write_text(json.dumps(report, indent=2, ensure_ascii=True) + "\n")
        write_markdown_report(args.report_out, report, args.bib_out)

        print(
            f"Normalized {normalized_publications} publication fields in {args.bib_out}"
        )
        print(
            f"Normalized {len(citekey_by_bibcode)} BibTeX citekeys in {args.bib_out}"
        )
        print(
            f"Rewrote {rewritten_citekeys} BibTeX citekeys in {args.bib_out}"
        )
        print(f"Wrote report to {args.report_out}")
        print(json.dumps(report["summary"], indent=2, sort_keys=True))
        return

    token = read_token(args)
    fetcher = ADSMissingReferenceFetcher(token=token, sleep_seconds=args.sleep)
    bibtex, report = fetcher.run()

    args.bib_out.write_text(bibtex)
    args.json_out.write_text(json.dumps(report, indent=2, ensure_ascii=True) + "\n")
    write_markdown_report(args.report_out, report, args.bib_out)

    print(f"Wrote {report['summary']['bibtex_entries']} BibTeX entries to {args.bib_out}")
    print(f"Wrote report to {args.report_out}")
    print(json.dumps(report["summary"], indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
