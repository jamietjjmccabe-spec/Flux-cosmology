#!/usr/bin/env python3
"""Fetch and preserve a Harvard Dataverse dataset manifest.

Default target: Overstreet et al. replication data, DOI 10.7910/DVN/O9IS7Y.
The tool preserves the complete API response, records checksums and metadata,
and writes JSON/CSV manifests. Filename heuristics are only triage aids; they
do not prove that shot-level contrast data are recoverable.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import time
import urllib.parse
import urllib.request
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

DEFAULT_DOI = "doi:10.7910/DVN/O9IS7Y"
DEFAULT_SERVER = "https://dataverse.harvard.edu"
RAW_TERMS = ("raw", "shot", "population", "port", "fluorescence", "image", "fringe", "ellipse", "counts")


def fetch(url: str, timeout: float, retries: int) -> bytes:
    req = urllib.request.Request(
        url,
        headers={
            "User-Agent": "Flux-MIP-V2.1-Manifest-Audit/1.0",
            "Accept": "application/json",
        },
    )
    error: Exception | None = None
    for attempt in range(retries + 1):
        try:
            with urllib.request.urlopen(req, timeout=timeout) as response:
                return response.read()
        except Exception as exc:  # network errors are reported, never hidden
            error = exc
            if attempt < retries:
                time.sleep(2**attempt)
    raise RuntimeError(f"Dataverse request failed: {error}")


def classify(entry: dict[str, Any]) -> dict[str, Any]:
    data_file = entry.get("dataFile", {}) or {}
    name = str(data_file.get("filename", "UNKNOWN"))
    description = str(entry.get("description", ""))
    directory = str(entry.get("directoryLabel", ""))
    content_type = str(data_file.get("contentType", ""))
    blob = " ".join((name, description, directory, content_type)).lower()
    hits = [term for term in RAW_TERMS if term in blob]
    checksum = data_file.get("checksum", {}) or {}
    return {
        "file_id": data_file.get("id"),
        "filename": name,
        "directory_label": directory,
        "description": description,
        "content_type": content_type,
        "filesize_bytes": int(data_file.get("filesize") or 0),
        "restricted": bool(entry.get("restricted", False)),
        "checksum_type": checksum.get("type", ""),
        "checksum_value": checksum.get("value", data_file.get("md5", "")),
        "raw_indicator_hits": hits,
        "classification": "possible raw-contrast candidate" if hits else "indeterminate/processed",
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--doi", default=DEFAULT_DOI)
    parser.add_argument("--server", default=DEFAULT_SERVER)
    parser.add_argument("--output-dir", type=Path, default=Path("overstreet_manifest"))
    parser.add_argument("--timeout", type=float, default=60.0)
    parser.add_argument("--retries", type=int, default=3)
    args = parser.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)
    encoded = urllib.parse.quote(args.doi, safe="")
    url = f"{args.server.rstrip('/')}/api/datasets/:persistentId/?persistentId={encoded}"
    raw = fetch(url, args.timeout, args.retries)
    api_path = args.output_dir / "dataverse_api_response.json"
    api_path.write_bytes(raw)
    payload = json.loads(raw.decode("utf-8"))
    latest = payload.get("data", {}).get("latestVersion", {}) or {}
    rows = [classify(item) for item in latest.get("files", []) or []]

    csv_path = args.output_dir / "file_manifest.csv"
    fieldnames = list(rows[0]) if rows else ["filename"]
    with csv_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    candidates = [row["filename"] for row in rows if row["raw_indicator_hits"]]
    summary = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "query_url": url,
        "doi": args.doi,
        "api_response_sha256": hashlib.sha256(raw).hexdigest(),
        "dataset_version": f"{latest.get('versionNumber', '')}.{latest.get('versionMinorNumber', '')}".strip("."),
        "release_time": latest.get("releaseTime"),
        "file_count": len(rows),
        "total_size_bytes": sum(row["filesize_bytes"] for row in rows),
        "files": rows,
        "step2_assessment": {
            "status": "provisionally inspect candidates" if candidates else "not established from manifest",
            "candidate_files": candidates,
            "rule": "Do not derive a contrast bound until file contents expose output-port populations, fringe images, or an equivalent per-shot contrast observable separated by source state.",
        },
    }
    summary_path = args.output_dir / "manifest_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary["step2_assessment"], indent=2))
    print(f"Wrote {api_path}, {csv_path}, and {summary_path}")


if __name__ == "__main__":
    main()
