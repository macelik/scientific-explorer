"""Build the offline research guide from saved explanations and research records.

This exports teaching material and existing descriptive summaries. It neither
reads individual-cell matrices nor fits or changes a scientific method.
"""

import csv
import hashlib
import json
from pathlib import Path
from urllib.parse import urlparse


HERE = Path(__file__).resolve().parent
RESEARCH = HERE.parent
FOLDERS = {
    "A": "a-cnv-methods",
    "B": "b-joint-segmentation",
    "C": "c-reference-strategies",
    "D": "d-representation-models",
    "E": "e-per-cell-evidence",
    "F": "f-technical-priors",
}
STAGES = {"inputs", "reference", "segment", "alleles", "interpret", "validate"}


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def kind(entry):
    code, category = entry["id"], entry["category"]
    if code == "A1":
        return "Current baseline"
    if code in {"A2", "A3", "A4", "B1"}:
        return "Published comparator"
    if category == "A" or code in {"B4", "B5"}:
        return "Cross-assay transfer" if code != "B5" else "Denoising concept"
    return {
        "B": "Proposed diagnostic",
        "C": "Reference strategy",
        "D": "Representation model",
        "E": "Statistical model",
        "F": "Technical context",
    }[category]


def build():
    inputs = [HERE / name for name in ("shell.html", "guide.css", "guide.js")]
    records = {}
    for category, folder in FOLDERS.items():
        path = RESEARCH / folder / "methods.json"
        inputs.append(path)
        for i, record in enumerate(json.loads(path.read_text())["methods"], 1):
            records[f"{category}{i}"] = record

    methods = []
    seen = set()
    for path in sorted((HERE / "content").glob("*.json")):
        inputs.append(path)
        for entry in json.loads(path.read_text()):
            code = entry["id"]
            assert code in records and code not in seen, code
            seen.add(code)
            assert entry["category"] == code[0]
            assert 3 <= len(entry["steps"]) <= 7, code
            assert all(s["title"] and s["text"] for s in entry["steps"])
            assert all(entry["example"][k] for k in ("title", "setup", "steps", "interpretation"))
            assert entry["adaptation"]["stage"] in STAGES
            assert all(entry["adaptation"][k] for k in ("available", "needed", "first_experiment", "success"))
            assert entry["cautions"] and entry["sources"]
            assert entry["demo"] in {"depth", "reference", "alleles", "zeros", "none"}
            assert all(urlparse(s["url"]).scheme in {"https", "http"} and s["title"] for s in entry["sources"])
            entry.update(record=records[code], originalName=records[code]["name"], kind=kind(entry))
            methods.append(entry)
    assert seen == set(records) and len(methods) == 36
    methods.sort(key=lambda m: (m["id"][0], int(m["id"][1:])))

    table = RESEARCH / "e-per-cell-evidence/artifacts/zero_support_by_stratum.tsv"
    inputs.append(table)
    with table.open() as handle:
        pilot = list(csv.DictReader(handle, delimiter="\t"))
    for row in pilot:
        for key, value in row.items():
            if key not in {"region", "stratum", "level"}:
                row[key] = float(value) if "." in value else int(value)
        assert 0 <= row["raw_whole_zero_cells"] <= row["n_cells"]
        assert 0 <= row["at_whole_zero_cells"] <= row["n_cells"]
    assert len(pilot) == 54
    assert len({(r["region"], r["stratum"], r["level"]) for r in pilot}) == 54

    data = {"methods": methods, "pilot": pilot, "pilotHash": sha(table), "edition": "2026-11-02"}
    # JSON sits in a script element; escape markup before embedding.
    encoded = json.dumps(data, ensure_ascii=False, separators=(",", ":")).replace("<", "\\u003c")
    html = (HERE / "shell.html").read_text()
    for marker, content in {
        "/*__STYLE__*/": (HERE / "guide.css").read_text(),
        "/*__DATA__*/": encoded,
        "/*__SCRIPT__*/": (HERE / "guide.js").read_text(),
    }.items():
        assert html.count(marker) == 1
        html = html.replace(marker, content)
    output = HERE / "index.html"
    output.write_text(html)
    manifest = {
        "edition": "2026-11-02",
        "entry_count": len(methods),
        "category_counts": {c: sum(m["category"] == c for m in methods) for c in FOLDERS},
        "pilot_rows": len(pilot),
        "inputs_sha256": {str(p.relative_to(RESEARCH)): sha(p) for p in inputs},
        "output_sha256": sha(output),
        "scope": "Offline educational export; no caller execution or new data fitting.",
    }
    (HERE / "build_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(f"Built {output.name}: {len(methods)} entries, {len(pilot)} saved pilot summaries, {output.stat().st_size:,} bytes")


if __name__ == "__main__":
    build()
