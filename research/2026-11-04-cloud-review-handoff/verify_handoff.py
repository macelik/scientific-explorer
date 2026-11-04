"""Verify an extracted cloud data handoff and the separately obtained pyEpi source."""
import argparse
import hashlib
import json
from pathlib import Path


def sha(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-root", type=Path, required=True)
    parser.add_argument("--package-dir", type=Path, required=True)
    args = parser.parse_args()
    manifest = json.loads(Path(__file__).with_name("DATA_MANIFEST.json").read_text())
    embedded = args.data_root / "MANIFEST.json"
    errors = []
    if not embedded.is_file() or json.loads(embedded.read_text()) != manifest:
        errors.append("Embedded manifest is missing or differs from this handoff.")
    for record in manifest["files"]:
        path = args.data_root / record["path"]
        if not path.is_file():
            errors.append("Missing data: " + record["path"])
        elif path.stat().st_size != record["size_bytes"] or sha(path) != record["sha256"]:
            errors.append("Data identity mismatch: " + record["path"])
    for name, expected in manifest["upstream"]["package_files_sha256"].items():
        path = args.package_dir / name
        if not path.is_file() or sha(path) != expected:
            errors.append("Upstream source mismatch: " + name)
    print(json.dumps({"status": "failed" if errors else "passed", "data_files": len(manifest["files"]),
                      "upstream_files": len(manifest["upstream"]["package_files_sha256"]),
                      "errors": errors, "scope": "Byte identity only; not scientific validation."}, indent=2))
    raise SystemExit(bool(errors))


if __name__ == "__main__":
    main()
