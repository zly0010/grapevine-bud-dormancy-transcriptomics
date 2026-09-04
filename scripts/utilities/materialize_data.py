"""Restore archived input layout by verified copying only. Never invokes analysis."""
from pathlib import Path
import argparse, csv, hashlib, shutil

def digest(path):
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(4 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--archive", required=True, type=Path)
    parser.add_argument("--destination", required=True, type=Path)
    args = parser.parse_args()
    archive, destination = args.archive.resolve(), args.destination.resolve()
    mapping = Path(__file__).resolve().parents[2] / "docs" / "DATA_LAYOUT.tsv"
    with mapping.open(encoding="utf-8", newline="") as f:
        rows = list(csv.DictReader(f, delimiter="\t"))
    # Validate the entire copy plan before creating files.
    pending = []
    for row in rows:
        source = (archive / row["package_path"]).resolve()
        target = (destination / row["workspace_path"]).resolve()
        if not source.is_relative_to(archive) or not target.is_relative_to(destination):
            raise ValueError("Path escapes selected directories")
        if digest(source) != row["sha256"]:
            raise ValueError("Archive hash mismatch: " + row["package_path"])
        if target.exists():
            if digest(target) != row["sha256"]:
                raise FileExistsError("Refusing to overwrite: " + row["workspace_path"])
        else:
            pending.append((source, target, row["sha256"]))
    for source, target, expected in pending:
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, target)
        if digest(target) != expected:
            raise ValueError("Copied file hash mismatch")
    print(f"Copied {len(pending)} files; no analysis was invoked.")

if __name__ == "__main__":
    main()
