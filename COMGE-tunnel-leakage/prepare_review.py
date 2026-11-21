"""Reassemble the complete review data in a new directory using only Python."""
import argparse
import hashlib
import json
import tempfile
from pathlib import Path, PurePosixPath

ROOT = Path(__file__).resolve().parent
MANIFEST_SHA256 = "103f2eea9dbcadbd60e7a7a3e00b2806fe7512b35cb381f8e55887e70928a507"


def sha256(data):
    return hashlib.sha256(data).hexdigest()


def package_path(directory, name):
    """Resolve a relative, forward-slash path from the package index."""
    relative = PurePosixPath(name)
    if relative.is_absolute() or ".." in relative.parts or ":" in name or "\\" in name:
        raise ValueError(f"Invalid package path: {name}")
    path = (directory / relative).resolve()
    if not path.is_relative_to(directory.resolve()):
        raise ValueError(f"Path outside the package: {name}")
    return path


def join_parts(record, destination):
    """Join binary blocks in order, or concatenate CSV rows with one header."""
    kind = record["format"]
    if kind not in ("csv", "binary"):
        raise ValueError(f"Unknown split-file format: {kind}")
    header = None
    with destination.open("xb") as output:
        for part in record["parts"]:
            data = package_path(ROOT / "delivery_parts", part["path"]).read_bytes()
            if sha256(data) != part["sha256"]:
                raise ValueError(f"Checksum mismatch for {part['path']}")
            if kind == "csv":
                first, newline, rows = data.partition(b"\n")
                if not newline:
                    raise ValueError(f"Missing CSV header in {part['path']}")
                if header is None:
                    header = first
                elif first != header:
                    raise ValueError(f"Inconsistent CSV header in {part['path']}")
                else:
                    data = rows
            output.write(data)


def prepare(destination):
    destination = destination.resolve()
    if destination.exists() or destination.is_relative_to(ROOT):
        raise ValueError("Choose a new output directory outside the study folder.")
    manifest = (ROOT / "manifest.json").read_bytes()
    if sha256(manifest) != MANIFEST_SHA256:
        raise ValueError("The study manifest does not match this release.")
    checksums = json.loads(manifest)
    split_files = json.loads((ROOT / "delivery_parts/index.json").read_bytes())
    if split_files.keys() - checksums.keys():
        raise ValueError("The split-file index contains an unknown study file.")

    destination.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(dir=destination.parent) as temp:
        study = Path(temp) / "study"
        study.mkdir()
        for name, expected_checksum in checksums.items():
            target = package_path(study, name)
            target.parent.mkdir(parents=True, exist_ok=True)
            if name in split_files:
                join_parts(split_files[name], target)
            else:
                source = "delivery_parts/README.study.md" if name == "README.md" else name
                target.write_bytes(package_path(ROOT, source).read_bytes())
            if sha256(target.read_bytes()) != expected_checksum:
                raise ValueError(f"Reassembled file has the wrong checksum: {name}")
        (study / "manifest.json").write_bytes(manifest)
        if destination.exists():
            raise FileExistsError(destination)
        study.rename(destination)
    print(f"Prepared {len(checksums) + 1} study files in {destination}.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, required=True, help="new study directory")
    prepare(parser.parse_args().output_dir)
