"""Fetch only selected calibration members of remote ZIP files with HTTP ranges.

Uses standard-library zipfile and validates every Content-Range, member CRC and
member SHA-256. The archive's published whole-file checksum is recorded, not
claimed verified when only some ranges were downloaded. Inputs stay under --work.
"""

from __future__ import annotations

import argparse
import io
import json
from pathlib import Path, PurePosixPath
import urllib.request
import zipfile

from ngesh.bundle import json_text


class RangeReader(io.RawIOBase):
    def __init__(self, url: str, size: int, budget: int = 128 * 1024 * 1024) -> None:
        self.url, self.size, self.budget = url, size, budget
        self.pos = self.fetched = 0
        self.cache_start = 0
        self.cache = b""

    def seekable(self) -> bool:
        return True

    def readable(self) -> bool:
        return True

    def tell(self) -> int:
        return self.pos

    def seek(self, offset: int, whence: int = 0) -> int:
        self.pos = (
            offset if whence == 0 else self.pos + offset if whence == 1 else self.size + offset
        )
        if self.pos < 0:
            raise ValueError("negative archive seek")
        return self.pos

    def fetch(self, start: int, count: int) -> bytes:
        count = min(count, self.size - start)
        if self.fetched + count > self.budget:
            raise ValueError("selective archive download exceeds explicit byte budget")
        request = urllib.request.Request(
            self.url,
            headers={
                "Range": f"bytes={start}-{start + count - 1}",
                "User-Agent": "ngesh-calibration/2",
            },
        )
        with urllib.request.urlopen(request, timeout=90) as response:
            expected = f"bytes {start}-{start + count - 1}/{self.size}"
            if response.status != 206 or response.headers.get("Content-Range") != expected:
                raise ValueError("server did not honor the requested archive byte range")
            data = response.read(count + 1)
        if len(data) != count:
            raise ValueError("truncated or oversized archive range")
        self.fetched += count
        return data

    def prefetch(self, info: zipfile.ZipInfo) -> None:
        self.cache_start = info.header_offset
        self.cache = self.fetch(
            info.header_offset, 30 + len(info.filename.encode("utf-8")) + 65535 + info.compress_size
        )

    def read(self, size: int = -1) -> bytes:
        size = min(self.size - self.pos, size if size >= 0 else self.size - self.pos)
        if size <= 0:
            return b""
        offset = self.pos - self.cache_start
        if offset >= 0 and offset + size <= len(self.cache):
            data = self.cache[offset : offset + size]
        else:
            data = self.fetch(self.pos, size)
        self.pos += len(data)
        return data


def get_json(url: str) -> dict:
    with urllib.request.urlopen(url, timeout=60) as response:
        return json.load(response)


def fetch_selected(work: Path, source: str) -> dict:
    import hashlib

    if source == "morphology":
        metadata_url = "https://zenodo.org/api/records/6579584"
        metadata = get_json(metadata_url)
        archives = [metadata["files"][0]]
        license_value = metadata["metadata"]["license"]
    else:
        metadata_url = "https://api.figshare.com/v2/articles/5854692"
        metadata = get_json(metadata_url)
        archives = [
            f
            for f in metadata["files"]
            if f["name"] in {"4_orthomcl_output.zip", "1_molecular_phylogenetics.zip"}
        ]
        license_value = metadata["license"]
    directory = work / source
    directory.mkdir(parents=True, exist_ok=True)
    manifest = {"source": metadata_url, "license": license_value, "archives": []}
    for archive in archives:
        url = archive.get("download_url") or archive["links"]["self"]
        remote = RangeReader(url, archive["size"])
        record = {
            "name": archive.get("name", archive.get("key")),
            "size": archive["size"],
            "published_checksum": archive.get("computed_md5", archive.get("checksum")),
            "members": [],
        }
        with zipfile.ZipFile(remote) as zipped:
            for info in zipped.infolist():
                name = info.filename
                wanted = (
                    (
                        ("/morph/mrbayes/" in name and name.lower().endswith("morph_overlap.nex"))
                        or (
                            "/094_caviidae/morph/TNT/" in name and name.endswith("nexus_xinact.nex")
                        )
                    )
                    if source == "morphology"
                    else (
                        name.endswith("orthomcl_clusters.txt")
                        or name.endswith("orthomcl_SeqIDs_index.txt")
                        or name.endswith("/343_2408OGs.treefile")
                    )
                )
                if not wanted:
                    continue
                parts = PurePosixPath(name)
                if parts.is_absolute() or ".." in parts.parts:
                    raise ValueError("unsafe archive path")
                dest = directory / Path(*parts.parts)
                remote.prefetch(info)
                data = zipped.read(info)  # zipfile checks decompression and CRC
                dest.parent.mkdir(parents=True, exist_ok=True)
                dest.write_bytes(data)
                record["members"].append(
                    {
                        "name": name,
                        "size": len(data),
                        "crc32": f"{info.CRC:08x}",
                        "sha256": hashlib.sha256(data).hexdigest(),
                    }
                )
        record["downloaded_bytes"] = remote.fetched
        manifest["archives"].append(record)
    if source == "yeast":
        names = next(f for f in metadata["files"] if f["name"].startswith("343taxa_speicies"))
        with urllib.request.urlopen(names["download_url"], timeout=30) as response:
            data = response.read()
        (directory / "species_names.txt").write_bytes(data)
        manifest["species_names"] = {
            "name": names["name"],
            "sha256": hashlib.sha256(data).hexdigest(),
        }
    (directory / "source-manifest.json").write_text(json_text(manifest), encoding="utf-8")
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", choices=["yeast", "morphology"])
    parser.add_argument("--work", type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(fetch_selected(args.work, args.source), indent=2))


if __name__ == "__main__":
    main()
