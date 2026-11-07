"""Measure the specified 1,000 x 50-taxon x 200-concept grid on one core."""

from pathlib import Path
import argparse
import json
import platform
import time

from ngesh.__main__ import run_grid


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--replicates", type=int, default=1000)
    args = parser.parse_args()
    start = time.perf_counter()
    cpu = time.process_time()
    result = run_grid(
        {
            "seed": "ngesh-pilot-benchmark-v1",
            "replicates": args.replicates,
            "preset": "lexical",
            "base": {"tree": {"n_taxa": 50}, "characters": {"n_characters": 200}},
        },
        args.out,
    )
    print(
        json.dumps(
            dict(
                result,
                wall_seconds=time.perf_counter() - start,
                cpu_seconds=time.process_time() - cpu,
                python=platform.python_version(),
                platform=platform.platform(),
            ),
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
