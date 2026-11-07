"""Small synthetic fixtures validate extraction without third-party test data."""

import io
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
import zipfile

from scripts.source_archives import RangeReader
from scripts.calibrate import cldf_matrix, compact_summary
from scripts.fit_presets import coverage_fit


class CalibrationToolTests(unittest.TestCase):
    def test_range_zip_and_budget(self):
        buffer = io.BytesIO()
        with zipfile.ZipFile(buffer, "w", zipfile.ZIP_DEFLATED) as z:
            z.writestr("matrix.nex", b"fixture " * 100)
        payload = buffer.getvalue()

        def request(req, timeout):
            start, end = map(int, req.get_header("Range").split("=")[1].split("-"))
            response = io.BytesIO(payload[start : end + 1])
            response.status = 206
            response.headers = {"Content-Range": f"bytes {start}-{end}/{len(payload)}"}
            return response

        with patch("urllib.request.urlopen", side_effect=request):
            remote = RangeReader("https://example.invalid/source.zip", len(payload))
            with zipfile.ZipFile(remote) as archive:
                info = archive.getinfo("matrix.nex")
                remote.prefetch(info)
                self.assertEqual(archive.read(info), b"fixture " * 100)
            with self.assertRaises(ValueError):
                RangeReader("https://example.invalid/source.zip", len(payload), budget=1).fetch(
                    0, 2
                )

    def test_reject_ignored_range(self):
        response = io.BytesIO(b"whole file")
        response.status = 200
        response.headers = {}
        with patch("urllib.request.urlopen", return_value=response), self.assertRaises(ValueError):
            RangeReader("https://example.invalid/archive", 10).read(2)

    def test_cldf_missing_and_partial_cognates(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp)
            (path / "forms.csv").write_text(
                "ID,Language_ID,Parameter_ID\nf1,A,g1\nf2,B,g2\nf3,A,g2\n"
            )
            (path / "cognates.csv").write_text(
                "ID,Form_ID,Cognateset_ID\nc1,f1,x\nc2,f1,y\nc3,f2,z\n"
            )
            matrix = cldf_matrix(path)
            self.assertEqual(len(matrix["A"][0]), 2)
            self.assertIsNone(matrix["A"][1])
            self.assertIsNone(matrix["B"][0])

    def test_beta_coverage_moments(self):
        matrix = {"a": [frozenset({0})] * 4, "b": [None] * 2 + [frozenset({0})] * 2}
        summary = compact_summary(matrix)
        params = coverage_fit([{"summary": summary}])
        a, b = params["coverage_a"], params["coverage_b"]
        self.assertAlmostEqual(a / (a + b), 0.75)
        self.assertAlmostEqual(a * b / ((a + b) ** 2 * (a + b + 1)), 0.0625)
