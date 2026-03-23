"""
Unit tests for bin/MutationsDictionary.py — position translation logic.
All tests use tmp_path and subprocess.run; no bioinformatics tools required.
"""
import csv
import subprocess
import sys
from pathlib import Path

# Path to the script under test
SCRIPT = Path(__file__).parents[2] / "nf_pipeline" / "bin" / "MutationsDictionary.py"
DICT_FIXTURE = Path(__file__).parents[1] / "fixtures" / "dictionary_mini.csv"


def _write_input(tmp_path, rows):
    """Write a minimal input CSV with a POSITION column."""
    p = tmp_path / "input.csv"
    with open(p, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["POSITION", "AA"])
        writer.writeheader()
        for row in rows:
            writer.writerow(row)
    return p


def _run(tmp_path, subtype, input_csv, base="H5"):
    output = tmp_path / "output.csv"
    result = subprocess.run(
        [
            sys.executable, str(SCRIPT),
            "--subtype", subtype,
            "--input", str(input_csv),
            "--dictionary", str(DICT_FIXTURE),
            "--base", base,
            "--output", str(output),
        ],
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, f"Script failed:\n{result.stderr}"
    return output


def _positions(output_csv):
    with open(output_csv) as f:
        return [row["POSITION"] for row in csv.DictReader(f)]


def test_h5_to_h3_known_position(tmp_path):
    inp = _write_input(tmp_path, [{"POSITION": "182", "AA": "S"}])
    out = _run(tmp_path, "H3", inp)
    assert _positions(out) == ["157"]


def test_h5_to_h1_known_position(tmp_path):
    inp = _write_input(tmp_path, [{"POSITION": "182", "AA": "S"}])
    out = _run(tmp_path, "H1", inp)
    assert _positions(out) == ["157"]


def test_h5_to_h7_known_position(tmp_path):
    inp = _write_input(tmp_path, [{"POSITION": "182", "AA": "S"}])
    out = _run(tmp_path, "H7", inp)
    assert _positions(out) == ["182"]


def test_h5_to_h9_known_position(tmp_path):
    inp = _write_input(tmp_path, [{"POSITION": "226", "AA": "Q"}])
    out = _run(tmp_path, "H9", inp)
    assert _positions(out) == ["226"]


def test_position_not_in_dictionary(tmp_path):
    inp = _write_input(tmp_path, [{"POSITION": "999", "AA": "X"}])
    out = _run(tmp_path, "H3", inp)
    assert _positions(out) == ["-"]


def test_whitespace_in_position(tmp_path):
    # Write raw CSV with leading/trailing space in POSITION value
    p = tmp_path / "input.csv"
    p.write_text("POSITION,AA\n 182 ,S\n")
    out = _run(tmp_path, "H3", p)
    assert _positions(out) == ["157"]


def test_all_positions_translated(tmp_path):
    rows = [
        {"POSITION": "100", "AA": "K"},
        {"POSITION": "182", "AA": "S"},
        {"POSITION": "226", "AA": "Q"},
        {"POSITION": "300", "AA": "T"},
    ]
    inp = _write_input(tmp_path, rows)
    out = _run(tmp_path, "H3", inp)
    assert _positions(out) == ["80", "157", "201", "275"]


def test_empty_input_produces_empty_output(tmp_path):
    p = tmp_path / "input.csv"
    p.write_text("POSITION,AA\n")
    out = _run(tmp_path, "H3", p)
    with open(out) as f:
        rows = list(csv.DictReader(f))
    assert rows == []
