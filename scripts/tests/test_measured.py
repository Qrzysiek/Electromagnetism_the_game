"""The documents' measured-value markers (DEC-22): each closed, each holding a number alone,
so that `scripts/measured.py update` replaces exactly the number."""
import re

import measured


def test_every_marker_is_closed_and_holds_a_number():
    number = re.compile(r"^−?\d+(\.\d+)?(e−?-?\d+)?$")
    for path in measured.documents():
        text = open(path, encoding="utf-8").read()
        assert text.count("<!--m:") == text.count("<!--/m-->"), path
        for key, value in measured.quotes(text):
            assert number.match(value), f"{path}: {key}: {value!r}"


def test_a_negative_record_is_typeset_with_a_minus_sign(tmp_path, monkeypatch):
    monkeypatch.setattr(measured, "RECORDS", str(tmp_path))
    (tmp_path / "x").write_text("-0.3", encoding="utf-8")
    (tmp_path / "y").write_text("1.2e-8", encoding="utf-8")
    assert measured.recorded("x") == "−0.3"
    assert measured.recorded("y") == "1.2e-8"
    assert measured.recorded("absent") is None
