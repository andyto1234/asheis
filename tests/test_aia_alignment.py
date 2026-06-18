from pathlib import Path

import astropy.units as u
from astropy.time import Time

from asheis.aia_alignment import (
    AIAAlignment,
    apply_aia_alignment,
    expected_aia_filename,
    find_nearest_cached_aia,
    require_jsoc_email,
)
from asheis.core import asheis


class DummyDate:
    def strftime(self, fmt):
        return "2024_01_02__03_04_05"


class DummyMap:
    measurement = "EIS Intensity"
    date = DummyDate()

    def __init__(self, meta=None):
        self.meta = dict(meta or {})
        self.shift = None
        self.saved_path = None

    def shift_reference_coord(self, dx, dy):
        shifted = DummyMap(self.meta)
        shifted.shift = (dx.to_value(u.arcsec), dy.to_value(u.arcsec))
        return shifted

    def save(self, path, overwrite=False):
        self.saved_path = Path(path)
        self.saved_path.write_text("dummy fits", encoding="utf-8")


def test_find_nearest_cached_aia_uses_tolerance(tmp_path):
    target = Time("2024-01-02T03:04:05", scale="utc")
    near = tmp_path / expected_aia_filename(Time("2024-01-02T03:04:17", scale="utc"))
    far = tmp_path / expected_aia_filename(Time("2024-01-02T03:09:05", scale="utc"))
    near.write_text("near", encoding="utf-8")
    far.write_text("far", encoding="utf-8")

    assert find_nearest_cached_aia(tmp_path, target, max_offset_sec=20) == near
    assert find_nearest_cached_aia(tmp_path, target, max_offset_sec=5) is None


def test_require_jsoc_email_uses_env(monkeypatch):
    monkeypatch.setenv("JSOC_EMAIL", "person@example.com")
    assert require_jsoc_email(None) == "person@example.com"


def test_apply_aia_alignment_updates_metadata_and_saves(tmp_path):
    alignment = AIAAlignment(
        dx_arcsec=1.5,
        dy_arcsec=-2.5,
        corr_max=0.91,
        x_shift_ref_pix=3.0,
        y_shift_ref_pix=4.0,
        eis_time="2024-01-02T03:04:05",
        aia_time="2024-01-02T03:04:17",
        aia_fits="/tmp/aia.fits",
        qa_plot="/tmp/qa.png",
    )

    aligned = apply_aia_alignment(
        DummyMap(),
        alignment,
        line="fe_16_262.98",
        outdir=tmp_path,
        overwrite=True,
    )

    assert aligned.shift == (1.5, -2.5)
    assert aligned.meta["aia_dx"] == 1.5
    assert aligned.meta["aia_dy"] == -2.5
    assert aligned.meta["aia_corr"] == 0.91
    assert aligned.meta["aia_ref_line"] == "fe_12_195.12"
    assert Path(aligned.meta["aia_aligned_fits"]).is_file()


def test_asheis_alignment_cache_reuses_computed_shift(monkeypatch, tmp_path):
    obj = asheis.__new__(asheis)
    obj.filename = "eis_l0_20240102_030405.data.h5"
    obj._aia_alignment_cache = {}
    target = DummyMap()

    calls = {"count": 0}
    alignment = AIAAlignment(
        dx_arcsec=1.0,
        dy_arcsec=2.0,
        corr_max=0.8,
        x_shift_ref_pix=0.0,
        y_shift_ref_pix=0.0,
        eis_time="2024-01-02T03:04:05",
        aia_time="2024-01-02T03:04:17",
        aia_fits="/tmp/aia.fits",
    )

    def fake_compute(*args, **kwargs):
        calls["count"] += 1
        return alignment

    monkeypatch.setattr("asheis.aia_alignment.compute_aia_alignment", fake_compute)

    first = obj._get_aia_alignment(
        target,
        "fe_12_195.12",
        tmp_path,
        False,
        2014,
        None,
        None,
        "person@example.com",
        60.0,
        False,
        False,
    )
    second = obj._get_aia_alignment(
        target,
        "fe_12_195.12",
        tmp_path,
        False,
        2014,
        None,
        None,
        "person@example.com",
        60.0,
        False,
        False,
    )

    assert first is second
    assert calls["count"] == 1

