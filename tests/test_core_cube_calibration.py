from pathlib import Path

import astropy.units as u
import numpy as np

import asheis.core as core
from asheis import asheis as public_asheis
from asheis.core import asheis


class FakeTemplate:
    central_wave = 195.119


class FakeCube:
    def __init__(self, source="counts"):
        self.source = source
        self.meta = {"calib_source": "Fake calibration"}
        self.smoothed = False

    def smooth_cube(self, rebin):
        self.smoothed = rebin
        return self


class FakeDate:
    def strftime(self, fmt):
        return "2024_01_01__00_00_00"


class FakeMap:
    measurement = "EIS Intensity"
    date = FakeDate()

    def __init__(self):
        self.data = np.full((2, 2), 2.0)
        self.meta = {}

    def save(self, path, overwrite=False):
        Path(path).parent.mkdir(parents=True, exist_ok=True)
        Path(path).write_text("fake fits", encoding="utf-8")


class FakeFitResult:
    def __init__(self):
        self.fit = {
            "wave_range": np.array([194.0, 196.0]),
            "err_int": np.full((2, 2, 1), 5.0),
        }
        self.meta = {
            "mod_index": {"crval2": 10.0},
            "slit_width": np.array([0.06, 0.06]),
        }

    def shift2wave(self, value, wave):
        return value

    def get_map(self, component, measurement):
        return FakeMap()

    def get_params(self, component, param_name):
        return np.full((2, 2), 195.12), np.full((2, 2), 0.01)


def make_asheis(tmp_path):
    obj = asheis.__new__(asheis)
    obj.filename = str(tmp_path / "eis_20240101_000000.data.h5")
    obj.ncpu = 1
    obj.rebin = False
    obj._aia_alignment_cache = {}
    obj.dict = {
        "fe_12_195.12": ["fe_12_195_119.2c.template.h5", 0, 6.2],
    }
    return obj


def patch_fit_dependencies(monkeypatch, read_cube_calls, calibrate_calls, save_dirs, ccd_offset=None):
    import eispac.instr

    monkeypatch.setattr(
        core.eispac.data,
        "get_fit_template_filepath",
        lambda template_name: template_name,
    )
    monkeypatch.setattr(core.eispac, "read_template", lambda template_path: FakeTemplate())

    def fake_read_cube(filename, window, apply_radcal=True):
        read_cube_calls.append(
            {
                "filename": filename,
                "window": window,
                "apply_radcal": apply_radcal,
            }
        )
        return FakeCube("preflight" if apply_radcal else "counts")

    def fake_calibrate_cube(cube, calibration):
        calibrate_calls.append(calibration)
        cube.meta["calib_source"] = f"calib {calibration}"
        return cube

    def fake_save_fit(fit_res, save_dir=None):
        save_dirs.append(Path(save_dir))
        return Path(save_dir) / "saved.fit.h5"

    monkeypatch.setattr(core.eispac, "read_cube", fake_read_cube)
    monkeypatch.setattr(core, "calibrate_cube", fake_calibrate_cube)
    monkeypatch.setattr(core.eispac, "fit_spectra", lambda cube, template, ncpu: FakeFitResult())
    monkeypatch.setattr(core.eispac, "save_fit", fake_save_fit)
    if ccd_offset is None:
        ccd_offset = lambda wavelength: 0.0 * u.pixel
    monkeypatch.setattr(eispac.instr, "ccd_offset", ccd_offset)


def test_public_package_import_exports_asheis_class():
    assert public_asheis is asheis


def test_default_fit_uses_2023_calibration(monkeypatch, tmp_path):
    obj = make_asheis(tmp_path)
    read_cube_calls = []
    calibrate_calls = []
    save_dirs = []
    patch_fit_dependencies(monkeypatch, read_cube_calls, calibrate_calls, save_dirs)

    obj.fit_data("fe_12_195.12", "int", False, tmp_path)

    assert read_cube_calls[0]["apply_radcal"] is False
    assert calibrate_calls == ["2023"]
    assert save_dirs[0].name == "fit_calib_2023"


def test_true_calibration_alias_uses_default_2023(monkeypatch, tmp_path):
    obj = make_asheis(tmp_path)
    read_cube_calls = []
    calibrate_calls = []
    save_dirs = []
    patch_fit_dependencies(monkeypatch, read_cube_calls, calibrate_calls, save_dirs)

    obj.fit_data("fe_12_195.12", "int", False, tmp_path, calib=True)

    assert read_cube_calls[0]["apply_radcal"] is False
    assert calibrate_calls == ["2023"]
    assert save_dirs[0].name == "fit_calib_2023"


def test_fit_accepts_eispac_numpy_ccd_offset(monkeypatch, tmp_path):
    obj = make_asheis(tmp_path)
    read_cube_calls = []
    calibrate_calls = []
    save_dirs = []
    patch_fit_dependencies(
        monkeypatch,
        read_cube_calls,
        calibrate_calls,
        save_dirs,
        ccd_offset=lambda wavelength: np.atleast_1d(float(np.asarray(wavelength).reshape(-1)[0])),
    )

    obj.fit_data("fe_12_195.12", "int", False, tmp_path, calib=2023)

    assert save_dirs[0].name == "fit_calib_2023"


def test_fit_passes_quantity_to_eispac_ccd_offset(monkeypatch, tmp_path):
    obj = make_asheis(tmp_path)
    read_cube_calls = []
    calibrate_calls = []
    save_dirs = []

    def strict_ccd_offset(wavelength):
        assert hasattr(wavelength, "unit")
        return 0.0 * u.pixel

    patch_fit_dependencies(
        monkeypatch,
        read_cube_calls,
        calibrate_calls,
        save_dirs,
        ccd_offset=strict_ccd_offset,
    )

    obj.fit_data("fe_12_195.12", "int", False, tmp_path, calib=2023)

    assert save_dirs[0].name == "fit_calib_2023"


def test_fit_falls_back_to_float_for_legacy_ccd_offset(monkeypatch, tmp_path):
    obj = make_asheis(tmp_path)
    read_cube_calls = []
    calibrate_calls = []
    save_dirs = []
    calls = []

    def legacy_ccd_offset(wavelength):
        calls.append(wavelength)
        if hasattr(wavelength, "unit"):
            raise u.UnitConversionError("legacy ccd_offset expects a bare wavelength")
        return np.atleast_1d(0.0)

    patch_fit_dependencies(
        monkeypatch,
        read_cube_calls,
        calibrate_calls,
        save_dirs,
        ccd_offset=legacy_ccd_offset,
    )

    obj.fit_data("fe_12_195.12", "int", False, tmp_path, calib=2023)

    assert any(hasattr(wavelength, "unit") for wavelength in calls)
    assert any(not hasattr(wavelength, "unit") for wavelength in calls)
    assert save_dirs[0].name == "fit_calib_2023"


def test_calibrated_fit_reads_counts_cube_and_uses_calibration_cache(monkeypatch, tmp_path):
    obj = make_asheis(tmp_path)
    read_cube_calls = []
    calibrate_calls = []
    save_dirs = []
    patch_fit_dependencies(monkeypatch, read_cube_calls, calibrate_calls, save_dirs)

    obj.fit_data("fe_12_195.12", "int", False, tmp_path, calib=2026)

    assert read_cube_calls[0]["apply_radcal"] is False
    assert calibrate_calls == ["2026"]
    assert save_dirs[0].name == "fit_calib_2026"


def test_preflight_fit_keeps_eispac_radcal_and_uses_preflight_cache(monkeypatch, tmp_path):
    obj = make_asheis(tmp_path)
    read_cube_calls = []
    calibrate_calls = []
    save_dirs = []
    patch_fit_dependencies(monkeypatch, read_cube_calls, calibrate_calls, save_dirs)

    obj.fit_data("fe_12_195.12", "int", False, tmp_path, calib=False)

    assert read_cube_calls[0]["apply_radcal"] is True
    assert calibrate_calls == []
    assert save_dirs[0].name == "fit_preflight"


def test_cached_fit_is_loaded_from_calibration_specific_path(monkeypatch, tmp_path):
    obj = make_asheis(tmp_path)
    path = obj._fit_cache_path("fe_12_195_119.2c.template.h5", "fe_12_195.12", 2023)
    path.parent.mkdir(parents=True)
    path.write_text("cached", encoding="utf-8")
    read_fit_calls = []

    monkeypatch.setattr(core.eispac, "read_fit", lambda fit_path: read_fit_calls.append(Path(fit_path)) or FakeFitResult())
    monkeypatch.setattr(
        core.eispac,
        "read_cube",
        lambda *args, **kwargs: (_ for _ in ()).throw(AssertionError("read_cube should not be called")),
    )

    obj.fit_data("fe_12_195.12", "int", False, tmp_path, calib=2023)

    assert read_fit_calls == [path]
    assert path.parent.name == "fit_calib_2023"


def test_get_intensity_mcmc_returns_fit_error_without_second_multiplier(tmp_path):
    obj = make_asheis(tmp_path)
    obj.fit_data = lambda *args, **kwargs: FakeFitResult()
    obj.directory_setup = lambda amap, line, outdir: "2024_01_01__00_00_00"

    data, error = obj.get_intensity(
        "fe_12_195.12",
        outdir=tmp_path,
        plot=False,
        mcmc=True,
        calib=2014,
    )

    np.testing.assert_allclose(data, np.full((2, 2), 2.0))
    np.testing.assert_allclose(error, np.full((2, 2), 5.0))


def test_get_width_uses_fits_safe_metadata_keys(tmp_path):
    obj = make_asheis(tmp_path)
    obj.fit_data = lambda *args, **kwargs: FakeFitResult()
    obj.directory_setup = lambda amap, line, outdir: "2024_01_01__00_00_00"

    width_map = obj.get_width(
        "fe_12_195.12",
        outdir=tmp_path,
        plot=False,
        width_only=True,
    )

    assert "slitwid" in width_map.meta
    assert "centerr" in width_map.meta
    assert "slit_width" not in width_map.meta
    assert "cent_error" not in width_map.meta
