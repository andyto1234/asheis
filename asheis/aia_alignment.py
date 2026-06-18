"""Align EIS maps to AIA 193 by correcting the EIS reference coordinate."""

from __future__ import annotations

import csv
import os
import re
import socket
import time
import warnings
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.request import urlopen

import astropy.units as u
import numpy as np
from astropy.table import QTable
from astropy.time import Time


AIA_START_TIME = Time("2010-05-05T11:17:34", scale="utc")
AIA_SERIES = "aia.lev1_euv_12s"
AIA_WAVELENGTH = 193 * u.angstrom
AIA_FILENAME_TIME_RE = re.compile(
    r"(?P<time>\d{4}-\d{2}-\d{2}T\d{6})Z.*193.*\.fits(?:\.gz)?$",
    re.IGNORECASE,
)
DRMS_TIME_RE = re.compile(
    r"(?P<date>\d{4}[.-]\d{2}[.-]\d{2})[_T](?P<time>\d{2}:\d{2}:\d{2})"
)


@dataclass(frozen=True)
class AIAAlignment:
    """Computed pointing correction from EIS Fe XII 195 to AIA 193."""

    dx_arcsec: float
    dy_arcsec: float
    corr_max: float
    x_shift_ref_pix: float
    y_shift_ref_pix: float
    eis_time: str
    aia_time: str
    aia_fits: str
    aia_ref_line: str = "fe_12_195.12"
    qa_plot: str = ""


def map_time(amap) -> Time:
    """Return a SunPy map date as an Astropy UTC Time."""
    return Time(amap.date).utc


def default_aia_cache_dir(outdir: str | os.PathLike[str]) -> Path:
    return Path(outdir).expanduser() / "aia_cache"


def aia_time_key(time_value: Time) -> str:
    return Time(time_value).utc.strftime("%Y-%m-%dT%H%M%SZ")


def expected_aia_filename(time_value: Time) -> str:
    return f"aia.lev1_euv_12s.{aia_time_key(time_value)}.193.image_lev1.fits"


def parse_aia_time_from_filename(path: Path) -> Time | None:
    match = AIA_FILENAME_TIME_RE.search(path.name)
    if not match:
        return None
    try:
        parsed = datetime.strptime(match.group("time"), "%Y-%m-%dT%H%M%S")
        return Time(parsed.replace(tzinfo=timezone.utc))
    except ValueError:
        return None


def seconds_between(left: Time, right: Time) -> float:
    return abs((Time(left).utc - Time(right).utc).to_value(u.s))


def find_nearest_cached_aia(
    cache_dir: Path,
    eis_time: Time,
    *,
    max_offset_sec: float,
) -> Path | None:
    """Find the nearest cached AIA 193 FITS file within a time tolerance."""
    if not cache_dir.exists():
        return None

    best: tuple[float, Path] | None = None
    for path in cache_dir.glob("*.fits*"):
        if not path.is_file() or path.stat().st_size <= 0:
            continue
        aia_time = parse_aia_time_from_filename(path)
        if aia_time is None:
            continue
        offset = seconds_between(aia_time, eis_time)
        if offset <= max_offset_sec and (best is None or offset < best[0]):
            best = (offset, path)
    return best[1] if best is not None else None


def require_jsoc_email(jsoc_email: str | None) -> str:
    email = jsoc_email or os.environ.get("JSOC_EMAIL")
    if not email:
        raise ValueError(
            "AIA download requires a registered JSOC email. Pass jsoc_email=... "
            "or set the JSOC_EMAIL environment variable."
        )
    return email


def validate_jsoc_email(email: str) -> None:
    import drms

    client = drms.Client()
    if not client.check_email(email):
        raise ValueError(
            f"JSOC email is not registered: {email}. Register it at "
            "http://jsoc.stanford.edu/ajax/register_email.html"
        )


def nearest_index(sorted_times: list[Time], target: Time) -> int:
    offsets = [seconds_between(item, target) for item in sorted_times]
    return int(np.argmin(offsets))


def parse_jsoc_trec(value: object) -> Time:
    text = str(value)
    match = DRMS_TIME_RE.search(text)
    if match:
        stamp = f"{match.group('date').replace('.', '-')} {match.group('time')}"
        parsed = datetime.strptime(stamp, "%Y-%m-%d %H:%M:%S")
        return Time(parsed.replace(tzinfo=timezone.utc))
    return Time(text, scale="utc")


def query_nearest_aia_record(
    eis_time: Time,
    *,
    jsoc_email: str,
    max_offset_sec: float,
):
    """Return the nearest JSOC AIA 193 row to an EIS observation time."""
    from sunpy.net import Fido, attrs as a

    start = eis_time - max_offset_sec * u.s
    end = eis_time + max_offset_sec * u.s
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        result = Fido.search(
            a.Time(start.to_datetime(), end.to_datetime()),
            a.jsoc.Series(AIA_SERIES),
            a.jsoc.Segment("image"),
            a.Wavelength(AIA_WAVELENGTH),
            a.jsoc.Notify(jsoc_email),
        )

    if len(result) == 0 or len(result[0]) == 0:
        raise FileNotFoundError(
            f"No AIA 193 JSOC record found within {max_offset_sec:.1f}s of {eis_time.isot}"
        )

    table = result[0]
    trec = [parse_jsoc_trec(value) for value in table["T_REC"]]
    index = nearest_index(trec, eis_time)
    offset = seconds_between(trec[index], eis_time)
    if offset > max_offset_sec:
        raise FileNotFoundError(
            f"Nearest AIA 193 record is {offset:.1f}s from EIS time, "
            f"exceeding max_aia_offset_sec={max_offset_sec:.1f}"
        )
    return table[[index]], trec[index]


def drms_time_label(time_value: Time) -> str:
    return Time(time_value).utc.strftime("%Y.%m.%d_%H:%M:%S_UTC")


def drms_aia_record(time_value: Time) -> str:
    return f"{AIA_SERIES}[{drms_time_label(time_value)}][193]{{image}}"


def wait_for_drms_export(request, *, poll_sec: float = 30.0, timeout_sec: float = 900.0) -> None:
    started = time.time()
    while True:
        if request.has_finished():
            if request.has_failed(skip_update=True):
                request._raise_on_error()
            return
        if time.time() - started > timeout_sec:
            raise TimeoutError("Timed out waiting for JSOC DRMS export")
        time.sleep(poll_sec)


def download_url_to_path(
    url: str,
    dest: Path,
    *,
    overwrite: bool,
    timeout_sec: float = 300.0,
    retries: int = 3,
    retry_wait_sec: float = 5.0,
) -> Path:
    if dest.is_file() and dest.stat().st_size > 0 and not overwrite:
        return dest

    dest.parent.mkdir(parents=True, exist_ok=True)
    last_error: Exception | None = None
    for attempt in range(1, retries + 1):
        tmp = dest.with_name(f".{dest.name}.{time.time_ns()}.part")
        try:
            with urlopen(url, timeout=timeout_sec) as response:
                with tmp.open("wb") as handle:
                    while True:
                        chunk = response.read(1024 * 1024)
                        if not chunk:
                            break
                        handle.write(chunk)
            tmp.replace(dest)
            return dest
        except (HTTPError, URLError, OSError, TimeoutError, socket.timeout) as exc:
            last_error = exc
            tmp.unlink(missing_ok=True)
            if attempt < retries:
                time.sleep(retry_wait_sec * attempt)

    if last_error is not None:
        raise last_error
    raise RuntimeError(f"Could not download {url}")


def download_aia_with_drms(
    aia_time: Time,
    cache_dir: Path,
    *,
    jsoc_email: str,
    overwrite: bool,
) -> Path:
    """Download an exact AIA 193 image segment through DRMS export URLs."""
    import drms

    expected = cache_dir / expected_aia_filename(aia_time)
    if expected.is_file() and expected.stat().st_size > 0 and not overwrite:
        return expected

    client = drms.Client(email=jsoc_email)
    request = client.export(
        drms_aia_record(aia_time),
        method="url_quick",
        protocol="as-is",
        email=jsoc_email,
    )
    wait_for_drms_export(request)
    urls = request.urls.reset_index(drop=True)
    if len(urls) == 0:
        raise FileNotFoundError(f"DRMS export returned no URL for {aia_time.isot}")

    path = download_url_to_path(
        str(urls.loc[0, "url"]),
        expected,
        overwrite=overwrite,
    )
    # Validate early so truncated or HTML error files do not enter the cache.
    import sunpy.map

    sunpy.map.Map(path)
    return path


def resolve_aia_fits(
    eis_time: Time,
    *,
    aia_fits: str | os.PathLike[str] | None,
    cache_dir: Path,
    jsoc_email: str | None,
    max_offset_sec: float,
    overwrite: bool,
) -> Path:
    """Resolve an AIA FITS path from explicit input, local cache, or JSOC."""
    if aia_fits is not None:
        path = Path(aia_fits).expanduser()
        if not path.is_file() or path.stat().st_size <= 0:
            raise FileNotFoundError(f"AIA FITS not found or empty: {path}")
        return path

    cached = None if overwrite else find_nearest_cached_aia(
        cache_dir,
        eis_time,
        max_offset_sec=max_offset_sec,
    )
    if cached is not None:
        return cached

    if eis_time < AIA_START_TIME:
        raise ValueError(f"EIS observation is before AIA availability: {eis_time.isot}")

    email = require_jsoc_email(jsoc_email)
    validate_jsoc_email(email)
    _record, aia_time = query_nearest_aia_record(
        eis_time,
        jsoc_email=email,
        max_offset_sec=max_offset_sec,
    )
    return download_aia_with_drms(
        aia_time,
        cache_dir,
        jsoc_email=email,
        overwrite=overwrite,
    )


def sanitize_pointing_table_for_cache(table: QTable) -> QTable:
    cached = QTable(table, copy=True)
    cached.meta.clear()
    for column in cached.itercols():
        if column.info.meta is not None:
            column.info.meta = {}
    return cached


def load_or_fetch_pointing_table(
    aia_time: Time,
    table_path: Path,
    *,
    margin_hours: float = 12.0,
) -> QTable:
    """Load a cached AIA pointing table, or fetch one that covers the AIA time."""
    from aiapy.calibrate.util import get_pointing_table

    start = aia_time - margin_hours * u.hour
    end = aia_time + margin_hours * u.hour
    if table_path.exists():
        try:
            table = QTable.read(table_path, format="ascii.ecsv")
        except Exception:
            table = None
        if table is not None and len(table) > 0:
            table_start = Time(table["T_START"][0])
            table_end = Time(table["T_STOP"][-1])
            if table_start <= start and table_end >= end:
                return table

    table = get_pointing_table(start, end)
    if len(table) == 0:
        raise RuntimeError(f"No AIA pointing rows found for {start.isot} -> {end.isot}")

    table = sanitize_pointing_table_for_cache(table)
    table_path.parent.mkdir(parents=True, exist_ok=True)
    table.write(table_path, format="ascii.ecsv", overwrite=True)
    return table


def finite_float_array(data: np.ndarray) -> np.ndarray:
    array = np.asarray(data, dtype=np.float64)
    finite = np.isfinite(array)
    if finite.all():
        return array
    fill = np.nanmedian(array[finite]) if finite.any() else 0.0
    return np.nan_to_num(array, nan=fill, posinf=fill, neginf=fill)


def preprocess_for_matching(data: np.ndarray) -> np.ndarray:
    array = finite_float_array(data)
    floor = np.nanpercentile(array, 1)
    array = np.clip(array - floor, 0, None)
    return np.sqrt(array)


def submap_aia_around_eis(aia_map, eis_map, margin_arcsec: float):
    bottom_left = eis_map.bottom_left_coord
    top_right = eis_map.top_right_coord
    bottom_left = bottom_left.__class__(
        bottom_left.Tx - margin_arcsec * u.arcsec,
        bottom_left.Ty - margin_arcsec * u.arcsec,
        frame=bottom_left.frame,
    )
    top_right = top_right.__class__(
        top_right.Tx + margin_arcsec * u.arcsec,
        top_right.Ty + margin_arcsec * u.arcsec,
        frame=top_right.frame,
    )
    return aia_map.submap(bottom_left, top_right=top_right)


def resample_aia_to_eis_scale(aia_map, eis_map):
    nx = (aia_map.scale.axis1 * aia_map.dimensions.x) / eis_map.scale.axis1
    ny = (aia_map.scale.axis2 * aia_map.dimensions.y) / eis_map.scale.axis2
    dims = u.Quantity(
        [max(int(np.round(nx.value)), 2), max(int(np.round(ny.value)), 2)],
        u.pix,
    )
    return aia_map.resample(dims)


def map_is_large_enough_for_template(image_map, template_map) -> bool:
    image_shape = np.asarray(image_map.data).shape
    template_shape = np.asarray(template_map.data).shape
    return image_shape[0] >= template_shape[0] and image_shape[1] >= template_shape[1]


def matching_aia_map(aia_map, eis_map, *, margin_arcsec: float):
    """Return an AIA map on the EIS plate scale that is large enough to match."""
    margins = [margin_arcsec, margin_arcsec * 2, margin_arcsec * 4, margin_arcsec * 8]
    last_map = None
    for margin in margins:
        aia_submap = submap_aia_around_eis(aia_map, eis_map, margin)
        candidate = resample_aia_to_eis_scale(aia_submap, eis_map)
        last_map = candidate
        if map_is_large_enough_for_template(candidate, eis_map):
            return candidate, aia_submap

    full_disk_candidate = resample_aia_to_eis_scale(aia_map, eis_map)
    if map_is_large_enough_for_template(full_disk_candidate, eis_map):
        return full_disk_candidate, aia_map

    image_shape = np.asarray(last_map.data).shape if last_map is not None else None
    raise ValueError(
        "AIA match image is smaller than the EIS template after expanding the cutout "
        f"(AIA shape={image_shape}, EIS shape={np.asarray(eis_map.data).shape})."
    )


def parabolic_turning_point(values: np.ndarray) -> float:
    denominator = values.dot([1, -2, 1])
    if denominator == 0:
        return 1.0
    return (-0.5 * values.dot([-1, 0, 1])) / denominator


def best_match_location(corr: np.ndarray) -> tuple[float, float]:
    y0, x0 = np.unravel_index(np.argmax(corr), corr.shape)
    y_min = max(0, y0 - 1)
    y_max = min(y0 + 2, corr.shape[0])
    x_min = max(0, x0 - 1)
    x_max = min(x0 + 2, corr.shape[1])
    local = corr[y_min:y_max, x_min:x_max]
    local_y, local_x = np.unravel_index(np.argmax(local), local.shape)
    y_peak = parabolic_turning_point(local[:, local_x]) if local.shape[0] == 3 else float(local_y)
    x_peak = parabolic_turning_point(local[local_y, :]) if local.shape[1] == 3 else float(local_x)
    return y_min + y_peak, x_min + x_peak


def coalign_eis_to_aia(eis_map, aia_map) -> tuple[object, float, float, float, float, float]:
    from skimage.feature import match_template

    corr = match_template(
        preprocess_for_matching(aia_map.data),
        preprocess_for_matching(eis_map.data),
    )
    if corr.size < 2:
        raise ValueError("match_template returned too few correlation values")

    y_shift_ref_pix, x_shift_ref_pix = best_match_location(corr)
    corr_max = float(np.nanmax(corr))
    eis_ref_pix = u.Quantity(eis_map.reference_pixel).to_value(u.pix)
    matched_x = x_shift_ref_pix + eis_ref_pix[0]
    matched_y = y_shift_ref_pix + eis_ref_pix[1]
    new_reference_coordinate = aia_map.wcs.pixel_to_world(matched_x, matched_y)
    dx = (new_reference_coordinate.Tx - eis_map.reference_coordinate.Tx).to_value(u.arcsec)
    dy = (new_reference_coordinate.Ty - eis_map.reference_coordinate.Ty).to_value(u.arcsec)
    corrected = eis_map.shift_reference_coord(dx * u.arcsec, dy * u.arcsec)
    return corrected, dx, dy, x_shift_ref_pix, y_shift_ref_pix, corr_max


def prepare_aia_map(aia_fits: Path, cache_dir: Path):
    from aiapy.calibrate import register, update_pointing
    import sunpy.map

    aia_map = sunpy.map.Map(aia_fits)
    pointing_table = load_or_fetch_pointing_table(
        map_time(aia_map),
        cache_dir / "aia_pointing_table.ecsv",
    )
    aia_map = update_pointing(aia_map, pointing_table=pointing_table)
    return register(aia_map)


def alignment_output_dir(amap, line: str, outdir: str | os.PathLike[str]) -> Path:
    measurement = amap.measurement.lower().split()[-1]
    return Path(outdir).expanduser() / "images" / measurement / line


def aligned_fits_path(amap, line: str, outdir: str | os.PathLike[str]) -> Path:
    date = amap.date.strftime("%Y_%m_%d__%H_%M_%S")
    measurement = amap.measurement.lower().replace(" ", "_").replace(".", "_")
    return alignment_output_dir(amap, line, outdir) / f"eis_{date}_{measurement}_aia_aligned.fits"


def qa_plot_path(reference_map, outdir: str | os.PathLike[str]) -> Path:
    date = reference_map.date.strftime("%Y_%m_%d__%H_%M_%S")
    return (
        Path(outdir).expanduser()
        / "images"
        / "intensity"
        / "fe_12_195.12"
        / f"eis_{date}_fe_12_195_12_aia_alignment_qa.png"
    )


def draw_alignment_qa(
    aia_map,
    original_eis_map,
    corrected_eis_map,
    qa_path: Path,
    *,
    dx_arcsec: float,
    dy_arcsec: float,
    corr_max: float,
) -> Path:
    import matplotlib.pyplot as plt

    qa_path.parent.mkdir(parents=True, exist_ok=True)
    fig = plt.figure(figsize=(10, 4.5), dpi=160)
    axes = [
        fig.add_subplot(1, 2, 1, projection=aia_map),
        fig.add_subplot(1, 2, 2, projection=aia_map),
    ]
    titles = ["Original EIS contours", "AIA-aligned EIS contours"]
    contour_maps = [original_eis_map, corrected_eis_map]

    for ax, title, contour_map in zip(axes, titles, contour_maps):
        aia_map.plot(axes=ax, clip_interval=(1, 99.8) * u.percent)
        contour_map.draw_contours(
            levels=[60, 75, 90, 97] * u.percent,
            axes=ax,
            colors="white",
            linewidths=0.8,
        )
        ax.set_title(title)
        ax.coords[0].set_axislabel("Solar X")
        ax.coords[1].set_axislabel("Solar Y")

    fig.suptitle(f"dx={dx_arcsec:.2f}\"  dy={dy_arcsec:.2f}\"  corr={corr_max:.3f}")
    fig.tight_layout()
    fig.savefig(qa_path, bbox_inches="tight")
    plt.close(fig)
    return qa_path


def write_manifest_row(manifest_csv: Path, row: dict[str, object]) -> None:
    fields = [
        "eis_time",
        "aia_time",
        "aia_fits",
        "dx_arcsec",
        "dy_arcsec",
        "corr_max",
        "qa_plot",
    ]
    rows: list[dict[str, object]] = []
    if manifest_csv.exists():
        with manifest_csv.open(newline="", encoding="utf-8") as handle:
            rows.extend(csv.DictReader(handle))
    rows.append({field: row.get(field, "") for field in fields})
    manifest_csv.parent.mkdir(parents=True, exist_ok=True)
    with manifest_csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def compute_aia_alignment(
    reference_eis_map,
    *,
    outdir: str | os.PathLike[str],
    aia_fits: str | os.PathLike[str] | None = None,
    aia_cache_dir: str | os.PathLike[str] | None = None,
    jsoc_email: str | None = None,
    max_aia_offset_sec: float = 60.0,
    overwrite: bool = False,
    qa_plot: bool = True,
    aia_margin_arcsec: float = 80.0,
) -> AIAAlignment:
    """Compute an EIS-to-AIA pointing correction from Fe XII 195 intensity."""
    from sunpy.coordinates import SphericalScreen

    cache_dir = (
        Path(aia_cache_dir).expanduser()
        if aia_cache_dir is not None
        else default_aia_cache_dir(outdir)
    )
    eis_time = map_time(reference_eis_map)
    aia_path = resolve_aia_fits(
        eis_time,
        aia_fits=aia_fits,
        cache_dir=cache_dir,
        jsoc_email=jsoc_email,
        max_offset_sec=max_aia_offset_sec,
        overwrite=overwrite,
    )

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        aia_map = prepare_aia_map(aia_path, cache_dir)
        with SphericalScreen(aia_map.observer_coordinate, only_off_disk=True):
            aia_match_map, aia_plot_map = matching_aia_map(
                aia_map,
                reference_eis_map,
                margin_arcsec=aia_margin_arcsec,
            )
        corrected_map, dx, dy, x_shift, y_shift, corr_max = coalign_eis_to_aia(
            reference_eis_map,
            aia_match_map,
        )

    qa_path = ""
    if qa_plot:
        qa_path = str(
            draw_alignment_qa(
                aia_plot_map,
                reference_eis_map,
                corrected_map,
                qa_plot_path(reference_eis_map, outdir),
                dx_arcsec=dx,
                dy_arcsec=dy,
                corr_max=corr_max,
            )
        )

    alignment = AIAAlignment(
        dx_arcsec=dx,
        dy_arcsec=dy,
        corr_max=corr_max,
        x_shift_ref_pix=x_shift,
        y_shift_ref_pix=y_shift,
        eis_time=eis_time.isot,
        aia_time=map_time(aia_map).isot,
        aia_fits=str(aia_path),
        qa_plot=qa_path,
    )
    write_manifest_row(
        cache_dir / "eis_aia193_manifest.csv",
        {
            "eis_time": alignment.eis_time,
            "aia_time": alignment.aia_time,
            "aia_fits": alignment.aia_fits,
            "dx_arcsec": f"{alignment.dx_arcsec:.6f}",
            "dy_arcsec": f"{alignment.dy_arcsec:.6f}",
            "corr_max": f"{alignment.corr_max:.6f}",
            "qa_plot": alignment.qa_plot,
        },
    )
    return alignment


def apply_aia_alignment(
    target_map,
    alignment: AIAAlignment,
    *,
    line: str,
    outdir: str | os.PathLike[str],
    overwrite: bool = False,
):
    """Apply a precomputed AIA alignment shift to any EIS SunPy map."""
    aligned = target_map.shift_reference_coord(
        alignment.dx_arcsec * u.arcsec,
        alignment.dy_arcsec * u.arcsec,
    )
    aligned.meta["AIADX"] = float(alignment.dx_arcsec)
    aligned.meta["AIADY"] = float(alignment.dy_arcsec)
    aligned.meta["AIACORR"] = float(alignment.corr_max)
    aligned.meta["AIAFITS"] = alignment.aia_fits
    aligned.meta["AIATIME"] = alignment.aia_time
    aligned.meta["AIAREF"] = alignment.aia_ref_line
    aligned.meta["AIAQA"] = alignment.qa_plot

    path = aligned_fits_path(aligned, line, outdir)
    path.parent.mkdir(parents=True, exist_ok=True)
    if overwrite or not path.exists():
        aligned.save(path, overwrite=overwrite)
    aligned.meta["aia_dx"] = float(alignment.dx_arcsec)
    aligned.meta["aia_dy"] = float(alignment.dy_arcsec)
    aligned.meta["aia_corr"] = float(alignment.corr_max)
    aligned.meta["aia_fits"] = alignment.aia_fits
    aligned.meta["aia_time"] = alignment.aia_time
    aligned.meta["aia_ref_line"] = alignment.aia_ref_line
    aligned.meta["aia_qa_plot"] = alignment.qa_plot
    aligned.meta["aia_aligned_fits"] = str(path)
    return aligned
