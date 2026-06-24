# asheis

A Python package for calibrating and analyzing EIS (EUV Imaging Spectrometer) data.

## Installation
```bash
git clone https://github.com/andyto1234/asheis.git
cd asheis
pip install -e .
```

## Usage

```python
from asheis.core import asheis
test = asheis('path/to/your/.data.h5')
fe12int = test.get_intensity('fe_12_195.12') # intensity
fe12doppler = test.get_velocity('fe_12_195.12') # doppler velocity
fe12ntv = test.get_width('fe_12_195.12') # non-thermal velocity
```

## Calibration

By default, `asheis` now calibrates the EIS spectral cube before fitting:

```python
fe12int_2014 = test.get_intensity("fe_12_195.12", calib=2014)
fe12int_2023 = test.get_intensity("fe_12_195.12", calib=2023)
fe12int_2026 = test.get_intensity("fe_12_195.12", calib=2026)
```

Use `calib=False` or `calib=None` to keep the standard EISPAC preflight
calibration. Calibrated fit files are cached separately in directories such as
`fit_calib_2014`, `fit_calib_2023`, `fit_calib_2026`, and `fit_preflight` next
to the input EIS data file.

## AIA alignment

Install the optional AIA dependencies when you want JSOC/AIA alignment:

```bash
pip install -e ".[aia]"
```

Then use the existing getters with `align_aia=True`:

```python
from asheis.core import asheis

test = asheis("path/to/eis.data.h5")

aligned_intensity = test.get_intensity("fe_16_262.98", align_aia=True)
aligned_velocity = test.get_velocity("fe_12_195.12", align_aia=True)
aligned_width = test.get_width("fe_13_202.04", align_aia=True)
```

The alignment shift is always measured from Fe XII 195.12 intensity against
AIA 193, then applied to the header/WCS of the requested line or product. If
`aia_fits` is not supplied, `asheis` checks the local AIA cache first and then
downloads the nearest JSOC AIA 193 FITS file. Set a registered JSOC email with:

```bash
export JSOC_EMAIL="you@example.com"
```

The returned object remains a SunPy map. Alignment metadata is added to the map
header, including `aia_dx`, `aia_dy`, `aia_corr`, `aia_fits`, `aia_time`,
`aia_ref_line`, `aia_qa_plot`, and `aia_aligned_fits`.
