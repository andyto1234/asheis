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
