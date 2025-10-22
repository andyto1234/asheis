import numpy as np
import datetime

from tqdm import tqdm
from eis_calibration.eis_calib_2023 import _calib_2023_ratio_only
from eis_calibration.eis_calib_2014 import eis_ea
import pickle

def get_line_wavelengths():
    """
    Extracts the available wavelengths from the line dictionary in asheis/core.py.
    Returns:
        line_names: list of line names (str)
        wavelengths: list of wavelengths (float)
    """
    # Copy-paste the line dictionary here or import it if possible
    line_dict = {
        "fe_8_185.21" : ["fe_08_185_213.1c.template.h5",0, 5.7],
        "fe_8_186.60" : ["fe_08_186_601.1c.template.h5",0, 5.7],
        "fe_9_188.50" : ["fe_09_188_497.3c.template.h5",0, 5.9],
        "fe_9_197.86" : ["fe_09_197_862.1c.template.h5",0, 5.9],
        "fe_10_184.54" : ["fe_10_184_536.1c.template.h5",0, 6.0],
        "fe_11_188.22" : ["fe_11_188_216.2c.template.h5",0, 6.1],
        "fe_11_188.30" : ["fe_11_188_299.2c.template.h5",1, 6.1],
        "fe_12_186.88" : ["fe_12_186_880.1c.template.h5",0, 6.2],
        "fe_12_195.12" : ["fe_12_195_119.2c.template.h5",0, 6.2],
        "fe_12_192.39" : ["fe_12_192_394.1c.template.h5",0, 6.2],
        "fe_12_196.64" : ["fe_12_196_640.2c.template.h5",0, 6.2],
        "fe_13_202.04" : ["fe_13_202_044.1c.template.h5",0, 6.2],
        "fe_13_203.83" : ["fe_13_203_826.2c.template.h5",2, 6.2],
        "fe_14_264.79" : ["fe_14_264_787.1c.template.h5",0, 6.3],
        "fe_14_270.52" : ["fe_14_270_519.2c.template.h5",1, 6.3],
        "fe_15_284.16" : ["fe_15_284_160.2c.template.h5",1, 6.3],
        "fe_16_262.98" : ["fe_16_262_984.1c.template.h5",0, 6.4],
        "fe_17_254.87" : ["fe_17_254_870.2c.template.h5",0, 6.7],
        "fe_22_253.10" : ["fe_22_253_170.1c.template.h5",0, 7.2],
        "fe_23_263.76" : ["fe_23_263_760.1c.template.h5",0, 7.2],
        "fe_24_255.10" : ["fe_24_255_100.2c.template.h5",1, 7.3],
        "ca_14_193.87" :["ca_14_193_874.6c.template.h5",1, 6.6],
        "ca_15_181.90" :["ca_15_181_900.1c.template.h5",0, 6.6],
        "ca_15_200.97" :["ca_15_200_972.2c.template.h5",0, 6.6],
        "ar_11_188.81" :["ar_11_188_806.3c.template.h5",2, 6.3],
        "ar_14_194.40" : ["ar_14_194_396.6c.template.h5",5, 6.5],
        "ar_14_191.40" : ["ar_14_191_404.2c.template.h5",1, 6.5],
        "si_10_258.37" :["si_10_258_375.1c.template.h5",0, 6.1],
        "s_10_264.23" : ["s__10_264_233.1c.template.h5",0, 6.2],
        "s_11_188.68" : ["s__11_188_675.3c.template.h5",1, 6.3],
        "s_13_256.69" : ["s__13_256_686.1c.template.h5",0, 6.4],
        "he_2_256.32" : ["he_02_256_317.2c.template.h5",0, 4.7],
    }
    line_names = []
    wavelengths = []
    for k in line_dict:
        try:
            wvl = float(k.split('_')[-1])
            line_names.append(k)
            wavelengths.append(wvl)
        except Exception:
            continue
    return line_names, np.array(wavelengths)

def daterange(start_date, end_date):
    """
    Generator for dates from start_date to end_date (inclusive).
    """
    for n in range(int((end_date - start_date).days) + 1):
        yield start_date + datetime.timedelta(n)

def generate_calib_factor_grid(start_date_str, end_date_str, save_path='calib_factor_grid.pkl'):
    """
    Precompute the calibration factor for each day and each available wavelength.
    Args:
        start_date_str: str, e.g. '2007-01-01'
        end_date_str: str, e.g. '2025-12-31'
        save_path: str, where to save the result (pickle file)
    Returns:
        dates: list of date strings
        wavelengths: list of wavelengths
        calib_factors: 2D numpy array [ndates, nwavelengths]
    """
    line_names, wavelengths = get_line_wavelengths()
    start_date = datetime.datetime.strptime(start_date_str, "%Y-%m-%d")
    end_date = datetime.datetime.strptime(end_date_str, "%Y-%m-%d")
    dates = []
    for single_date in daterange(start_date, end_date):
        dates.append(single_date.strftime("%Y-%m-%dT00:00:00"))
    ndates = len(dates)
    nwvl = len(wavelengths)
    calib_factors = np.zeros((ndates, nwvl))
    for i, date in tqdm(enumerate(dates), desc="Processing dates", total=ndates):
        for j, wvl in enumerate(wavelengths):
            try:
                calib_factors[i, j] = _calib_2023_ratio_only(date, wvl)
            except Exception as e:
                print(f"Failed for date {date}, wavelength {wvl}: {e}")
                calib_factors[i, j] = np.nan
    # Save to file
    with open(save_path, 'wb') as f:
        pickle.dump({'dates': dates, 'wavelengths': wavelengths, 'calib_factors': calib_factors}, f)
    print(f"Saved calibration factor grid to {save_path}")
    return dates, wavelengths, calib_factors

def load_calib_factor_grid(save_path='calib_factor_grid.pkl'):
    """
    Load the precomputed calibration factor grid.
    Returns:
        dates: list of date strings
        wavelengths: list of wavelengths
        calib_factors: 2D numpy array [ndates, nwavelengths]
    """
    with open(save_path, 'rb') as f:
        data = pickle.load(f)
    return data['dates'], data['wavelengths'], data['calib_factors']

def interp_calib_factor(query_date, query_wavelength, dates, wavelengths, calib_factors):
    """
    Interpolate the calibration factor for a given date and wavelength.
    Args:
        query_date: str, e.g. '2022-05-01T00:00:00'
        query_wavelength: float
        dates: list of date strings
        wavelengths: list of wavelengths
        calib_factors: 2D numpy array [ndates, nwavelengths]
    Returns:
        interpolated calibration factor (float)
    """
    # Convert dates to ordinal for interpolation
    date_ordinals = np.array([datetime.datetime.strptime(d, "%Y-%m-%dT%H:%M:%S").toordinal() for d in dates])
    query_ordinal = datetime.datetime.strptime(query_date, "%Y-%m-%dT%H:%M:%S").toordinal()
    # 2D interpolation: first over date, then over wavelength
    from scipy.interpolate import RegularGridInterpolator
    interp_func = RegularGridInterpolator((date_ordinals, wavelengths), calib_factors, bounds_error=False, fill_value=None)
    return interp_func([[query_ordinal, query_wavelength]])[0]

# Example usage:
# 1. Generate and save the grid (run once)
# dates, wavelengths, calib_factors = generate_calib_factor_grid('2007-04-01', '2025-12-31', save_path='calib_factor_grid.pkl')
#
# 2. Load and interpolate
# dates, wavelengths, calib_factors = load_calib_factor_grid('calib_factor_grid.pkl')
# factor = interp_calib_factor('2022-05-01T00:00:00', 195.12, dates, wavelengths, calib_factors)
# print(f"Interpolated calibration factor: {factor}")
