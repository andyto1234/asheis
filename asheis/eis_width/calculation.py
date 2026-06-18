import numpy as np
import sunpy.map
from asheis.eis_width.util import get_line_info, eis_element2mass

def _instrument_width_array(slit_width, shape):
    """Return slit width values broadcast to the EIS map shape."""
    yy, xx = shape
    slit_width = np.asarray(slit_width, dtype=float)
    if slit_width.shape == shape:
        return slit_width
    if slit_width.size == 1:
        return np.full(shape, float(slit_width.reshape(-1)[0]))
    if slit_width.ndim == 1 and slit_width.size == yy:
        return np.repeat(slit_width[:, np.newaxis], xx, axis=1)
    if slit_width.ndim == 1 and slit_width.size == xx:
        return np.repeat(slit_width[np.newaxis, :], yy, axis=0)
    try:
        return np.broadcast_to(slit_width, shape).astype(float)
    except ValueError as exc:
        raise ValueError(
            f"Cannot broadcast slit_width with shape {slit_width.shape} to width map shape {shape}"
        ) from exc


def _get_ntv(width_data, line_id, slit_width, cent):
    """Calculate non-thermal velocity from EIS width data. This is heavily based on the 
    wrapper asheis.
    
    Parameters
    ----------
    width : EISMap
        Width map data from EIS
    line : str
        Spectral line (e.g. 'fe_12_195.12')
        
    Returns
    -------
    v_nt : numpy.ndarray
        Non-thermal velocity map in km/s
    """
    speed_of_light = 2.9979e5  # Speed of light in km/s
    line = str(line_id)
    print(line)
    yy, xx = width_data.shape   

    # Get instrumental width array from width map
    inst_wid_arr = _instrument_width_array(slit_width, width_data.shape)
        
    # Get thermal width
    ref_wvl = np.median(cent)
    t_max = get_line_info(line)['T_MAX']
    print(t_max)
    th_wvl = ref_wvl
    th_temp = 10**(t_max)
    mass = eis_element2mass(line)
    print(mass)
    therm_wid = np.sqrt(2 * th_temp * 1.6022E-12 / 11604.0 / mass) / 1.0E5
    
    # Calculate thermal FWHM in Angstroms
    thermal_fwhm = np.sqrt(4*np.log(2))*th_wvl*therm_wid/speed_of_light
    
    # Converting width to FWHM
    width_fwhm = width_data*2*np.sqrt(2*np.log(2))

    # Subtracting instrumental width and thermal width from the sigma FWHM 
    dl_nt_2 = width_fwhm**2 - np.array(inst_wid_arr)**2 - np.array(thermal_fwhm)**2
    
    # Setting negative values to NaN
    dl_nt_2[dl_nt_2 <= 0] = np.nan
    
    # Converting FWHM to non-thermal velocity
    v_nt = speed_of_light*np.sqrt(dl_nt_2)/th_wvl/np.sqrt(4*np.log(2))
    
    return v_nt

def _calculate_non_thermal_velocity_map(width_map):
    v_nt = _get_ntv(width_map.data, width_map.meta['line_id'], width_map.meta['slit_width'], width_map.meta['cent'])
    return sunpy.map.Map(v_nt, width_map.meta)
