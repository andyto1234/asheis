import numpy as np
import sunpy.map
from asheis.eis_width.util import get_line_info, eis_element2mass

def get_ntv(width):
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
    line = str(width.meta['line_id'])
    print(line)
    yy, xx = width.data.shape   

    # Get instrumental width array from width map
    inst_wid_arr = np.empty([yy, xx])
    slit_wid = width.meta['slit_width']
    for k in range(0, xx):
        inst_wid_arr[0:, k] = slit_wid
        
    # Get thermal width
    cent = width.meta['cent']
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
    width_fwhm = width.data*2*np.sqrt(2*np.log(2))

    # Subtracting instrumental width and thermal width from the sigma FWHM 
    dl_nt_2 = width_fwhm**2 - np.array(inst_wid_arr)**2 - np.array(thermal_fwhm)**2
    
    # Setting negative values to NaN
    dl_nt_2[dl_nt_2 <= 0] = np.nan
    
    # Converting FWHM to non-thermal velocity
    v_nt = speed_of_light*np.sqrt(dl_nt_2)/th_wvl/np.sqrt(4*np.log(2))
    
    return v_nt

def _calculate_non_thermal_velocity_map(width_map):
    v_nt = get_ntv(width_map)
    return sunpy.map.Map(v_nt, width_map.meta)
