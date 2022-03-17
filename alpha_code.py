import eispac
import numpy as np
from numpy import array
from pathlib import Path
from scipy import integrate
from scipy.interpolate import interp1d
import multiprocessing as mp
from functools import partial
from os import cpu_count
from tqdm import tqdm

def alpha(fit_res, data_cube, iy, ix):
    # Get fit info for pixel.
    ex_coords = [iy, ix] # [Y,X] array coords in units of [pixels]
    fit_x, fit_y = fit_res.get_fit_profile(coords=ex_coords,
                                               num_wavelengths=100) #Total fit; if multiple Guassian then Gaussian 1 + Guassian 2 + Background
    c0_fit_x, c0_fit_y = fit_res.get_fit_profile(component=0,
                                         coords=ex_coords, num_wavelengths=100) #Gaussian 1
    c1_fit_x, c1_fit_y = fit_res.get_fit_profile(component=1,
                                         coords=ex_coords, num_wavelengths=100) #Background

    # Get fit params for pixels.
    cen_vals, cen_errs = fit_res.get_params(component = 0, param_name = 'centroid')
    fwhm_vals, fwhm_errs = fit_res.get_params(component = 0, param_name = 'width')
    peak_vals, peak_errs = fit_res.get_params(component = 0, param_name = 'peak')

    # Get observed intensities for pixel.
    win_wvl = data_cube.wavelength[iy,ix,:]
    obs_spec = data_cube.data[iy,ix,:]

    # Get centroid value and convert velocities in km/s to wavelength.
    vel1w = cen_vals[iy, ix] /((vel1/299792) + 1)
    vel2w = cen_vals[iy, ix] /((vel2/299792) + 1)
    
    # Calc begin/end wavelengths for trimmed curves
    dw = 0.15
    min_w = cen_vals[iy, ix] - dw
    max_w = cen_vals[iy, ix] + dw
    
    
    # Case four:  Calc alpha for observed and fitted curves with trimming curves and using calculated Gaussian from fits params and no background.
    
    
    # Calculate Gaussian curve from fit params.
    mk_wave = np.linspace(min_w[0], max_w[0], 500)
    if peak_vals[iy, ix] == 0 or fwhm_vals[iy, ix] == 0:
        return -100, iy, ix
    else:
        gcurve = peak_vals[iy, ix] * np.exp(-((mk_wave - cen_vals[iy, ix])**2)/(2*(fwhm_vals[iy, ix])**2))

        bg = c1_fit_y[0]
        # Interpolate observed intesity
        interp_obs_i = interp1d(win_wvl, obs_spec-bg, kind='linear')
        new_obs_spec = interp_obs_i(mk_wave)


        # Calculate cumulative area for observed and fitted curves.
        obs_curve = integrate.cumtrapz(new_obs_spec, mk_wave)
        gauss_curve = integrate.cumtrapz(gcurve, mk_wave)

        # Interpolate along each cumulative area curves using cubic spline and adjusting length of arrays by removing the last index.
        f_obs5 = interp1d(mk_wave[0:-1], obs_curve, kind ='cubic')
        f_gauss5 = interp1d(mk_wave[0:-1], gauss_curve, kind ='cubic')

        # Calculate area in velocity interval for each curve.
        obs_area5 = f_obs5(vel1w) - f_obs5(vel2w)
        gauss_area5 = f_gauss5(vel1w) - f_gauss5(vel2w)

        # Calculate alpha
        alpha5 = (obs_area5 - gauss_area5)/((obs_area5 + gauss_area5)/2)
        return alpha5[0], iy, ix

    

def alpha_map(fit_res, data_cube, ncpu): #make it a function of the emission line
    ncpu = cpu_count() if ncpu == 'max' else ncpu
    x=len(data_cube.data[0, :, 0])
    y=len(data_cube.data[:, 0, 0])
    a_map=np.empty((y, x))

    pool = mp.Pool(mp.cpu_count())
    function_base = partial(alpha, fit_res, data_cube)
    for yy in tqdm(range(0, y-1), position=0, leave=True, desc=f"Total Progress - {data_cube.meta['date_obs'][-1]}"):
        function = partial(function_base, yy)
        ix = range(0, x-1)
        results = pool.map(function, ix)
        for result in results:
            a_map[result[1],result[2]] = result[0]
    pool.close()
    pool.join() 

    return a_map

# Set velocity range in bluewing (km/s)
vel1 = 70
vel2 = 170
if __name__ == "__main__":
    alpha()
    alpha_map()
# Get observed (NRL h5 files) info for entire raster.
# data_filename = 'Tests_folder/eis_20071207_001826.data.h5' 
# header_filename = 'Tests_folder/eis_20071207_001826.head.h5'
# tmplt_filename = 'fe_13_202_044.1c.template.h5'
# tmplt_dir = eispac.data.get_fit_template_filepath(tmplt_filename)
# tmplt = eispac.read_template(tmplt_dir)
# data_cube = eispac.read_cube(data_filename, tmplt.central_wave)

# # Get fit info for entire raster (fitted data cube for wavelength)
# fit_res = eispac.read_fit('Tests_folder/eis_20071207_001826.fe_13_202_044.1c-0.fit.h5')

# a_map = alpha_map() #this will be the 2D alpha array

# np.save('alpha_array', a_map)
