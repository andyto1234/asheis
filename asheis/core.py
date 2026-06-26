# a lot of information are gathered through Will Barnes Hinode-14 sunpy tutorial - It will be a good manor to cite the sunpy team
"""
A Python package for analyzing Hinode/EIS data with a focus on non-thermal velocity calculations and data alignment.

This package provides tools for:
- Loading and processing Hinode/EIS FITS files
- Calculating non-thermal velocities from spectral line widths
- Automatically aligning data to the Fe XII 195.119Å window
- Handling multiple iron ion species from Fe VIII to Fe XXIV
- Applying various calibration routines (2014, 2023, and 2026 calibrations)

The package builds upon the eispac package and sunpy framework, extending their functionality
for specialized EIS data analysis. Many of the techniques implemented here are based on the
sunpy team's work, particularly from Will Barnes' Hinode-14 sunpy tutorial.

References:
    - Will Barnes' Hinode-14 sunpy tutorial
    - The sunpy team (https://sunpy.org)
    - Hinode/EIS documentation

Example:
    >>> from asheis import asheis
    >>> eis_data = asheis('path/to/eis/file.fits')
"""

from pathlib import Path
import re
import eispac
import sunpy
from matplotlib import colors
import matplotlib.pyplot as plt
from datetime import datetime
import numpy as np
# from alpha_code import alpha, alpha_map
import platform
from astropy.visualization import ImageNormalize, quantity_support
from eis_calibration.eis_cube_calib import calibrate_cube
from asheis.eis_width.calculation import _calculate_non_thermal_velocity_map
from asheis.eis_density.density_config import DENSITY_DIAGNOSTICS
import os
import astropy.units as u

def load_plotting_routine():
    fig = plt.figure()
    fig.set_dpi(150)
    BIGGER_SIZE=9
    SMALLER_SIZE=8
    plt.rc('font', size=BIGGER_SIZE)  # controls default text sizes
    plt.rc('axes', titlesize=BIGGER_SIZE)  # fontsize of the axes title
    plt.rc('axes', labelsize=SMALLER_SIZE)  # fontsize of the x and y labels
    plt.rc('xtick', labelsize=SMALLER_SIZE)  # fontsize of the tick labels
    plt.rc('ytick', labelsize=SMALLER_SIZE)  # fontsize of the tick labels
    plt.rc('legend', fontsize=BIGGER_SIZE)  # legend fontsize
    plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

def load_axes_labels():
    plt.xlabel('x (arcsec)')
    plt.ylabel('y (arcsec)')

class asheis:
    '''
    The fitting routine auto align array and map to the 195.119 window 
    '''
    def __init__(self, filename, ncpu='max', rebin=False):
        self.filename = filename
        """
        The dictionary is in the format of:
        {
            "line_name" : ["template_name", "centroid location", "peak logT"]
        }
        """

        self.dict = {
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
            "mg_5_276.58" : ["mg_05_276_579.1c.template.h5",0, 5.5], # wrong temperature
            "mg_6_268.99" : ["mg_06_268_986.1c.template.h5",0, 5.5], # wrong temperature
            "mg_6_270.40": ["mg_06_270_394.2c.template.h5",0, 5.5], # wrong temperature
            "mg_7_276.15" : ["mg_07_276_153.1c.template.h5",0, 5.5], # wrong temperature
            "mg_7_278.39" : ["my_Mg_VII_template-2c.h5", 0, 5.8],  # Adjust component (0 or 1) and peak logT
            "mg_7_280.74" : ["mg_07_280_737.1c.template.h5",0, 5.8], # wrong temperature
            "o_4_279.63" : ["o__04_279_631.1c.template.h5",0, 5.5], # wrong temperature
            "o_4_279.93" : ["o__04_279_933.1c.template.h5",0, 5.5], # wrong temperature
            "o_5_192.91": ["o__05_192_906.2c.template.h5",1, 5.5], # wrong temperature
            "o_5_248.46": ["o__05_248_456.1c.template.h5",0, 5.5], # wrong temperature
            "o_6_184.12" : ["o__06_184_117.1c.template.h5",0, 5.5], # wrong temperature
            "ca_14_183.46" :["ca_14_183_457.1c.template.h5",0, 6.6],
            "ca_14_193.87" :["ca_14_193_874.6c.template.h5",1, 6.6],
            "ca_15_181.90" :["ca_15_181_900.1c.template.h5",0, 6.6],
            "ca_15_200.97" :["ca_15_200_972.2c.template.h5",0, 6.6],
            "ar_11_188.81" :["ar_11_188_806.3c.template.h5",2, 6.3],
            "ar_14_187.96" : ["ar_14_187_964.1c.template.h5",0, 6.5],
            "ar_14_194.40" : ["ar_14_194_396.6c.template.h5",5, 6.5],
            "ar_14_191.40" : ["ar_14_191_404.2c.template.h5",1, 6.5],
            "si_10_258.37" :["si_10_258_375.1c.template.h5",0, 6.1],
            "s_10_264.23" : ["s__10_264_233.1c.template.h5",0, 6.2],
            "s_11_188.68" : ["s__11_188_675.3c.template.h5",1, 6.3],
            "s_13_256.69" : ["s__13_256_686.1c.template.h5",0, 6.4],
            "he_2_256.32" : ["he_02_256_317.2c.template.h5",0, 4.7],
        }
        self.ncpu = ncpu
        self.rebin = rebin
        self._aia_alignment_cache = {}
        if platform.system() == "Linux":
            self.dens_dir = '/disk/solar17/st3/density'
        elif platform.system() == "Darwin":
            self.dens_dir = '/Users/andysh.to/Script/idl_code/density'


    def check_window(self, line):
        print(f'Checking {line}')
        template_name=self.dict[f'{line}'][0]
        template = eispac.read_template(eispac.data.get_fit_template_filepath(template_name))
        cube = eispac.read_cube(self.filename, window=template.central_wave)
        return cube

    def read_template(self, template_name):
        template = eispac.read_template(eispac.data.get_fit_template_filepath(template_name))
        return template

    def _calibration_name(self, calib):
        if calib is False or calib is None:
            return None
        if calib is True:
            calib = 2023

        calibration_name = str(calib)
        if calibration_name not in {"2014", "2023", "2026"}:
            raise ValueError("calib must be 2014, 2023, 2026, True, False, or None.")

        return calibration_name

    def _fit_cache_dir(self, calib):
        calibration_name = self._calibration_name(calib)
        if calibration_name is None:
            cache_label = "preflight"
        else:
            cache_label = f"calib_{calibration_name}"

        return Path(self.filename).expanduser().parent / f"fit_{cache_label}"

    def _fit_cache_path(self, template_name, line, calib):
        cache_dir = self._fit_cache_dir(calib)
        fit_name = Path(
            f'{self.filename}'.replace("data.h5",template_name).replace(".template",f"-{self.dict[f'{line}'][1]}.fit")
        ).name
        return cache_dir / fit_name
    
    def fit_data(self,line,product,refit, outdir, calib=2023):
        from eispac.instr import ccd_offset

        def ccd_offset_pixel(wavelength):
            wave = u.Quantity(wavelength, u.AA).to(u.AA)
            try:
                offset = ccd_offset(wave)
            except u.UnitConversionError:
                offset = ccd_offset(wave.to_value(u.AA))
            if hasattr(offset, "to_value"):
                values = np.asarray(offset.to_value(u.pixel), dtype=float)
            else:
                values = np.asarray(offset, dtype=float)
            if values.size != 1:
                raise ValueError(f"Expected one CCD offset value for {wavelength}, got shape {values.shape}")
            return float(values.reshape(-1)[0])

        template_name=self.dict[f'{line}'][0]
        # print(self.filename.replace("data.h5",template_name))

        if template_name == 'my_Mg_VII_template-2c.h5':
            template = eispac.read_template(Path(__file__).parent / 'custom_templates/my_Mg_VII_template-2c.h5')
            template_name = 'mg_07_278.39.2c.template.h5'
        elif template_name != 'fe_13_203_826.2c.template.h5':
            template = eispac.read_template(eispac.data.get_fit_template_filepath(template_name))
        else:
            template = eispac.read_template(Path(__file__).parent / 'eis_density/fe_13_203_830.3c.template.h5')
            template_name = 'fe_13_203_830.3c.template.h5'

        calibration_name = self._calibration_name(calib)
        path = self._fit_cache_path(template_name, line, calib)
        # if self rebin != False:
        print(path)
        if path.is_file() == False or refit==True or self.rebin!=False:
            if calibration_name is None:
                cube = eispac.read_cube(self.filename, window=template.central_wave, apply_radcal=True)
            else:
                cube = eispac.read_cube(self.filename, window=template.central_wave, apply_radcal=False)
                cube = calibrate_cube(cube, calibration_name)
                print(f'---------------------Calibrated cube using {cube.meta["calib_source"]}---------------------')
            if self.rebin != False:
                print('Rebinning')
                cube = cube.smooth_cube(self.rebin)
            fit_res = eispac.fit_spectra(cube, template, ncpu=self.ncpu)
            disp = ccd_offset_pixel(195.119*u.AA) - ccd_offset_pixel(fit_res.fit['wave_range'].mean()*u.AA)
            fit_res.meta['mod_index']['crval2'] = float(fit_res.meta['mod_index']['crval2'] - disp)
            path.parent.mkdir(parents=True, exist_ok=True)
            save_filepaths = eispac.save_fit(fit_res, save_dir=path.parent)
        else:
            fit_res=eispac.read_fit(path)
        # shift the width and centroid to the 195.119 window
        if product == 'width':
            fit_res.fit['params'][:,:,2+3*self.dict[f'{line}'][1]] = fit_res.shift2wave(fit_res.fit['params'][:,:,2+3*self.dict[f'{line}'][1]],wave=195.119)
            fit_res.fit['perror'][:,:,2+3*self.dict[f'{line}'][1]] = fit_res.shift2wave(fit_res.fit['perror'][:,:,2+3*self.dict[f'{line}'][1]],wave=195.119)
            fit_res.fit['params'][:,:,2+3*self.dict[f'{line}'][1]][fit_res.fit['perror'][:,:,2+3*self.dict[f'{line}'][1]] > fit_res.fit['params'][:,:,2+3*self.dict[f'{line}'][1]]] = np.nan  # filter out width error > width
        if product in fit_res.fit:
            fit_res.fit[f'{product}'] = fit_res.shift2wave(fit_res.fit[f'{product}'],wave=195.119)
        fit_res.fit['err_int'] = fit_res.shift2wave(fit_res.fit['err_int'],wave=195.119)

        return fit_res
    
    
    def directory_setup(self, amap, line, outdir):
        # Set up directory and save fit
        date = amap.date.strftime("%Y_%m_%d__%H_%M_%S")
        if outdir == '':
            outdir = os.getcwd()

        Path(f'{outdir}/images/fits/').mkdir(parents=True, exist_ok=True)
        Path(f'{outdir}/images/{amap.measurement.lower().split()[-1]}/{line}/').mkdir(parents=True, exist_ok=True)
        # print(date)
        # print(amap.measurement.lower().split())
        # print(f'{date}')
        amap.save(f"{outdir}/images/{amap.measurement.lower().split()[-1]}/{line}/eis_{date}_{'_'.join(amap.measurement.lower().split())}.fits", overwrite=True)
        return date
    
    def plot_map(self, date, amap, line, outdir= os.getcwd(), colorbar=False, savefig=True, **kwargs):
        load_plotting_routine()
        amap.plot(**kwargs)
        if colorbar==True: plt.colorbar() 
        load_axes_labels()
        # plt.savefig(f'{date}/eis_{m.measurement.lower().replace(" ","_").replace(".","_")}.png')
        if savefig==True: plt.savefig(f'{outdir}/images/{amap.measurement.lower().split()[-1]}/{line}/eis_{date}_{amap.measurement.lower().replace(" ","_").replace(".","_")}.png')
        # plt.savefig(f'images/{amap.measurement.lower().split()[-1]}/eis_{date}_{amap.measurement.lower().replace(" ","_").replace(".","_")}.png')

    def _aia_alignment_cache_key(
        self,
        outdir,
        aia_fits,
        aia_cache_dir,
        jsoc_email,
        max_aia_offset_sec,
        qa_plot,
    ):
        return (
            str(Path(self.filename).expanduser()),
            str(Path(outdir).expanduser()),
            str(Path(aia_fits).expanduser()) if aia_fits is not None else "",
            str(Path(aia_cache_dir).expanduser()) if aia_cache_dir is not None else "",
            jsoc_email or os.environ.get("JSOC_EMAIL", ""),
            float(max_aia_offset_sec),
            bool(qa_plot),
        )

    def _get_aia_alignment(
        self,
        target_map,
        line,
        outdir,
        refit,
        calib,
        aia_fits,
        aia_cache_dir,
        jsoc_email,
        max_aia_offset_sec,
        overwrite_alignment,
        qa_plot,
    ):
        from asheis.aia_alignment import compute_aia_alignment

        key = self._aia_alignment_cache_key(
            outdir,
            aia_fits,
            aia_cache_dir,
            jsoc_email,
            max_aia_offset_sec,
            qa_plot,
        )
        if not overwrite_alignment and key in self._aia_alignment_cache:
            return self._aia_alignment_cache[key]

        if line == "fe_12_195.12" and target_map.measurement.lower().split()[-1] == "intensity":
            reference_map = target_map
        else:
            reference_map = self.get_intensity(
                "fe_12_195.12",
                outdir=outdir,
                refit=refit,
                plot=False,
                mcmc=False,
                calib=calib,
                align_aia=False,
            )

        alignment = compute_aia_alignment(
            reference_map,
            outdir=outdir,
            aia_fits=aia_fits,
            aia_cache_dir=aia_cache_dir,
            jsoc_email=jsoc_email,
            max_aia_offset_sec=max_aia_offset_sec,
            overwrite=overwrite_alignment,
            qa_plot=qa_plot,
        )
        self._aia_alignment_cache[key] = alignment
        return alignment

    def _apply_aia_alignment_to_map(
        self,
        amap,
        line,
        outdir,
        refit=False,
        calib=2023,
        aia_fits=None,
        aia_cache_dir=None,
        jsoc_email=None,
        max_aia_offset_sec=60.0,
        overwrite_alignment=False,
        qa_plot=True,
    ):
        from asheis.aia_alignment import apply_aia_alignment

        alignment = self._get_aia_alignment(
            amap,
            line,
            outdir,
            refit,
            calib,
            aia_fits,
            aia_cache_dir,
            jsoc_email,
            max_aia_offset_sec,
            overwrite_alignment,
            qa_plot,
        )
        return apply_aia_alignment(
            amap,
            alignment,
            line=line,
            outdir=outdir,
            overwrite=overwrite_alignment,
        )

    def get_intensity(
        self,
        line,
        outdir=os.getcwd(),
        refit=False,
        plot=True,
        mcmc=False,
        calib=2023,
        align_aia=False,
        aia_fits=None,
        aia_cache_dir=None,
        jsoc_email=None,
        max_aia_offset_sec=60.0,
        overwrite_alignment=False,
        qa_plot=True,
    ):
        if align_aia and mcmc:
            raise ValueError("align_aia=True is not supported with mcmc=True because mcmc returns arrays.")
        fit_res = self.fit_data(line,'int',refit, outdir, calib=calib) # Get fitdata
        m = fit_res.get_map(self.dict[f'{line}'][1],measurement='intensity') # From fitdata get map
        date = self.directory_setup(m,line,outdir) # Creating directories
        if align_aia:
            m = self._apply_aia_alignment_to_map(
                m,
                line,
                outdir,
                refit=refit,
                calib=calib,
                aia_fits=aia_fits,
                aia_cache_dir=aia_cache_dir,
                jsoc_email=jsoc_email,
                max_aia_offset_sec=max_aia_offset_sec,
                overwrite_alignment=overwrite_alignment,
                qa_plot=qa_plot,
            )
        if plot == True: self.plot_map(date, m, line, outdir) # Plot maps
        if mcmc:
            m_error = fit_res.fit['err_int'][:,:,self.dict[f'{line}'][1]]
            return m.data, m_error
        else:
            return m
        
    def get_velocity(
        self,
        line,
        outdir=os.getcwd(),
        vmin=-10,
        vmax=10,
        refit=False,
        plot=True,
        calib=2023,
        align_aia=False,
        aia_fits=None,
        aia_cache_dir=None,
        jsoc_email=None,
        max_aia_offset_sec=60.0,
        overwrite_alignment=False,
        qa_plot=True,
    ):
        fit_res = self.fit_data(line,'vel',refit, outdir, calib=calib)
        m = fit_res.get_map(component = self.dict[f'{line}'][1],measurement='velocity')
        date = self.directory_setup(m,line,outdir)
        m.plot_settings['norm'] = ImageNormalize(vmin=vmin,vmax=vmax) # adjusting the velocity saturation
        if align_aia:
            m = self._apply_aia_alignment_to_map(
                m,
                line,
                outdir,
                refit=refit,
                calib=calib,
                aia_fits=aia_fits,
                aia_cache_dir=aia_cache_dir,
                jsoc_email=jsoc_email,
                max_aia_offset_sec=max_aia_offset_sec,
                overwrite_alignment=overwrite_alignment,
                qa_plot=qa_plot,
            )
        if plot == True: self.plot_map(date, m, line, outdir, colorbar=True)
        return m
    
    def get_width(
        self,
        line,
        outdir=os.getcwd(),
        refit=False,
        plot=True,
        width_only=False,
        calib=2023,
        align_aia=False,
        aia_fits=None,
        aia_cache_dir=None,
        jsoc_email=None,
        max_aia_offset_sec=60.0,
        overwrite_alignment=False,
        qa_plot=True,
    ):
        fit_res = self.fit_data(line,'width', refit, outdir, calib=calib)
        cent, cent_error = fit_res.get_params(component = self.dict[f'{line}'][1], param_name = 'centroid')
        m = fit_res.get_map(component = self.dict[f'{line}'][1],measurement='width')
        m.meta['slitwid'] = [float(x) for x in fit_res.meta['slit_width']]  # Store as list of floats
        m.meta['cent'] = cent.tolist()  # Store as list of floats
        m.meta['centerr'] = cent_error.tolist()  # Store as list of floats
        m.meta['measrmnt'] = 'ntv'
        date = self.directory_setup(m,line,outdir)
            # Process the map based on width_only flag
        if width_only:
            final_map = m
        else:
            final_map = _calculate_non_thermal_velocity_map(m)
            final_map.meta['measrmnt']='NTV'
            final_map.meta['bunit']='km/s'
        if align_aia:
            final_map = self._apply_aia_alignment_to_map(
                final_map,
                line,
                outdir,
                refit=refit,
                calib=calib,
                aia_fits=aia_fits,
                aia_cache_dir=aia_cache_dir,
                jsoc_email=jsoc_email,
                max_aia_offset_sec=max_aia_offset_sec,
                overwrite_alignment=overwrite_alignment,
                qa_plot=qa_plot,
            )
        if plot == True: self.plot_map(date, final_map, line, outdir, colorbar=True)
        return final_map
        
    def get_density(self, diagnostic='fe_13', outdir=os.getcwd(), refit=False, plot=True, mcmc=False, **kwargs):
        from scipy.io import readsav

        """
        Calculate electron density using various EIS line pair diagnostics.
        
        Parameters
        ----------
        diagnostic : str
            The density diagnostic to use. Options:
            - 'fe_13': Fe XIII 203.83/202.04 ratio
            - 'fe_12': Fe XII 196.64/195.12 ratio
            (Add more as needed)
        """
        # Dictionary of density diagnostics
        if diagnostic not in DENSITY_DIAGNOSTICS:
            raise ValueError(f"Invalid diagnostic. Choose from: {list(DENSITY_DIAGNOSTICS.keys())}")

        diag_info = DENSITY_DIAGNOSTICS[diagnostic]
        
        # Load density calibration data
        density_file = Path(__file__).parent / f'eis_density/{diag_info["sav_file"]}'
        sav_data = readsav(str(density_file))
        density_ratios = sav_data['smooth_rat']
        density_values = sav_data['smooth_den']

        # Get intensity maps for line pair
        m_nom = self.get_intensity(diag_info['nom_line'], outdir, plot=False, **kwargs)
        m_denom = self.get_intensity(diag_info['denom_line'], outdir, plot=False, **kwargs)
        obs_ratio = m_nom.data / m_denom.data

        # Calculate density using lookup table
        density_map = np.zeros_like(obs_ratio)
        for i in range(obs_ratio.shape[0]):
            for j in range(obs_ratio.shape[1]):
                obs = obs_ratio[i, j]
                closest_index = np.argmin(np.abs(density_ratios - obs))
                density_map[i, j] = density_values[closest_index]

        # Create density map
        m = sunpy.map.Map(density_map, m_nom.meta)
        m.meta['measrmnt'] = 'density'
        m.meta['bunit'] = 'cm^-3'  # Set proper density units
        m.meta['line_id'] = diag_info['description']
        m.plot_settings['norm'] = ImageNormalize(vmin=diag_info['vmin'], vmax=diag_info['vmax'])

        if plot:
            date = self.directory_setup(m, diagnostic, outdir)
            self.plot_map(date, m, diagnostic, outdir, colorbar=True)

        return m.data if mcmc else m



    def get_composition(self, linepair, outdir = os.getcwd(), vmin=0, vmax=4, **kwargs):
        '''
        This quick look composition code is incomplete and probably doesn't work. Especially be careful of shift2wave code.
        '''
        line_databases = {
            "SiS": ['si_10_258', 's_10_264', 'Si X-S X'],
            "CaAr": ['ca_14_193', 'ar_14_194', 'Ca XIV-Ar XIV'],
            "FeS": ['fe_16_262', 's_13_256', 'Fe XVI-S XIII'],
            "SAr": ['s_13_256', 'ar_14_194', 'S XIII-Ar XIV'],
            "SAr2": ['s_11_188', 'ar_11_188', 'S XI-Ar XI'],
            # "FeAr": ['fe_16_262', 'ar_11_188', 'Fe XVI-Ar XI'], wrong dignostic
            # Add more line pairs as needed
        }

        if linepair not in line_databases:
            print('No line database can be found. Add your line in code.')
            return

        lines = line_databases[linepair]

        m_nom = self.get_intensity(lines[0], outdir)
        m_denom = self.get_intensity(lines[1], outdir)
        m_ref = self.get_intensity('fe_12_195', outdir, plot=False)

        m_fip = sunpy.map.Map(m_nom.data / m_denom.data, m_ref.meta)
        m_fip.meta['line_id'] = lines[2]
        m_fip.meta['measrmnt'] = 'composition'
        m_fip.meta.pop('bunit', None)

        date = self.directory_setup(m_fip, lines[2])
        
        # m_fip.save(f"images/{m_fip.measurement.lower().split()[-1]}/{lines[2]}/eis_{date}_{'_'.join(m_fip.measurement.lower().split())}.fits", overwrite=True)
        self.plot_map(date, m_fip, lines[2], outdir, vmin=vmin, vmax=vmax, norm=colors.Normalize(), cmap='CMRmap')
        
if __name__ == '__main__':
    load_plotting_routine()
    load_axes_labels()
    # def get_alpha(self, line,vmin=0,vmax=0.3,cmap='viridis'):
    #     template_name=self.dict[f'{line}'][0]
    #     template = eispac.read_template(eispac.data.get_fit_template_filepath(template_name))
    #     fit_res = self.fit_data(line,'int')
    #     data_cube = eispac.read_cube(self.filename, window=template.central_wave)
    #     m = fit_res.get_map(self.dict[f'{line}'][1],measurement='intensity')
    #     m.meta['measrmnt'] = 'Alpha'
    #     # mapvalue = alpha_map(fit_res, data_cube, self.ncpu)
    #     m = sunpy.map.Map(alpha_map(fit_res, data_cube, self.ncpu), m.meta)
    #     m.plot_settings['norm'] = ImageNormalize(vmin=vmin,vmax=vmax)
    #     m.plot_settings['cmap']=plt.get_cmap(cmap)
    #     date = self.directory_setup(m)
    #     self.plot_map(date, m, colorbar=True)
