import sunkit_image.enhance as enhance
import scipy.ndimage as ndimage
from aiapy.calibrate import normalize_exposure, register, update_pointing
import sunpy
from astropy.visualization import ImageNormalize, quantity_support
import numpy as np
from pathlib import Path
from asheis import load_plotting_routine, load_axes_labels
from matplotlib import pyplot as plt
from datetime import datetime

def directory_setup_aia(amap):
    date = datetime.strptime(amap.meta['date-obs'],"%Y-%m-%dT%H:%M:%S.%f")
    wvln = amap.meta['wavelnth']
    fullpath = f'images/mgn/{date.year}_{date.month:02d}_{date.day:02d}/{wvln}/fits/'
    filename = f'mgn_{wvln}A_{date.year}_{date.month:02d}_{date.day:02d}__{date.hour:02d}_{date.minute:02d}_{date.second:02d}'
    Path(fullpath).mkdir(parents=True, exist_ok=True)
    return fullpath, filename

class ashaia:
    def __init__(self, filename):
        self.filename = filename
        self.map = sunpy.map.Map(filename)

    def plot_map(self, date, amap, colorbar=False, savefig=True, **kwargs):
        load_plotting_routine()
        amap.plot()
        if colorbar==True: plt.colorbar() 
        load_axes_labels()
        # plt.savefig(f'{date}/eis_{m.measurement.lower().replace(" ","_").replace(".","_")}.png')
        if savefig==True: plt.savefig(f'images/{amap.measurement.lower().split()[-1]}/eis_{date}_{amap.measurement.lower().replace(" ","_").replace(".","_")}.png')
        # plt.savefig(f'images/{amap.measurement.lower().split()[-1]}/eis_{date}_{amap.measurement.lower().replace(" ","_").replace(".","_")}.png')

    def aiaprep(self, normalise=None):
        prep_map = register(self.map)            
        if self.map.detector == 'AIA':
            prep_map = update_pointing(self.map)
            prep_map = register(prep_map)
            if normalise is not None:
                prep_map = normalize_exposure(prep_map)
        return prep_map
    
    def euv_super_res(self, bottom_left=None, width=None, height=None, savefile = False):
        passband = str(self.map.wavelength).split(" ")[0].split(".")[0]
        h_param = {
            '171': 0.915,
            '193': 0.94,
            '211': 0.94,
            '131': 0.915,
            '94': 0.94,
            '304': 0.95,
            '335': 0.8,
            '1600': 0.9,
            '1700': 0.9,
            '4500': 0.9,
            '6173': 0.8
        }
        smth = 3 if h_param[passband] != '94' else 5
        if bottom_left != None:
            m_pre_mgn = self.map.submap(bottom_left, width=width, height=height)
        else:
            m_pre_mgn = self.map
        m_mgn_data = enhance.mgn(m_pre_mgn.data.astype(float), h=h_param[passband], sigma=[2.5,5,10,20,40,80])
        m_mgn_data = ndimage.gaussian_filter(m_mgn_data, sigma=(smth), order=0)
        m_mgn = sunpy.map.Map(m_mgn_data, self.map.meta)
        q1 = int(len(m_mgn_data) * 0.005 / 2)
        m_mgn.plot_settings['norm'] = ImageNormalize(vmin=np.sort(m_mgn_data.ravel())[int(q1)],vmax=np.sort(m_mgn_data.ravel())[int(q1)]+1)
        fullpath, filename = directory_setup_aia(m_mgn)
        if savefile != False:
            m_mgn.save(f'{fullpath}{filename}.fits', overwrite=True)

        return m_mgn

