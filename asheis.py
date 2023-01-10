# a lot of information are gathered through Will Barnes Hinode-14 sunpy tutorial - It will be a good manor to cite the sunpy team

from pathlib import Path
import eispac
import sunpy
from matplotlib import colors
import matplotlib.pyplot as plt
from datetime import datetime
from astropy.visualization import ImageNormalize, quantity_support
from alpha_code import alpha, alpha_map


def load_plotting_routine():
    fig = plt.figure()
    fig.set_dpi(300)
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
    def __init__(self, filename, ncpu='max'):
        self.filename = filename
        self.dict = {
            "fe_12_195" : ["fe_12_195_119.2c.template.h5",0],
            "ar_14_194" : ["ar_14_194_396.6c.template.h5",5],
            "ca_14_193" : ["ar_14_194_396.6c.template.h5",1],
            "si_10_258" : ["si_10_258_375.1c.template.h5",0],
            "s_10_264" : ["s__10_264_233.1c.template.h5",0],
            "fe_13_202" : ["fe_13_202_044.1c.template.h5",0],
            "fe_13_203" : ["fe_13_203_826.2c.template.h5",1]
        }
        self.ncpu = ncpu
        
    def fit_data(self,line,product):
        from eispac.instr import ccd_offset
        template_name=self.dict[f'{line}'][0]
        # print(self.filename.replace("data.h5",template_name))
        path = Path(f'{self.filename}'.replace("data.h5",template_name).replace(".template",f"-{self.dict[f'{line}'][1]}.fit"))
        if path.is_file() == False:
            template = eispac.read_template(eispac.data.get_fit_template_filepath(template_name))
            cube = eispac.read_cube(self.filename, window=template.central_wave)
            fit_res = eispac.fit_spectra(cube, template, ncpu=self.ncpu)
            fit_res.fit[f'{product}'] = fit_res.shift2wave(fit_res.fit[f'{product}'],wave=195.119)
            disp = ccd_offset(195.119) - ccd_offset(float(fit_res.fit['line_ids'][self.dict[f'{line}'][1]].split(' ')[-1]))
            fit_res.meta['mod_index']['crval2']-=disp
            save_filepaths = eispac.save_fit(fit_res)
        else:
            fit_res=eispac.read_fit(path)

        return fit_res
    
    def directory_setup(self, amap):
        date = amap.date.strftime("%Y_%m_%d__%H_%M_%S")
        Path(f'images/fits/').mkdir(parents=True, exist_ok=True)
        Path(f'images/{amap.measurement.lower().split()[-1]}/').mkdir(parents=True, exist_ok=True)
        amap.save(f"images/fits/eis_{date}_{'_'.join(amap.measurement.lower().split())}.fits", overwrite=True)
        return date
    
    def plot_map(self, date, amap, colorbar=False, savefig=True, **kwargs):
        load_plotting_routine()
        amap.plot(**kwargs)
        if colorbar==True: plt.colorbar() 
        load_axes_labels()
        # plt.savefig(f'{date}/eis_{m.measurement.lower().replace(" ","_").replace(".","_")}.png')
        if savefig==True: plt.savefig(f'images/{amap.measurement.lower().split()[-1]}/eis_{date}_{amap.measurement.lower().replace(" ","_").replace(".","_")}.png')
        # plt.savefig(f'images/{amap.measurement.lower().split()[-1]}/eis_{date}_{amap.measurement.lower().replace(" ","_").replace(".","_")}.png')

    def get_intensity(self, line):
        fit_res = self.fit_data(line,'int') # Get fitdata
        m = fit_res.get_map(self.dict[f'{line}'][1],measurement='intensity') # From fitdata get map
        date = self.directory_setup(m) # Creating directories
        self.plot_map(date, m) # Plot maps
        return m
        
    def get_velocity(self, line, vmin=-10,vmax=10):
        fit_res = self.fit_data(line,'vel')
        m = fit_res.get_map(component = self.dict[f'{line}'][1],measurement='velocity')
        date = self.directory_setup(m)
        m.plot_settings['norm'] = ImageNormalize(vmin=vmin,vmax=vmax) # adjusting the velocity saturation
        self.plot_map(date, m, colorbar=True)
        
    def get_width(self, line):
        fit_res = self.fit_data(line,'vel')
        m = fit_res.get_map(component = self.dict[f'{line}'][1],measurement='width')
        date = self.directory_setup(m)
        self.plot_map(date, m, colorbar=True)

    def get_alpha(self, line,vmin=0,vmax=0.3,cmap='viridis'):
        template_name=self.dict[f'{line}'][0]
        template = eispac.read_template(eispac.data.get_fit_template_filepath(template_name))
        fit_res = self.fit_data(line,'int')
        data_cube = eispac.read_cube(self.filename, window=template.central_wave)
        m = fit_res.get_map(self.dict[f'{line}'][1],measurement='intensity')
        m.meta['measrmnt'] = 'Alpha'
        # mapvalue = alpha_map(fit_res, data_cube, self.ncpu)
        m = sunpy.map.Map(alpha_map(fit_res, data_cube, self.ncpu), m.meta)
        m.plot_settings['norm'] = ImageNormalize(vmin=vmin,vmax=vmax)
        m.plot_settings['cmap']=plt.get_cmap(cmap)
        date = self.directory_setup(m)
        self.plot_map(date, m, colorbar=True)
        
    def get_composition(self, linepair, vmin=0,vmax=4):
        '''
        This quick look composition code is incomplete and probably doesn't work. Especially be careful of shift2wave code.
        '''
        if linepair == "SiS":
            lines=['si_10_258','s_10_264','Si X-S X 1']
        elif linepair == "CaAr":
            lines=['ca_14_193','ar_14_194','Ca XIV-Ar XIV 1'] 
        else:
            print('No line database can be found. Add your line in code.')
        
        m_nom = self.get_intensity(lines[0])
        m_denom = self.get_intensity(lines[1])
        m_fip = m_nom 
        m_fip.data = m_nom.data/m_denom.data
        m_fip.meta['line_id'] = lines[2]
        m_fip.meta['measrmnt'] = 'composition'
        date = self.directory_setup(m_fip) 
        m_fip.save(f"images/fits/eis_{'_'.join(m_SiS.measurement.lower().split())}.fits", overwrite=True)
        self.plot_map(date, m_fip,vmin=vmin, vmax=vmax, norm=colors.Normalize(), cmap='CMRmap')
        
if __name__ == '__main__':
    load_plotting_routine()
    load_axes_labels()