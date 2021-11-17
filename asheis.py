# a lot of information are gathered through Will Barnes Hinode-14 sunpy tutorial - It will be a good manor to cite the sunpy team

from pathlib import Path
import eispac
import sunpy
from matplotlib import colors
import matplotlib.pyplot as plt
from datetime import datetime


def plotting_routine():
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

def axes_labels():
    plt.xlabel('x (arcsec)')
    plt.ylabel('y (arcsec)')


class asheis:
    
    def __init__(self, filename):
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

        
    def get_intensity(self, line):
        template_name=self.dict[f'{line}'][0]
        template = eispac.read_template(eispac.data.get_fit_template_filepath(template_name))
        cube = eispac.read_cube(self.filename, window=template.central_wave)
        fit_res = eispac.fit_spectra(cube, template, ncpu='max')
        fit_res.fit['int'] = fit_res.shift2wave(fit_res.fit['int'],wave=195.119)
        m = fit_res.get_map(self.dict[f'{line}'][1],measurement='intensity')
        date = m.date.strftime("%Y_%m_%d__%H_%M_%S")
        # Path(f'{date}/fitted_data/').mkdir(parents=True, exist_ok=True)
        # m.save(f"{date}/fitted_data/eis_{'_'.join(m.measurement.lower().split())}.fits", overwrite=True)
        Path(f'fitted_data/').mkdir(parents=True, exist_ok=True)
        m.save(f"fitted_data/eis_{'_'.join(m.measurement.lower().split())}.fits", overwrite=True)
        plotting_routine()
        m.plot()
        axes_labels()
        plt.savefig(f'{date}/eis_{m.measurement.lower().replace(" ","_").replace(".","_")}.png')
        
    def get_velocity(self, line):
        template_name=self.dict[f'{line}'][0]
        template = eispac.read_template(eispac.data.get_fit_template_filepath(template_name))
        cube = eispac.read_cube(self.filename, window=template.central_wave)
        fit_res = eispac.fit_spectra(cube, template, ncpu='max')
        fit_res.fit['int'] = fit_res.shift2wave(fit_res.fit['vel'],wave=195.119)
        m = fit_res.get_map(component = self.dict[f'{line}'][1],measurement='velocity')
        date = m.date.strftime("%Y_%m_%d__%H_%M_%S")
        # Path(f'{date}/fitted_data/').mkdir(parents=True, exist_ok=True)
        # m.save(f"{date}/fitted_data/eis_{'_'.join(m.measurement.lower().split())}.fits", overwrite=True)
        Path(f'fitted_data/').mkdir(parents=True, exist_ok=True)
        m.save(f"fitted_data/eis_{'_'.join(m.measurement.lower().split())}.fits", overwrite=True)
        plotting_routine()
        m.plot()
        plt.colorbar()
        axes_labels()
        # plt.savefig(f'{date}/eis_{m.measurement.lower().replace(" ","_").replace(".","_")}.png')
        plt.savefig(f'eis_{m.measurement.lower().replace(" ","_").replace(".","_")}.png')
        
    def get_width(self, line):
        template_name=self.dict[f'{line}'][0]
        template = eispac.read_template(eispac.data.get_fit_template_filepath(template_name))
        cube = eispac.read_cube(self.filename, window=template.central_wave)
        fit_res = eispac.fit_spectra(cube, template, ncpu='max')
        fit_res.fit['int'] = fit_res.shift2wave(fit_res.fit['vel'],wave=195.119)
        m = fit_res.get_map(component = self.dict[f'{line}'][1],measurement='width')
        date = m.date.strftime("%Y_%m_%d__%H_%M_%S")
        # Path(f'{date}/fitted_data/').mkdir(parents=True, exist_ok=True)
        # m.save(f"{date}/fitted_data/eis_{'_'.join(m.measurement.lower().split())}.fits", overwrite=True)
        Path(f'fitted_data/').mkdir(parents=True, exist_ok=True)
        m.save(f"fitted_data/eis_{'_'.join(m.measurement.lower().split())}.fits", overwrite=True)
        plotting_routine()
        m.plot()
        plt.colorbar()
        axes_labels()
        # plt.savefig(f'{date}/eis_{m.measurement.lower().replace(" ","_").replace(".","_")}.png')
        plt.savefig(f'eis_{m.measurement.lower().replace(" ","_").replace(".","_")}.png')


    def get_composition(self, linepair, vmin=0,vmax=4):
        if linepair == "SiS":
            lines=['si_10_258','s_10_264','Si X-S X 1']
        elif linepair == "CaAr":
            lines=['ca_14_193','ar_14_194','Ca XIV-Ar XIV 1'] 
        else:
            print('No line database can be found. Add your line in code.')
        
        template_names=[self.dict[lines[0]][0],self.dict[lines[1]][0]]
        templates = [eispac.read_template(eispac.data.get_fit_template_filepath(t)) for t in template_names]
        for t in templates:
            cube = eispac.read_cube(self.filename, window=t.central_wave)
            fit_res = eispac.fit_spectra(cube, t, ncpu='max')
            fit_res.fit['int'] = fit_res.shift2wave(fit_res.fit['int'],wave=195.119)
            m = fit_res.get_map(component = self.dict[lines[0]][1], measurement='intensity')
            date = m.date.strftime("%Y_%m_%d__%H_%M_%S")
            # Path(f'{date}/fitted_data/').mkdir(parents=True, exist_ok=True)
            # m.save(f"{date}/fitted_data/eis_{'_'.join(m.measurement.lower().split())}.fits", overwrite=True)
            Path(f'fitted_data/').mkdir(parents=True, exist_ok=True)
            m.save(f"fitted_data/eis_{'_'.join(m.measurement.lower().split())}.fits", overwrite=True)
            lines.append(f"eis_{'_'.join(m.measurement.lower().split())}.fits")
        m_Si = sunpy.map.Map(f'{date}/fitted_data/{lines[-2]}')
        m_S = sunpy.map.Map(f'{date}/fitted_data/{lines[-1]}')
        m_SiS = m_Si
        m_SiS.meta['line_id'] = lines[2]
        m_SiS.save(f"{date}/fitted_data/eis_{'_'.join(m_SiS.measurement.lower().split())}.fits", overwrite=True)
        m_SiS = sunpy.map.Map(m_Si.data/m_S.data, m_Si.meta)
        # m_SiS.peek(vmax=4,cmap='RdYlBu')
        plotting_routine()
        m_SiS.plot(vmin=vmin, vmax=vmax, norm=colors.Normalize(), cmap='CMRmap')
        # plt.figure()
        plt.colorbar()
        axes_labels()
        # plt.savefig(f'{date}/eis_{m_SiS.measurement.lower().replace(" ","_").replace(".","_")}.png')
        plt.savefig(f'eis_{m.measurement.lower().replace(" ","_").replace(".","_")}.png')
