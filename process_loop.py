from asheis import asheis
import glob
from asheis import asheis
files = glob.glob('/disk/solar17/st3/python/data_eis/*.data.h5')

for file in files:
    myfits = asheis(file)
    
    try:
        myfits.get_intensity('fe_12_195')
    except:
        pass

    try:
        myfits.get_velocity('fe_12_195')
    except:
        pass
    
    try:
        myfits.get_width('fe_12_195')
    except:
        pass
    
    try:
        myfits.get_composition('SiS')
    except:
        pass