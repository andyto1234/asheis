from asheis import asheis
import glob
files = glob.glob('/disk/solar9/st3/data_eis/*.data.h5')

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