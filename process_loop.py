from asheis import asheis
import glob


if __name__ == '__main__':
    # [0,1000] : '/disk/solar9/st3/data_eis/*.data.h5'
    # [1000,2000] : '/disk/solar8/st3/data_eis/*.data.h5'
    files = glob.glob('/disk/solar8/st3/data_eis/*.data.h5')

    for file in files:
        myfits = asheis(file, ncpu="40")
        
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