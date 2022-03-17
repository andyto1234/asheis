import astropy.units as u
from astropy.time import TimeDelta

from sunpy.net.dataretriever import GenericClient
from sunpy.time import TimeRange
from sunpy.util.scraper import Scraper

__all__ = ['EISPACClient']

class EISPACClient(GenericClient):
    """
    Provides access to Level 1 .h5 Extremeultraviolet Imaging Spectrometer (EIS NRL) data.

    To use this client you must request Level 1 .h5 header and data.
    It is hosted by `NRL <https://eis.nrl.navy.mil/level1/hdf5/>`__.

    Examples
    --------

    >>> from sunpy.net import Fido, attrs as a
    >>> results = Fido.search(a.Time("2016/1/1", "2016/1/2"),
    ...                       a.Instrument('EIS'), a.Provider('NRL'), a.Source('Hinode'))  #doctest: +REMOTE_DATA
    >>> results  #doctest: +REMOTE_DATA
    <sunpy.net.fido_factory.UnifiedResponse object at ...>
    Results from 1 Provider:
    <BLANKLINE>
    2 Results from the EVEClient:
    Source: https://lasp.colorado.edu/eve/data_access/evewebdata/quicklook/L0CS/SpWx/
    <BLANKLINE>
           Start Time               End Time        Instrument ... Provider Level
    ----------------------- ----------------------- ---------- ... -------- -----
    2016-01-01 00:00:00.000 2016-01-01 23:59:59.999        EVE ...     LASP     0
    2016-01-02 00:00:00.000 2016-01-02 23:59:59.999        EVE ...     LASP     0
    <BLANKLINE>
    <BLANKLINE>

    """
    baseurl = (r'https://eis.nrl.navy.mil/level1/hdf5/'
               r'%Y/%m/%d/eis_%Y%m%d_%H%M%S.(\w){4}.h5')
    pattern = '{}/{year:4d}/{month:2d}/{day:2d}/eis_{year:4d}{month:2d}{day:2d}_{hour:2d}{minute:2d}{second:2d}.{}'

    @property
    def info_url(self):
        return 'https://eis.nrl.navy.mil/level1/hdf5/'

    @classmethod
    def register_values(cls):
        from sunpy.net import attrs
        adict = {attrs.Instrument: [('EIS', 'Extreme ultraviolet Variability Experiment, which is part of the NASA Solar Dynamics Observatory mission.')],
                 attrs.Source: [('Hinode', 'The Solar Dynamics Observatory.')],
                 attrs.Provider: [('NRL', 'The Laboratory for Atmospheric and Space Physics.')],
                 attrs.Level: [('1', 'EVE: The specific EVE client can only return Level 0C data. Any other number will use the VSO Client.')]}
        return adict
