import cfscrape
import bs4
import datetime, calendar
import eispac
import multiprocessing as mp
from functools import partial


years = [2017,2018,2019,2020]
months = [1,2,3,4,5,6]


function = partial(eispac.download.download_hdf5_data, local_top='/disk/solar9/st3/data_eis')

places = []
file=open('/disk/solar17/st3/python/asheis/data.txt','r')
for line in file:
    currentPlace = line[:-1]
    places.append(currentPlace)

# for place in places:
pool = mp.Pool(processes=7)
mylist = [str(x) for y,x in enumerate(places)]
pool.map(function, mylist)

# for year in years:
#     for month in months:
#         num_days = calendar.monthrange(year, month)[1]
#         days=[day+1 for day in range(num_days)]
#         month = '{:0>2}'.format(month)
#         for day in days:
#             day = '{:0>2}'.format(day)
#             url = f'https://eis.nrl.navy.mil/level1/hdf5/{year}/{month}/{day}/'
#             scraper = cfscrape.create_scraper()
#             html = scraper.get(url).content
#             soup = bs4.BeautifulSoup(html, 'html.parser')
#             iframes = soup.find_all(text=True)
#             iframes = [i for i in iframes if 'data.h5' in i]
#             pool = mp.Pool(processes=5)
#             mylist = [str(x) for y,x in enumerate(iframes)]
#             pool.map(function, mylist)
