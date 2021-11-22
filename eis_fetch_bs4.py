import cfscrape
import bs4
import datetime, calendar
import eispac

year = 2020
months = [1,2,3,4,5,6,7,8,9,10,11,12]

for month in months:
    num_days = calendar.monthrange(year, month)[1]
    days=[day+1 for day in range(num_days)]
    month = '{:0>2}'.format(month)
    for day in days:
        day = '{:0>2}'.format(day)
        url = f'https://eis.nrl.navy.mil/level1/hdf5/2021/{month}/{day}/'
        scraper = cfscrape.create_scraper()
        html = scraper.get(url).content
        soup = bs4.BeautifulSoup(html, 'html.parser')
        iframes = soup.find_all(text=True)
        for data in iframes:
            a = eispac.download.download_hdf5_data(data)