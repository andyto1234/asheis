import eispac
import cfscrape
import bs4

url = 'https://eis.nrl.navy.mil/level1/hdf5/2021/11/03/'
scraper = cfscrape.create_scraper()
html = scraper.get(url).content

soup = bs4.BeautifulSoup(html, 'html.parser')
iframes = soup.find_all(text=True)
iframes = [i for i in iframes if 'data.h5' in i]

for i in iframes:
    a = eispac.download.download_hdf5_data(i)