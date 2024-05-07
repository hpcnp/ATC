import ephem
from astroquery.simbad import Simbad
from datetime import date, datetime, timedelta

result_table = Simbad.query_object("m1")
result_table.pprint()

lat = 33.697472      # Latitude    
long = -117.72314     # Longitude
elevation = 22.0           # Elevation (meters)

longitude = ephem.degrees(str(long))
latitude = ephem.degrees(str(lat))
today_date = datetime.now()

obj = ephem.FixedBody()
obj._ra = '01:52:50.4'
obj._dec = '36:08:46'

observer = ephem.Observer()
observer.date = today_date
observer.lon = longitude
observer.lat = latitude

obj.compute(observer)

print(f"Observer: {observer}")
print(f"RA: {obj.ra}, DEC: {obj.dec}")
print(f"AZ: {obj.az}, ALT: {obj.alt}")