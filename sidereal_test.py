from astropy.coordinates import EarthLocation
from astropy.time import Time
from astropy import units as u

observing_location = EarthLocation(lat=33.697432*u.deg, lon=-117.72314*u.deg)
observing_time = Time(datetime.datetime.utcnow(), scale='utc', location=observing_location)
LST = observing_time.sidereal_time('mean')

print(f"Local Sidereal Time = {LST}")