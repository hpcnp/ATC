import numpy as np
import astropy.units as u
from astropy.coordinates import SkyCoord, AltAz, EarthLocation, SiderealTime

def radec_to_altaz(ra, dec, latitude, longitude, time):
    """
    Convert right ascension and declination to altitude and azimuth.
    
    Parameters:
    ra (float): Right ascension in degrees.
    dec (float): Declination in degrees.
    latitude (float): Observer's latitude in degrees.
    longitude (float): Observer's longitude in degrees.
    time (astropy.time.Time): Observation time.
    
    Returns:
    alt (float): Altitude in degrees.
    az (float): Azimuth in degrees.
    """
    # Convert input to Astropy coordinates
    coord = SkyCoord(ra=ra*u.degree, dec=dec*u.degree, frame='icrs')
    
    # Get observer location and local sidereal time
    observer = EarthLocation(lat=latitude*u.degree, lon=longitude*u.degree)
    lst = SiderealTime.from_time(time, observer.longitude).value * u.hour
    
    # Convert to altaz frame
    altaz_frame = AltAz(location=observer, obstime=time)
    altaz = coord.transform_to(altaz_frame)
    
    return altaz.alt.value, altaz.az.value

def altaz_to_radec(alt, az, latitude, longitude, time):
    """
    Convert altitude and azimuth to right ascension and declination.
    
    Parameters:
    alt (float): Altitude in degrees.
    az (float): Azimuth in degrees.
    latitude (float): Observer's latitude in degrees.
    longitude (float): Observer's longitude in degrees.
    time (astropy.time.Time): Observation time.
    
    Returns:
    ra (float): Right ascension in degrees.
    dec (float): Declination in degrees.
    """
    # Convert input to Astropy coordinates
    altaz_frame = AltAz(alt=alt*u.degree, az=az*u.degree, location=EarthLocation(lat=latitude*u.degree, lon=longitude*u.degree), obstime=time)
    
    # Convert to ICRS frame
    coord = altaz_frame.transform_to('icrs')
    
    return coord.ra.value, coord.dec.value
