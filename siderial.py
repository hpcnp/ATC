import numpy as np
from astropy.time import Time

def utc_to_lst(utc_time, longitude):
    """
    Convert UTC time to Local Sidereal Time (LST).
    
    Parameters:
    utc_time (astropy.time.Time): UTC time.
    longitude (float): Observer's longitude in degrees (positive for East, negative for West).
    
    Returns:
    lst (float): Local Sidereal Time in hours.
    """
    # Convert UTC time to Julian Day
    jd = utc_time.jd
    
    # Calculate Greenwich Mean Sidereal Time (GMST)
    t = (jd - 2451545.0) / 36525.0
    gmst = 280.46061837 + 360.98564736629 * t + 0.000387933 * t**2 - t**3/38710000
    gmst = gmst % 360
    
    # Convert GMST to Local Sidereal Time (LST)
    lst = (gmst + longitude) % 360
    lst /= 15  # Convert to hours
    
    return lst

def lst_to_utc(lst, longitude):
    """
    Convert Local Sidereal Time (LST) to UTC time.
    
    Parameters:
    lst (float): Local Sidereal Time in hours.
    longitude (float): Observer's longitude in degrees (positive for East, negative for West).
    
    Returns:
    utc_time (astropy.time.Time): UTC time.
    """
    # Convert LST to Greenwich Mean Sidereal Time (GMST)
    gmst = (lst * 15) % 360
    
    # Calculate Julian Day from GMST
    t = (gmst - 280.46061837) / 360.98564736629
    jd = 2451545.0 + t * 36525.0
    
    # Convert Julian Day to UTC time
    utc_time = Time(jd, format='jd', scale='utc')
    
    return utc_time

# Example usage
utc_time = Time("2023-04-20 20:00:00", format="iso", scale="utc")
longitude = -122.4194  # Observer's longitude in degrees (West)

lst = utc_to_lst(utc_time, longitude)
print(f"Local Sidereal Time: {lst:.2f} hours")

new_utc_time = lst_to_utc(lst, longitude)
print(f"UTC time: {new_utc_time.iso}")
