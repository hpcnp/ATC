import math
from datetime import datetime, timedelta

def calculate_moon_ra_dec(year, month, day, hour, minute, second, latitude, longitude):
    """
    Calculates the right ascension (RA) and declination (Dec) of the moon for a given date, time, and location.
    
    Parameters:
    year (int): The 4-digit year.
    month (int): The month (1-12).
    day (int): The day of the month (1-31).
    hour (int): The hour of the day (0-23).
    minute (int): The minute of the hour (0-59).
    second (int): The second of the minute (0-59).
    latitude (float): The latitude of the location in degrees.
    longitude (float): The longitude of the location in degrees.
    
    Returns:
    float, float: The right ascension (RA) and declination (Dec) of the moon in degrees.
    """
    # Convert the date and time to a datetime object
    date_time = datetime(year, month, day, hour, minute, second)
    
    # Calculate the Julian Day Number (JDN)
    jdn = 367 * year - int(7 * (year + int((month + 9) / 12)) / 4) + int(275 * month / 9) + day + 1721013.5 + (date_time.hour + date_time.minute / 60.0 + date_time.second / 3600.0) / 24.0
    
    # Calculate the Geocentric Ecliptic Longitude (L) and Latitude (B) of the moon
    T = (jdn - 2451545.0) / 36525.0
    L = 218.3164477 + 481267.88123421 * T - 0.0015786 * T**2 + T**3 / 538841.0 - T**4 / 65194000.0
    B = 5.13188 + 3.54648 * T - 0.00003 * T**2 - 0.00000367 * T**3
    
    # Convert the Geocentric Ecliptic Longitude and Latitude to Equatorial Coordinates
    epsilon = 23.439291 - 0.0130042 * T
    ra = math.atan2(math.cos(math.radians(epsilon)) * math.sin(math.radians(L)) - math.tan(math.radians(B)) * math.sin(math.radians(epsilon)), math.cos(math.radians(L)))
    dec = math.asin(math.sin(math.radians(B)) * math.cos(math.radians(epsilon)) + math.cos(math.radians(B)) * math.sin(math.radians(epsilon)) * math.sin(math.radians(L)))
    
    # Convert the RA and Dec to degrees
    ra_deg = math.degrees(ra)
    dec_deg = math.degrees(dec)
    
    return ra_deg, dec_deg

# Example usage
year = 2024
month = 4
day = 16
hour = 10
minute = 30
second = 0
latitude = 40.73  # New York City
longitude = -74.00  # New York City

ra, dec = calculate_moon_ra_dec(year, month, day, hour, minute, second, latitude, longitude)
print(f"On {month}/{day}/{year} at {hour:02d}:{minute:02d}:{second:02d}, the moon's right ascension is {ra:.2f} degrees and its declination is {dec:.2f} degrees.")
