#!/usr/bin/python3
# 
#  Pull in libraries
#from datetime import date
from datetime import date, datetime, timedelta
import numpy as np
import scipy
from astropy.coordinates import SkyCoord, EarthLocation, AltAz
from astropy.time import Time
import astropy.units as u
import astroplan
import math
import requests
import ephem
import pytz
from astronomica import *
import pandas as pd
#from astral import Astral
#from astral.geocoder import NominatimGeocoder

# ======================================================================================================================================
# Orient the code to where and when it is
#
latitude = 33.697472      # Latitude    
longitude = -117.72314     # Longitude
elevation = 22.0           # Elevation (meters)
bortle_class = 7          # integer see https://www.lightpollutionmap.info/
scopes_file = "/Users/christopherporter/Desktop/ASTRO/ATC/scopes.dat"
cameras_file = "/Users/christopherporter/Desktop/ASTRO/ATC/cameras.dat"
astro_db_file = "/Users/christopherporter/Desktop/ASTRO/ATC/astro_db.dat"
LOG = "/Users/christopherporter/Desktop/ASTRO/ATC/atc.log"
DEBUG = "YES"
# 
# Preferences
max_obj_magnitude = 14.0    # minimum brightness
min_obj_size = 9.0          # min object size (arcmin)
max_obj_size = 30.0         # max object size
min_time_up = 2.5       # hours that the object will be above min height
min_altitude = 30.0     # minimum altitude above horizon (deg)
max_moon_pct = 100.0    # maximum moon phase.  Set to 100.0 to disable
min_moon_angle = 25.0   # min angle between object and the moon (deg)

# ======================================================================================================================================
#
# Open LOG for appending
log_file = open(LOG, "a")

#
# Pull in Time info, different formats
from datetime import datetime
today_date = datetime.now()

#
# Open the header in the log file
log_file.write("============================  Start Run of Astro Target Cooser [ATC] ============================")
log_file.write(f"                                {today_date}\n")

year = today_date.year
month = today_date.month
day = today_date.day
hour = today_date.hour
minute = today_date.minute
second = today_date.second
if month < 10:
   date_yyyymmdd = str(year) + "-" + "0" + str(month) + "-" + str(day)
else:
   date_yyyymmdd = str(year) + "-" + str(month) + "-" + str(day)

# DEBUG : print(f"DATE (short):",date_yyyymmdd,"DATE:",date,"  YEAR:",year,"  MONTH:",month,"  DAY:",day,"  HOUR:",hour,"  MINUTE:",minute)
#
# MOON CALCULATIONS =====================================================================================================================================
#
#   Phase

def calculate_moon_phase(year, month, day, hour, latitude, longitude):
    """
    Calculates the phase of the moon for a given date, time, and location.
    
    Parameters:
    year (int): The 4-digit year.
    month (int): The month (1-12).
    day (int): The day of the month (1-31).
    hour (float): The hour of the day (0-23).
    latitude (float): The latitude of the location in degrees.
    longitude (float): The longitude of the location in degrees.
    
    Returns:
    str: The phase of the moon.
    """
    print("  ... Calculating moon phase")
    # Calculate the Julian Day Number
    jdn = int(365.25 * (year + 4716)) + int(30.6001 * (month + 1)) + day + (hour / 24) - 1524.5
    
    # Calculate the moon's age in days
    moon_age = (jdn - 2451549.5) / 29.530588853
    
    # Determine the moon phase
    if moon_age < 1:
        moon_phase = "New Moon"
    elif moon_age < 7.38:
        moon_phase = "Waxing Crescent"
    elif moon_age < 7.38 + 7.38:
        moon_phase = "First Quarter"
    elif moon_age < 14.77:
        moon_phase = "Waxing Gibbous"
    elif moon_age < 14.77 + 7.38:
        moon_phase = "Full Moon"
    elif moon_age < 22.15:
        moon_phase = "Waning Gibbous"
    elif moon_age < 22.15 + 7.38:
        moon_phase = "Last Quarter"
    else:
        moon_phase = "Waning Crescent"

    moon_pct = moon_age / 29.53
    
    return moon_phase, moon_pct 

moon_phase, moon_pct = calculate_moon_phase(year, month, day, hour, latitude, longitude)
# DEBUG print(f"On {month}/{day}/{year} at {hour:02.2f} hours, the moon will be in the {moon_phase} phase at {moon_pct:02.1f} % full.")
#
#   Moon Location in the sky  ---------------------------------------------------------------------------------------------------------

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
    print("  ... Calculating moon location")

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


ra, dec = calculate_moon_ra_dec(year, month, day, hour, minute, second, latitude, longitude)
# DEBUG print(f"On {month}/{day}/{year} at {hour:02d}:{minute:02d}:{second:02d}, the moon's right ascension is {ra:.2f} degrees and its declination is {dec:.2f} degrees.")

#
# Astronomy time available tonight - astro sunset to sunrise =============================================================================================k

import datetime
from astronomica import *

def get_astronomical_night():
    # Get the local sunset time as a datetime object
    #
    print("  ...Calculating astronomical night and sunrise")

    sun = Sun(latitude, longitude)
    sunset_time = sun.get_local_sunset_time(datetime.datetime.now())
    sunset_time_str = sunset_time.strftime("%H:%M")
    #print(f"\n\nThe local sunset time is {sunset_time_str}")

    # Calculate duration of the day in hours and minutes
    #day_duration = datetime.timedelta(hours=24) - datetime.timedelta(hours=sunset_time.hour, minutes=sunset_time.minute)
    sunrise_time = sun.get_local_sunrise_time(datetime.datetime.now())
    day_duration =  datetime.timedelta(hours=sunset_time.hour, minutes=sunset_time.minute) - datetime.timedelta(hours=sunrise_time.hour, minutes=sunrise_time.minute)
    print(f"      DEBUG: day_duration: {day_duration}")

    day_duration_hours, day_duration_minutes = divmod(day_duration.seconds // 60, 60)

    # Calculate the duration of astronomical twilight as a timedelta object
    twilight_duration = datetime.timedelta(minutes=day_duration.total_seconds() * 0.2666 / 60)

    # Calculate the duration of astronomical night as a timedelta object
    night_duration =  day_duration - twilight_duration

    # Calculate the time of astronomical night as a datetime object
    night_time = sunset_time + (twilight_duration / 2.4)

    # Format the time of astronomical night as a string in the format "HH:MM"
    night_time_str = night_time.strftime("%H:%M")

    print(f"    The start time of astronomical night tonight is {night_time_str}")

    tomorrow_date = datetime.date.today() + datetime.timedelta(days=1)
    sunrise_time = sun.get_local_sunrise_time(tomorrow_date)
    sunrise_time_str = sunrise_time.strftime("%H:%M")

    dawn_time = sunrise_time - (twilight_duration / 2.4)
    dawn_time_str = dawn_time.strftime("%H:%M")
    
    #print(f"\nTotal night duration is: {night_duration}")
    #print(f"Twilight duration is: {twilight_duration}")

    print(f"    The time of astronomical night ends tomrrow is {dawn_time_str}")
    print(f"      DEBUG: The local sunrise time tomorrow is {sunrise_time_str}")

get_astronomical_night()

#
# READ INPUT FILES ===================================================================================================================
#
import pandas as pd
import math

# READ INPUT FILES ===================================================================================================================
#
#print("SCOPES")
print("  ...Reading SCOPES.dat")
scopes_dataframe = pd.read_csv('/Users/christopherporter/Desktop/ASTRO/ATC/SCOPES.dat')
#print(scopes_dataframe)
#print("\n\n")

# Scope Fields
#  1- Label, 
#  2- Focal Length (mm), 
#  3- aperture (mm), 
#  4- Reducer (factor), 
#  5- Reducer On Scope (YES/NO)

print("  ...Reading CAMERAS.dat")
#print("CAMERAS")
cameras_dataframe = pd.read_csv('/Users/christopherporter/Desktop/ASTRO/ATC/CAMERAS.dat')
#print(cameras_dataframe)
#print("\n\n")

# Camera Fields
#  1- Label, 
#  2- Vertical Pixels, 
#  3- Horizontal Pixels, 
#  4- Pixel Size (Microns), 
#  5- Type (Mono/OSC), 
#  6- QE (%), 
#  7- Cooled (YES,NO)
#  8- Imaging (YES,NO)

print("  ..Reading objects DB")
objects_dataframe = pd.read_csv('/Users/christopherporter/Desktop/ASTRO/ATC/astro_db.dat')
print("  ..Done Reading objects DB")
# print out the stats about the dataframe
#objects_dataframe.info()

# DB Fields
#  1- Label 1, 
#  2- Label 2, 
#  3- Constellation, 
#  4- Magnitude, 
#  5- Season, 
#  6- Size (arcmin), 
#  7- RA, 
#  8- Dec, 
#  9- Type, 
#  10- Notes

# =======================================================================================================================================================
# Calculate Filed of View for each scope and camera combination
# Loop through cameras
for camera, row in cameras_dataframe.iterrows():
    if str(cameras_dataframe.loc[camera,'imaging']) == "NO":
        print(f"      DEBUG: Camera: {cameras_dataframe.loc[camera,'label']} - Skipping, not a DSO imaging camera")
        continue
    # Loop through scopes
    for scope, row in scopes_dataframe.iterrows():
        # FOV (in arc-min) = 3436 * D / L
        # D is the sensor dimension in mm
        # L is the FL in mm
        vert_sensor_size = 0.001 * cameras_dataframe.loc[camera,'pix_size'] * cameras_dataframe.loc[camera,'vert_pix']
        horz_sensor_size = 0.001 * cameras_dataframe.loc[camera,'pix_size'] * cameras_dataframe.loc[camera,'horz_pix']
        field_of_view =  (3436.0 * math.sqrt(vert_sensor_size**2 + horz_sensor_size**2)) / scopes_dataframe.loc[scope,'focal_length_mm']
        print(f"      DEBUG: Camera: {cameras_dataframe.loc[camera,'label']},  Scope: {scopes_dataframe.loc[scope,'label']}, field of view (arcmin): {field_of_view:.3f}")

        # Arcsec per pixel =(PS/FL)∗206.265
        # PS is the pixel size (microns)
        # FL is the focal length (mm)
        arcsec_per_pixel = (206.265 * cameras_dataframe.loc[camera,'pix_size']) / scopes_dataframe.loc[scope,'focal_length_mm']
        print(f"      DEBUG: Camera: {cameras_dataframe.loc[camera,'label']},  Scope: {scopes_dataframe.loc[scope,'label']},     arcsec per pixel: {arcsec_per_pixel:.3f}")

        #tmp_data = {cameras_dataframe.loc[camera,0]: [scopes_dataframe.loc[scope,0], field_of_view, arcsec_per_pixel]}

 # =========================================================================================================================================
 # FILTER - go throough the DB of objects and eliminate those that are not possible or don't satisfy criteria
 #
 #    Filters: 
 #       1) higher than min magnitude
 #       2) higher than max size
 #       3) lower than min size
 #       4) not visible long enough
 #       5) min angle to the moon not reached
 #       6) min elevation insufficient
 #       7) moon phase too high 
 for object, row in objects_dataframe.iterrows():
    # magnitude
    object_magnitude = objects_dataframe.loc[object,'magnitude']
    if object_magnitude > max_obj_magnitude:
        continue

    # max size
    object_size = objects_dataframe.loc[object,'size']
    if object_size > max_obj_size:
        continue
    if object_size < min_obj_size:
        continue

    # max moon phase
    if moon_pct > max_moon_pct:
        continue

    # min moon angle to the object

    observer = ephem.Observer()
    observer.lat = str(latitude)
    observer.lon = str(longitude)

    moon = ephem.Moon()
    moon.compute(observer)

    moon_RA = (moon.ra * 12) / ephem.pi
    moon_DEC = (moon.dec * 180) / ephem.pi

    obj_label = objects_dataframe.loc[object,'label']
    obj_RA = objects_dataframe.loc[object,'RA']
    obj_DEC = objects_dataframe.loc[object,'DEC']

    print(f"   DEBUG: From DB Object: {obj_label}, RA: {obj_RA}, DEC: {obj_RA}")

    obj_center = SkyCoord.from_name(str(obj_label))
    obj_RA = obj_center.ra   # in deg
    obj_DEC = obj_center.dec # in deg

    print(f"   DEBUG: From AstroPy.SkyCoord Object: {obj_label}, RA: {obj_RA}, DEC: {obj_RA}")

    # θ = arccos(sin(Dec1) * sin(Dec2) + cos(Dec1) * cos(Dec2) * cos(RA1 - RA2))
    angle_to_moon = math.acos(sin(moon_DEC) * math.sin(obj_DEC) + math.cos(moon_DEC) * math.cos(obj_DEC) * math.cos(moon_RA - obj_RA))
    print(f"   DEBUG: Angle between {obj_label} and Moon is {angle_to_moon}")
    
    if angle_to_moon < min_moon_angle:
        continue