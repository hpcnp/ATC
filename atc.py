#!/usr/bin/python3
# 
#  Pull in libraries
#from datetime import date
from datetime import date, datetime, timedelta
import re
from re import split
import numpy as np
import scipy
from scipy.interpolate import interp2d
from astropy.coordinates import SkyCoord, EarthLocation, AltAz, ICRS, get_sun, get_moon
from astropy.time import Time
import astropy.units as u
import astroplan
import math
import requests
import ephem
import pytz
from astronomica import *
import pandas as pd
from tabulate import tabulate 
#from astral import Astral
#from astral.geocoder import NominatimGeocoder

# ======================================================================================================================================
# Orient the code to where and when it is
#
latitude = 33.697472      # Latitude    
longitude = -117.72314     # Longitude
elevation = 22.0           # Elevation (meters)
bortle_class = 7          # integer see https://www.lightpollutionmap.info/
scopes_file = "/Users/christopherporter/Desktop/ASTRO/ATC/SCOPES.dat"
cameras_file = "/Users/christopherporter/Desktop/ASTRO/ATC/CAMERAS.dat"
astro_db_file = "/Users/christopherporter/Desktop/ASTRO/ATC/astro_db4.dat"
LOG = "/Users/christopherporter/Desktop/ASTRO/ATC/atc.log"
DEBUG = 0               # DEBUG 0,1,2,3   0 - no debug output, 1-3 more and more verbose output
# 
# Preferences
max_obj_magnitude = 15.0    # minimum brightness
min_obj_size = 20.0          # min object size (arcmin)
max_obj_size = 100.0         # max object size (arcmin)
min_time_up = 3.5       # hours that the object will be above min height
min_altitude = 30.0     # minimum altitude above horizon (deg)
max_moon_pct = 100.0    # maximum moon phase.  Set to 100.0 to disable
max_moon_pct_broad = 25.0    # maximum moon phase for broadband target.  Set to 100.0 to disable
min_moon_angle = 25.0   # min angle between object and the moon (deg)
timezone_offset = 7     # local timezone offset to GMT
FOV_ideal_frac = 0.75   # the ideal fraction of the field of view (used for recommendations)
obj_count = 15          # number of objects in the output table
preferred_types = ['galaxy', 'emission nebula', 'reflection nebula', 'bright nebula', 'supernova']
discarded_types = ['star', 'open clus', 'asterism', 'diffuse', 'existant', 'open cluster']
lp_constants = ['1','3','25','40','100'] # Filters: none, LRGB/UVIR, 12nm, 7nm, 3nm
narrow_types = ['nebula', 'supernova']
output_RA = 'hms'       # Output RA in hours:minutes:seconds [hms] or decimal degrees [dd]
output_DEC = 'dms'      # Output DEC in degrees:minutes:seconds [dms] or decimal degrees [dd]
output_RECS = 'YES'     # Recommend which scope, reducer, filter, exposure range to use per target

# ======================================================================================================================================
#
# Open LOG for appending
log_file = open(LOG, "a")

# ======================================================================================================================================
#  Functions 

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#   Convert decimal degress to hour-minute-second
#   Format:  HH:MM:SS.SSS
def decdeg_to_hms(coord_degrees, pretty_print=None, ndp=2):
    """Convert from decimal degrees to degrees, minutes, seconds."""
    import math

    total_seconds = (coord_degrees / 360.0) * (24.0 * 3600.0)
    coord_hours = int(total_seconds / 3600.0)

    rem_seconds = total_seconds - coord_hours * 3600.0
    coord_minutes = int(rem_seconds / 60.0)
    coord_seconds = total_seconds - coord_hours * 3600.0 - coord_minutes * 60.0
    #coord_hms = str(coord_hours) + ":" + str(coord_minutes) + ":" + str(coord_seconds)
    coord_hms = '{}:{}:{:.3f}'.format(coord_hours, coord_minutes, coord_seconds)

    if DEBUG >= 2:
        log_file.write(f"   Convert DecDeg to HMS: Deg: {coord_degrees} whcih is {coord_hours} hours, {coord_minutes} min, and {coord_seconds} sec\n")
    if DEBUG >= 3:
        print(f"     Convert DecDeg to HMS: Deg: {coord_degrees} whcih is {coord_hours} hours, {coord_minutes} min, and {coord_seconds} sec\n")

    if pretty_print:
        if pretty_print=='latitude':
            hemi = 'N' if d>=0 else 'S'
        elif pretty_print=='longitude':
            hemi = 'E' if d>=0 else 'W'
        else:
            hemi = '?'
        return '{d:d}° {coord_minutes:coord_degrees}′ {coord_seconds:.{ndp:coord_degrees}f}″ {hemi:1s}'.format(
                    coord_degress=abs(coord_degrees), coord_minutes=m, coord_seconds=seconds, hemi=hemi, ndp=ndp)
    return coord_hours, coord_minutes, coord_seconds, coord_hms

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#   Convert decimal degrees to degrees-minutes-seconds
#   Format:  DD:MM:SS.SSS
def decdeg_to_dms(coord_degrees, pretty_print=None, ndp=3):
    """Convert from decimal degrees to degrees, minutes, seconds."""

    coord_minutes, coord_seconds = divmod(abs(coord_degrees)*3600, 60)
    coord_degrees, coord_minutes = divmod(coord_minutes, 60)

    if coord_degrees < 0:
        coord_degrees = -coord_degrees
    coord_degrees, coord_minutes = int(coord_degrees), int(coord_minutes)

    #coord_dms = str(coord_degrees) + ":" + str(coord_minutes) + ":" + str(coord_seconds)
    coord_dms = '{}:{}:{:.3f}'.format(coord_degrees, coord_minutes, coord_seconds)

    if pretty_print:
        if pretty_print=='latitude':
            hemi = 'N' if d>=0 else 'S'
        elif pretty_print=='longitude':
            hemi = 'E' if d>=0 else 'W'
        else:
            hemi = '?'
        return '{d:d}° {coord_minutes:coord_degrees}′ {coord_seconds:.{ndp:coord_degrees}f}″ {hemi:1s}'.format(
                    coord_degress=abs(coord_degrees), coord_minutes=m, coord_seconds=seconds, hemi=hemi, ndp=ndp)
    return coord_degrees, coord_minutes, coord_seconds, coord_dms

def decdeg_to_rad(degrees):
    import math
    radians = degrees * math.pi / 180.0

    return radians

def rad_to_decdeg(coord_radians):
    import math
    coord_degrees = coord_radians * 180.0 / math.pi 

    return coord_degrees

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#++++++++++++
#   TO DO
#++++++++++++
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#   Convert degress-minutes-seconds to decimal degrees
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#   Convert hour-minute-second to decimal degrees
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# Pull in Time info, different formats
from datetime import datetime
today_date = datetime.now()

#
# Open the header in the log file
log_file.write("============================  Start Run of Astro Target Chooser [ATC] ============================\n")
log_file.write(f"                                {today_date}\n")

year = today_date.year
month = today_date.month
day = today_date.day
hour = today_date.hour
minute = today_date.minute
second = today_date.second

observing_location = EarthLocation(lat=latitude*u.deg, lon=longitude*u.deg)
#observing_time = Time(datetime.datetime.utcnow(), scale='utc', location=observing_location)
#observing_time = Time(datetime.now(datetime.UTC), location=observing_location)
#local_siderial_time = observing_time.sidereal_time('mean')

if month < 10:
   date_yyyymmdd = str(year) + "-" + "0" + str(month) + "-" + str(day)
else:
   date_yyyymmdd = str(year) + "-" + str(month) + "-" + str(day)

if DEBUG >= 2:
    log_file.write(f"   DATE (short): {date_yyyymmdd} DATE: {date}  YEAR: {year}  MONTH: {month}  DAY: {day}  HOUR: {hour}  MINUTE: {minute}\n")
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
if DEBUG >= 1:
    log_file.write(f"   On {month}/{day}/{year} at {hour:02.2f} hours, the moon will be in the {moon_phase} phase at {moon_pct:02.1f} % full.\n")
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
    print("  ... Calculating moon location for tonight")

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
# Astronomy time available tonight - astro sunset to sunrise =============================================================================================

import datetime
from astronomica import *

def get_astronomical_night():
    # Get the local sunset time as a datetime object
    #
    print("  ...Calculating astronomical night and sunrise")

    sun = Sun(latitude, longitude)
    sunset_time = sun.get_local_sunset_time(datetime.datetime.now())
    sunset_time_str = sunset_time.strftime("%H:%M")
    if DEBUG >=2:
        log_file.write(f"   The local sunset time is {sunset_time_str}\n")

    # Calculate duration of the day in hours and minutes
    #day_duration = datetime.timedelta(hours=24) - datetime.timedelta(hours=sunset_time.hour, minutes=sunset_time.minute)
    sunrise_time = sun.get_local_sunrise_time(datetime.datetime.now())
    day_duration =  datetime.timedelta(hours=sunset_time.hour, minutes=sunset_time.minute) - datetime.timedelta(hours=sunrise_time.hour, minutes=sunrise_time.minute)
    if DEBUG >= 1:
        log_file.write(f"   Day_duration: {day_duration}\n")

    day_duration_hours, day_duration_minutes = divmod(day_duration.seconds // 60, 60)

    # Calculate the duration of astronomical twilight as a timedelta object
    twilight_duration = datetime.timedelta(minutes = day_duration.total_seconds() * 0.2666 / 60.0 )

    # Calculate the duration of astronomical night as a timedelta object
    night_duration =  day_duration - twilight_duration

    # Calculate the time of astronomical night as a datetime object
    night_time = sunset_time + (twilight_duration / 2.4)

    # Format the time of astronomical night as a string in the format "HH:MM"
    night_time_str = night_time.strftime("%H:%M")

    print(f"     The start time of astronomical night tonight is {night_time_str}")
    if DEBUG >= 1:
        log_file.write(f"     The start time of astronomical night tonight is {night_time_str}\n")

    tomorrow_date = datetime.date.today() + datetime.timedelta(days=1)
    sunrise_time = sun.get_local_sunrise_time(tomorrow_date)
    sunrise_time_str = sunrise_time.strftime("%H:%M")

    dawn_time = sunrise_time - (twilight_duration / 2.4)
    dawn_time_str = dawn_time.strftime("%H:%M")
    
    #print(f"\nTotal night duration is: {night_duration}")
    #print(f"Twilight duration is: {twilight_duration}")

    print(f"     The time of astronomical night ends tomrrow is {dawn_time_str}")
    if DEBUG >= 1:
        log_file.write(f"   The time of astronomical night ends tomrrow is {dawn_time_str}\n")
        log_file.write(f"   The local sunrise time tomorrow is {sunrise_time_str}\n")

get_astronomical_night()

#
# READ INPUT FILES ===================================================================================================================
#
import pandas as pd
import math

# READ INPUT FILES ===================================================================================================================
#
print("  ...Reading SCOPES.dat")
scopes_dataframe = pd.read_csv(scopes_file)

if DEBUG >= 3:
    log_file.write("SCOPES\n")
    print(scopes_dataframe)
    log_file.write("\n")

# Scope Fields
#  1- Label, 
#  2- Focal Length (mm), 
#  3- aperture (mm), 
#  4- Reducer (factor), 
#  5- Reducer On Scope (YES/NO)

print("  ...Reading CAMERAS.dat")
cameras_dataframe = pd.read_csv(cameras_file)

if DEBUG >= 3:
    log_file.write("CAMERAS\n")
    print(cameras_dataframe)
    log_file.write("\n")

# Camera Fields
#  1- Label, 
#  2- Vertical Pixels, 
#  3- Horizontal Pixels, 
#  4- Pixel Size (Microns), 
#  5- Type (Mono/OSC), 
#  6- QE (%), 
#  7- Cooled (YES,NO)
#  8- Imaging (YES,NO)
#  9- Full well (keV)
# 10- Unity Gain
# 11- Read noise low
# 12- Read noise high

print("  ...Reading objects DB")
objects_dataframe = pd.read_csv('/Users/christopherporter/Desktop/ASTRO/ATC/astro_db3.dat')
print("  ....Done Reading objects DB")

if DEBUG >= 3:
    # print out the stats about the dataframe
    log_file.write("OBJECTS\n")
    objects_dataframe.info()
    log_file.write("\n")

# DB Fields
#  0- Label 1, 
#  1- Label 2, 
#  2- Constellation, 
#  3- Magnitude, 
#  4- Season, 
#  5- Size (arcmin), 
#  6- RA, 
#  7- Dec, 
#  8- Type, 
#  9- Notes

# =======================================================================================================================================================
# Calculate Field of View for each scope and camera combination
# Loop through cameras
for camera, row in cameras_dataframe.iterrows():
    if str(cameras_dataframe.loc[camera,'imaging']) == "NO":
        if DEBUG >= 1:
            log_file.write(f"   Camera: {cameras_dataframe.loc[camera,'label']} - Skipping, not a DSO imaging camera\n")
        if DEBUG >= 2:
            print(f"     Camera: {cameras_dataframe.loc[camera,'label']},  - Skipping, not a DSO imaging camera")

        continue
    # Loop through scopes
    for scope, row in scopes_dataframe.iterrows():

        # FOV (in arc-min) = 3436 * D / L
        # D is the sensor dimension in mm
        # L is the FL in mm
        vert_sensor_size = 0.001 * cameras_dataframe.loc[camera,'pix_size'] * cameras_dataframe.loc[camera,'vert_pix']
        horz_sensor_size = 0.001 * cameras_dataframe.loc[camera,'pix_size'] * cameras_dataframe.loc[camera,'horz_pix']
        field_of_view =  (3436.0 * math.sqrt(vert_sensor_size**2 + horz_sensor_size**2)) / scopes_dataframe.loc[scope,'focal_length_mm']

        if DEBUG >=1:
            log_file.write(f"   Camera: {cameras_dataframe.loc[camera,'label']},  Scope: {scopes_dataframe.loc[scope,'label']}, field of view (arcmin): {field_of_view:.3f}\n")
        if DEBUG >=2:
            print(f"     Camera: {cameras_dataframe.loc[camera,'label']},  Scope: {scopes_dataframe.loc[scope,'label']}, field of view (arcmin): {field_of_view:.3f}")

        # Arcsec per pixel =(PS/FL)∗206.265
        # PS is the pixel size (microns)
        # FL is the focal length (mm)
        arcsec_per_pixel = (206.265 * cameras_dataframe.loc[camera,'pix_size']) / scopes_dataframe.loc[scope,'focal_length_mm']

        if DEBUG >= 1:
            log_file.write(f"   Camera: {cameras_dataframe.loc[camera,'label']},  Scope: {scopes_dataframe.loc[scope,'label']},     arcsec per pixel: {arcsec_per_pixel:.3f}\n")
            log_file.write("\n")
        if DEBUG >=2:
            print(f"     Camera: {cameras_dataframe.loc[camera,'label']},  Scope: {scopes_dataframe.loc[scope,'label']},  arcsec per pixel: {arcsec_per_pixel:.3f}")

        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # REVISIT
        # Calculations for exposure based on bortle class, scope, and camera
        #if scopes_dataframe.loc[scope,'reducer_on_scope'] == 'YES':
        #    f_stop = (scopes_dataframe.loc[scope,'reducer_factor'] * scopes_dataframe.loc[scope,'focal_length_mm']) / scopes_dataframe.loc[scope,'aperture_mm']
        #else:
        #    f_stop = (scopes_dataframe.loc[scope,'focal_length_mm']) / scopes_dataframe.loc[scope,'aperture_mm']

        #for camera, row in cameras_dataframe.iterrows():
        #    if camaras_dataframe.loc[camera,'imaging'] == 'NO':
        #        continue
        #    lp_mod_param = 1.0 + (cameras_dataframe.loc[camera,'QE'] - 0.5) / ( )

        #min_exposure = (10.0 * cameras_dataframe.loc[camera,'read_noise_l']**2)
        #max_exposure = (10.0 * cameras_dataframe.loc[camera,'read_noise_h']**2)



 # =========================================================================================================================================
 # FILTER - go throough the DB of objects and eliminate those that are not possible or don't satisfy criteria
 #
 #    Object Filters: 
 #       1) higher than min magnitude
 #       2) higher than max size
 #       3) lower than min size
 #       4) not visible long enough
 #       5) min angle to the moon not reached
 #       6) min elevation insufficient
 #       7) moon phase too high 

objects_passed = 0

#
# Create a DataFrame for storing filtered objects and their scores
filtered_obj_dataframe = pd.DataFrame({
    'label1': [],
    'label2': [],
    'constellation': [],
    'magnitude': [],
    'season': [],
    'size': [],
    'ra': [],
    'dec': [],
    'type': [],
    'notes': [],
    'score1': [],
    'score2': [],
    'score3': [],
    'score4': [],
    'score5': [],
    'totalscore': []
})

#
# Create a DataFrame for storing recommendations for scope, camera, reducer, filter, exposure
reccommendations_df = pd.DataFrame({
    'label1': [],
    'total_score': [],
    'scope': [],
    'camera': [],
    'reducer': [],
    'obj_size': [],
    'field_of_view': [],
    'fov_frac': [],
    'fov_err': [],
    'filter': [],
    'min_exp': [],
    'max_exp': []
})

print("  ...Looping through object DB and eliminating those that violate constraints")

for object, row in objects_dataframe.iterrows():
    # ---------------------------------------------------------------------------------------------------------------------------------
    #   Apparent magnitude filter

    if DEBUG >= 2:
        log_file.write(f"   Filter: Checking minimum apparent brightness (magnitude < {max_obj_magnitude})\n")
    if DEBUG >= 3:
        print(f"      Filter: Checking minimum apparent magnitude < {max_obj_magnitude}\n")

    object_magnitude = objects_dataframe.loc[object,'magnitude']
    obj_label = objects_dataframe.loc[object,'label1']

    obj_score1 = (max_obj_magnitude - float(object_magnitude)) * 2.0

    if float(object_magnitude) > max_obj_magnitude:
        if DEBUG >= 1:
            log_file.write(f"   ELIMINATED {obj_label} for being too feint at magnitude {object_magnitude}\n")
        if DEBUG >= 2:
            print(f"      ELIMINATED {obj_label} for being too feint at magnitude {object_magnitude}")

        continue

    # ---------------------------------------------------------------------------------------------------------------------------------
    #   Size filter - operating in arcminutes

    if DEBUG >= 2:
        log_file.write(f"   Filter: Checking size of object between minimum ({min_obj_size}) and maximum {max_obj_size} arcmin.\n")
    if DEBUG >= 3:
        print(f"      Filter: Checking size of object between minimum {min_obj_size} and maximum {max_obj_size} arcmin.\n")

    object_size = objects_dataframe.loc[object,'size']

    obj_score2 = ( float(object_size) / max_obj_size ) * 10.0
    if float(object_size) > max_obj_size:
        if DEBUG >= 1:
            log_file.write(f"   ELIMINATED {obj_label} for being too large at {object_size} arcmin.\n")
        if DEBUG >= 2:
            print(f"      ELIMINATED {obj_label} for being too large at {object_size} arcmin.")
        continue

    if float(object_size) < min_obj_size:
        if DEBUG >= 1:
            log_file.write(f"   ELIMINATED {obj_label} for being too small at {object_size} arcmin.\n")
        if DEBUG >= 2:
            print(f"      ELIMINATED {obj_label} for being too small at {object_size} arcmin.")
        continue

    # ---------------------------------------------------------------------------------------------------------------------------------
    #   max moon phase filter
    if DEBUG >= 2:
        log_file.write(f"   Filter: Checking maximum moon phase ({max_moon_pct}%).\n")
    if DEBUG >= 3:
        print(f"      Filter: Checking moon phase ({max_moon_pct}).\n")

    obj_score5 = 7.5 * math.cos((math.pi / 2.0) * (moon_pct / 100.0))

    if moon_pct > 40.0:
        # narrow band targets get bumps up in score
        # broad band targets get bumps down in score
        for n_type in narrow_types:
            obj_type = objects_dataframe.loc[object,'type']

            if (n_type.lower() in obj_type.lower()):
                obj_score5 = 10.0 * math.cos((math.pi / 2.0) * (moon_pct / 100.0)) 
                break

    if moon_pct > max_moon_pct:
        if DEBUG >= 1:
            log_file.write(f"   ELIMINATED {obj_label} because the moon phase {moon_pct} % is higher than threshold {max_moon_pct}%\n")
        if DEBUG >= 2:
            print(f"      ELIMINATED {obj_label} because the moon phase {moon_pct} % is higher than threshold {max_moon_pct}%")
        # ********* 
        # This should be a break statement because the object doesn't change the moon phase
        # Add something for narrowband vs broadband imaging
        continue

    # ---------------------------------------------------------------------------------------------------------------------------------
    #  min moon angle to the object filter

    if DEBUG >= 2:
        log_file.write(f"   Filter: Checking Minimum moon angle({min_moon_angle} deg).\n")
    if DEBUG >= 3:
        print(f"      Filter: Checking Minimum moon angle ({min_moon_angle} deg)\n")

    observer = ephem.Observer()
    observer.lat = str(latitude)
    observer.lon = str(longitude)

    moon = ephem.Moon()
    moon.compute(observer)

    moon_RA = (moon.ra * 12) / ephem.pi
    moon_DEC = (moon.dec * 180) / ephem.pi

    try:
        obj_label = objects_dataframe.loc[object,'label1']
        obj_center = SkyCoord.from_name(str(obj_label))
        obj_RA = obj_center.ra.value   # in deg
        obj_DEC = obj_center.dec.value # in deg

        if DEBUG >= 1:
            log_file.write(f"   From AstroPy.SkyCoord Object: {obj_label}, RA: {obj_RA:.3f} deg., DEC: {obj_DEC:.3f} deg.\n")
        if DEBUG >= 3:
            print(f"    From AstroPy.SkyCoord Object: {obj_label}, RA: {obj_RA:.3f} deg., DEC: {obj_DEC:.3f} deg.")

    except:
        obj_label = objects_dataframe.loc[object,'label1']
        obj_RA = objects_dataframe.loc[object,'RA']
        obj_DEC = objects_dataframe.loc[object,'DEC']

        if DEBUG >= 1:
            log_file.write(f"   From DB Object: {obj_label}, RA: {obj_RA} deg., DEC: {obj_DEC:.4f} deg.\n")
        if DEBUG >= 2:
            print(f"    From DB Object: {obj_label}, RA: {obj_RA} deg., DEC: {obj_DEC:.4f} deg.")

    #
    # Python math.<trig_function> uses radians instead of degrees - so convert the quantities to radians
    obj_RA_rad = ( math.pi / 180.0 ) * obj_RA
    obj_DEC_rad = ( math.pi / 180.0 ) * obj_DEC
    moon_RA_rad = ( math.pi / 180.0 ) * moon_RA
    moon_DEC_rad = ( math.pi / 180.0 ) * moon_DEC

    # θ = arccos(sin(Dec1) * sin(Dec2) + cos(Dec1) * cos(Dec2) * cos(RA1 - RA2))
    angle_to_moon_rad = math.acos(math.sin(moon_DEC_rad) * math.sin(obj_DEC_rad) + math.cos(moon_DEC_rad) * math.cos(obj_DEC_rad) * math.cos(moon_RA_rad - obj_RA_rad))
    angle_to_moon = (180.0 / math.pi) * angle_to_moon_rad # convert back to degrees

    # calculatue the score - where 90 deg is a 10.
    obj_score3 = 10.0 * math.cos(angle_to_moon_rad)

    if DEBUG >= 2:
        log_file.write(f"   From Ephem: Moon, RA: {moon_RA:.4f} deg., DEC: {moon_DEC:.4f} deg., Fullness: {moon_pct:.1f}%\n")
        log_file.write(f"     Angle between {obj_label} and Moon is {angle_to_moon:.3f}\n")
    if DEBUG >= 3:
        print(f"      Angle between {obj_label} and Moon is {angle_to_moon:.3f}")
    
    if angle_to_moon < min_moon_angle:
        if DEBUG >= 1:
            log_file.write(f"   ELIMINATED: Angle between the moon and the object {angle_to_moon} deg. is too low. Threshold is {min_moon_angle} deg\n")
        if DEBUG >= 2:
            print(f"      ELIMINATED: Angle between the moon and the object {angle_to_moon} deg. is too low. Threshold is {min_moon_angle} deg")
        continue

    # ---------------------------------------------------------------------------------------------------------------------------------
    #   Object type filter

    if DEBUG >= 2:
        log_file.write(f"   Filter: Checking that the type of object is of interest\n")
    if DEBUG >= 3:
        print("      Filter: Checking that the type of object is of interest")

    # ID a narrowband vs. broadband target (e.g. is it a nebula or not)
    obj_type = objects_dataframe.loc[object,'type']
    if obj_type.find('Nebula') != -1:
        obj_type_filter = 'narrowband'
    else:
        obj_type_filter = 'broadband'

    # Reject broadband targets when the moon is too full
    if obj_type_filter == 'broadband' and moon_pct >= max_moon_pct_broad:
        if DEBUG >= 1:
            log_file.write(f"    ELIMINATED {obj_label} because it is a broadband target and the moon is too full {moon_pct}%\n")
        if DEBUG >= 2:
            print(f"      ELIMINATED {obj_label} because it is a broadband target and the moon is too full {moon_pct}%\n")
        continue 

    skip_obj = "NO"
    for dtype in discarded_types:
        if (dtype.lower() in obj_type.lower()):
            if DEBUG >= 1:
                log_file.write(f"    Discarded {obj_label} because it is a discard type\n")
            if DEBUG >= 2:
                print(f"      ELIMINATED {obj_label} with type: {obj_type} because this is a discarded type {dtype}")
            skip_obj = "YES"
            break
    if (skip_obj == "YES"):
        continue

    # ---------------------------------------------------------------------------------------------------------------------------------
    # Duration above minimum altitude

    if DEBUG >= 2:
        log_file.write(f"   Filter: Checking minimum visible time ({min_time_up} hours) above minimum altitude ({min_altitude} deg.)\n")
    if DEBUG >= 3:
        print(f"      Filter: Duration above minimum ({min_time_up}) altitude")

    hours_visible = 0

    for i in range(20+timezone_offset, 29+timezone_offset):
        if i > 23:
            mil_hour = i - 24
        else:
            mil_hour = i
        mil_time = str(mil_hour) + ":00:00"

        if i > 23:
            query_date = str(year) + "-" + str(month) + "-" + str(day+1)
        else:
            query_date = str(year) + "-" + str(month) + "-" + str(day)

        query_time = str(query_date) + " " + str(mil_time)

        # Calculate the altitude of the object at this time
        try:
            observer_location = EarthLocation(lat=latitude*u.deg, lon=longitude*u.deg)
            time = Time(query_time)
            aspy_object = SkyCoord.from_name(obj_label)
            alt_az = aspy_object.transform_to(AltAz(obstime=time, location=observer_location))
            altitude = alt_az.alt.value # in deg
            azimuth = alt_az.az.value # in deg
        except:
            # Use ephem to convert RA,DEC, Date, Time, and Location to Alt-Az

            if DEBUG >= 1:
                log_file.write(f"       Alternate Method for calculating Alt-Az\n")
                #print(f"       Alternate method for calculating Alt-Az")

            obj = ephem.FixedBody()
            obj._ra = str(obj_RA)
            obj._dec = str(obj_DEC)

            time = Time(query_time)
            observer = ephem.Observer()
            observer.date = query_time
            observer.lon = longitude
            observer.lat = latitude
            observer.elevation = elevation 

            obj.compute(observer)

            # Change Altitude from Ephem in DD:MM:SS to DD.DDDDD
            tmp_alt = str(obj.alt)
            alt_deg, alt_min, alt_sec = tmp_alt.split(':')

            if float(alt_deg) >= 0.0:
                altitude = float(alt_deg) + float(alt_min)/60.0 + float(alt_sec)/3600.0
            else:
                altitude = float(alt_deg) - float(alt_min)/60.0 - float(alt_sec)/3600.0

            # Change Azimuth from Ephem in DD:MM:SS to DD.DDDDD
            tmp_az = str(obj.az)
            az_deg, az_min, az_sec = tmp_az.split(':')

            #
            # Not needed (yet) but calcualte local sidereal time and hour angle
            # local_sidereal_time = 
            # hour_angle =

            if float(az_deg) >= 0.0:
                altitude = float(az_deg) + float(az_min)/60.0 + float(az_sec)/3600.0
            else:
                altitude = float(az_deg) - float(az_min)/60.0 - float(az_sec)/3600.0


        if altitude >= min_altitude:
            hours_visible = hours_visible + 1
            if DEBUG >= 2:
                log_file.write(f"       Above min ALT {min_altitude} at {time} GMT. Added to hours_visible {hours_visible}\n")
    
    obj_score4 = hours_visible + 1.0

    if hours_visible < min_time_up:
        if DEBUG >= 1:
            log_file.write(f"    Object RA: {obj_RA}, Object DEC: {obj_DEC}, Object ALT: {altitude}, Object AZ: {azimuth}\n")
            log_file.write(f"     ELIMINATED: {obj_label} not visible ({hours_visible} hrs) for long enough {min_time_up} hrs\n")
        if DEBUG >= 2:
            log_file.write(f"       Loop: At {time} {obj_label} is at: Altitude: {altitude:3f} Azimuth: {azimuth:3f} hours_visible: {hours_visible}\n")
        continue

    # ---------------------------------------------------------------------------------------------------------------------------------
    #   Count the passed objects
    
    objects_passed = objects_passed + 1
    print(f"  Object passing filter: {obj_label}")

    # ---------------------------------------------------------------------------------------------------------------------------------
    #   append to filtered_obj_dataframe

    # DB Fields
    #  0- Label 1, 
    #  1- Label 2, 
    #  2- Constellation, 
    #  3- Magnitude, 
    #  4- Season, 
    #  5- Size (arcmin), 
    #  6- RA, 
    #  7- Dec, 
    #  8- Type, 
    #  9- Notes
    #  10- Score1:  max_magnitude - object_magnitude 
    #  11- Score2:  10% of the fraction of the max_object_size
    #  12- Score3:  10x the cosine of the angle to the moon (so 90 deg. is a "10")
    #  13- Score4:  Hours visible tonight + 1 
    #  14- Score5:  Moon percent full > 40% preference to narrowband targets, otherwise no preference
    #  15- Reccommended scope
    #  16- Recommended reducer 
    #  17- Recommended camera
    #  18- Recommended filter

    obj_label2 = objects_dataframe.loc[object,'label2']
    obj_const  = objects_dataframe.loc[object,'constellation']
    obj_season = objects_dataframe.loc[object,'season']
    obj_notes  = objects_dataframe.loc[object,'notes']
    score_sum = obj_score1 + obj_score2 + obj_score3 + obj_score4 + obj_score5

    if output_RA == 'hms':
        obj_RA_hours, obj_RA_minutes, obj_RA_seconds, obj_RA_hms = decdeg_to_hms(obj_RA)
        obj_RA = obj_RA_hms

    if output_DEC == 'dms':
        obj_DEC_deg, obj_DEC_minutes, obj_DEC_seconds, obj_DEC_dms = decdeg_to_dms(obj_DEC)
        obj_DEC = obj_DEC_dms

    # Create a dictionary with the data for the new row
    new_row = {'label1': obj_label, 'label2': obj_label2, 'constellation': obj_const, 'magnitude': object_magnitude, 'season': obj_season,
               'size': object_size, 'ra': obj_RA, 'dec': obj_DEC, 'type': obj_type, 'notes': obj_notes, 'score1': obj_score1, 'score2': obj_score2,
               'score3': obj_score3, 'score4': obj_score4, 'score5': obj_score5, 'totalscore': score_sum}

    # Inserting the new row
    filtered_obj_dataframe.loc[len(filtered_obj_dataframe)] = new_row

    # Reset the index
    filtered_obj_dataframe = filtered_obj_dataframe.reset_index(drop=True)

    if DEBUG >= 1:
        log_file.write(f"   SCORES: Mag: {obj_score1:.2f}, Size: {obj_score2:.1f}, Moon angle: {obj_score3:.2f}, Moon full: {obj_score5:.1f}, Hours: {obj_score4}, TOTAL: {score_sum:.5f}\n")
    if DEBUG >= 2:
        print(f"   SCORES: Mag: {obj_score1:.2f}, Size: {obj_score2:.1f}, Moon angle: {obj_score3:.2f}, Moon full: {obj_score5:.1f}, Hours: {obj_score4}, TOTAL: {score_sum:.3f}")

    # ---------------------------------------------------------------------------------------------------------------------------------
    # Reccommend Scope, Camera, Filter, Reducer, Exposure
    #    create a data frame with the object label and the above
    #

    if output_RECS == 'YES':
        # Loop through camera and scope to check the object size against the field of view
        for camera, row in cameras_dataframe.iterrows():
            if str(cameras_dataframe.loc[camera,'imaging']) == "NO":
                continue

            if DEBUG >= 2:
                log_file.write(f"   Recommendation loop - camera: {cameras_dataframe.loc[camera,'label']} and scope:{scopes_dataframe.loc[scope, 'label']}\n")
            if DEBUG >= 3:
                print(f"     Recommendation loop - camera: {cameras_dataframe.loc[camera,'label']} and scope:{scopes_dataframe.loc[scope, 'label']}\n")

            # Loop through scopes
            for scope, row in scopes_dataframe.iterrows():
                # FOV (in arc-min) = 3436 * D / L
                # D is the sensor dimension in mm
                # L is the FL in mm
                vert_sensor_size = 0.001 * cameras_dataframe.loc[camera,'pix_size'] * cameras_dataframe.loc[camera,'vert_pix']
                horz_sensor_size = 0.001 * cameras_dataframe.loc[camera,'pix_size'] * cameras_dataframe.loc[camera,'horz_pix']
                field_of_view =  (3436.0 * math.sqrt(vert_sensor_size**2 + horz_sensor_size**2)) / scopes_dataframe.loc[scope,'focal_length_mm']
                scope_f_stop = scopes_dataframe.loc[scope,'focal_length_mm'] / scopes_dataframe.loc[scope, 'aperture_mm']
                FOV_frac = float(object_size) / float(field_of_view)
                ideal_fov_err1 = abs(FOV_frac - FOV_ideal_frac)

                #print(f"DEBUG: Field_of_view = {field_of_view}")
                #print(f"DEBUG: Field_of_view fraction = {FOV_frac}")
                #print(f"DEBUG: Ideal Field_of_view error 1 = {ideal_fov_err1}")

                # START exposure range
                if obj_type_filter == 'broadband':
                    lp_filter_factor = 3.0
                else:
                    lp_filter_factor = 100.0

                lp_multiplier = (1.0 + (cameras_dataframe.loc[camera, 'QE'] / 100.0 - 0.5)) / lp_filter_factor

                x = [4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0]
                y = [9, 8, 7, 6, 5, 4, 3, 2, 1]

                # bortle  9     8     7     6    5    4    3     2    1   
                z = [[175.0, 62.0, 22.0, 10.8, 5.3, 2.5, 1.2, 0.98, 0.8],  # f4.0
                     [112.0, 40.0, 14.0, 7.2, 3.7, 1.7, 0.81, 0.64, 0.51], # f5.0
                     [78.0,  27.0,  9.6, 5.0, 2.6, 1.2, 0.56, 0.45, 0.36], # f6.0
                     [57.0,  20.0,  7.1, 3.7, 1.9, 0.9, 0.41, 0.33, 0.26], # f7.0
                     [47.0,  17.0,  5.9, 3.0, 1.6, 0.7, 0.34, 0.27, 0.22], # f8.0
                     [38.0,  13.0,  4.6, 2.4, 1.2, 0.6, 0.26, 0.21, 0.17], # f9.0
                     [28.0,  10.0,  3.4, 1.7, 0.85,0.4, 0.19, 0.16, 0.13]] # f10.0

                z = np.array(z).T

                f = interp2d(x, y, z)

                                # fstop, bortle
                interp_lp_value = f(x = scope_f_stop, y = bortle_class) 

                modified_lp_value = interp_lp_value * lp_multiplier
                
                min_exposure_time = float((10.0 * (cameras_dataframe.loc[camera, 'read_noise_l'])**2) / modified_lp_value)
                max_exposure_time = float((10.0 * (cameras_dataframe.loc[camera, 'read_noise_h'])**2) / modified_lp_value)
                # stop here
                # END exposure range

                for jj in range(1, 3):
                    if jj == 1:
                        new_row_rec = {'label1': obj_label, 'total_score': score_sum, 'scope': scopes_dataframe.loc[scope, 'label'], 'camera': cameras_dataframe.loc[camera, 'label'], 
                                       'reducer': 'NO', 'field_of_view': field_of_view, 'fov_frac': FOV_frac, 'fov_err': ideal_fov_err1, 'filter': obj_type_filter,
                                       'min_exp': min_exposure_time, 'max_exp': max_exposure_time}
                    else:
                        field_of_view_ff =  field_of_view / scopes_dataframe.loc[scope, 'reducer_factor']
                        FOV_frac_ff = float(object_size) / float(field_of_view_ff)
                        ideal_fov_err2 = abs(FOV_frac_ff - FOV_ideal_frac)
                        new_row_rec = {'label1': obj_label, 'total_score': score_sum, 'scope': scopes_dataframe.loc[scope, 'label'], 'camera': cameras_dataframe.loc[camera, 'label'], 
                                       'reducer': 'YES', 'field_of_view': field_of_view_ff, 'fov_frac': FOV_frac_ff, 'fov_err': ideal_fov_err2, 'filter': obj_type_filter,
                                       'min_exp': min_exposure_time, 'max_exp': max_exposure_time} 

                    if jj == 2:
                        if ideal_fov_err1 <= ideal_fov_err2:
                            new_row_rec = {'label1': obj_label, 'total_score': score_sum, 'scope': scopes_dataframe.loc[scope, 'label'], 'camera': cameras_dataframe.loc[camera, 'label'], 
                                       'reducer': 'NO', 'obj_size': object_size, 'field_of_view': field_of_view, 'fov_frac': FOV_frac, 'fov_err': ideal_fov_err1, 'filter': obj_type_filter,
                                       'min_exp': min_exposure_time, 'max_exp': max_exposure_time}
                        else:
                            interp_lp_value = f(x = (scope_f_stop * scopes_dataframe.loc[scope, 'reducer_factor']), y = bortle_class) 
                            modified_lp_value = interp_lp_value * lp_multiplier
                
                            min_exposure_time = float((10.0 * (cameras_dataframe.loc[camera, 'read_noise_l'])**2) / modified_lp_value)
                            max_exposure_time = float((10.0 * (cameras_dataframe.loc[camera, 'read_noise_h'])**2) / modified_lp_value)

                            new_row_rec = {'label1': obj_label, 'total_score': score_sum, 'scope': scopes_dataframe.loc[scope, 'label'], 'camera': cameras_dataframe.loc[camera, 'label'], 
                                       'reducer': 'YES', 'obj_size': object_size, 'field_of_view': field_of_view_ff, 'fov_frac': FOV_frac_ff, 'fov_err': ideal_fov_err2, 'filter': obj_type_filter,
                                       'min_exp': min_exposure_time, 'max_exp': max_exposure_time}

                        reccommendations_df.loc[len(reccommendations_df)] = new_row_rec
                        reccommendations_df = reccommendations_df.reset_index(drop=True)

        # Sort the filtered dataframe by total score descending
        sorted_recs_df = reccommendations_df.sort_values(by=['fov_err'], ascending=True)

print(f"\nTOTAL Objects passing filter: {objects_passed}\n")

if objects_passed >= 1:
    # Sort the filtered dataframe by total score descending
    sorted_filtered_df = filtered_obj_dataframe.sort_values(by=['totalscore'], ascending=False)

    #
    # Print up to the first <obj_count> lines of the sorted table
    print(tabulate((sorted_filtered_df[['label1','constellation','magnitude','size','ra','dec','type','score4','totalscore']].head(obj_count)), headers='keys', tablefmt='psql', floatfmt=".3f", showindex=False))

    if output_RECS == 'YES':
        print("\nReccommended Equipment per Object:")
        print(tabulate((sorted_recs_df[['label1','total_score','scope','camera','reducer','obj_size','field_of_view','fov_frac','fov_err','filter','min_exp','max_exp']].head(obj_count)), headers='keys', tablefmt='psql',floatfmt=".3f", showindex=False))