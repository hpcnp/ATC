# Import datetime module
import datetime
from astronomica import *

def get_astronomical_night():
    # Get the local sunset time as a datetime object
    #
    latitude = 33.697472      # Latitude    
    longitude = -117.72314     # Longitude

    sun = Sun(latitude, longitude)
    sunset_time = sun.get_local_sunset_time(datetime.datetime.now())
    sunset_time_str = sunset_time.strftime("%H:%M")
    print(f"\n\nThe local sunset time is {sunset_time_str}")

    # Calculate duration of the day in hours and minutes
    #day_duration = datetime.timedelta(hours=24) - datetime.timedelta(hours=sunset_time.hour, minutes=sunset_time.minute)
    sunrise_time = sun.get_local_sunrise_time(datetime.datetime.now())
    day_duration =  datetime.timedelta(hours=sunset_time.hour, minutes=sunset_time.minute) - datetime.timedelta(hours=sunrise_time.hour, minutes=sunrise_time.minute)
    #print(f" DEBUG: day_duration: {day_duration}")

    day_duration_hours, day_duration_minutes = divmod(day_duration.seconds // 60, 60)


    # Calculate the duration of astronomical twilight as a timedelta object
    twilight_duration = datetime.timedelta(minutes=day_duration.total_seconds() * 0.2666 / 60)

    # Calculate the duration of astronomical night as a timedelta object
    night_duration =  day_duration - twilight_duration

    # Calculate the time of astronomical night as a datetime object
    night_time = sunset_time + (twilight_duration / 2.4)

    # Format the time of astronomical night as a string in the format "HH:MM"
    night_time_str = night_time.strftime("%H:%M")

    print(f"The start time of astronomical night tonight is {night_time_str}")

    tomorrow_date = datetime.date.today() + datetime.timedelta(days=1)
    sunrise_time = sun.get_local_sunrise_time(tomorrow_date)
    sunrise_time_str = sunrise_time.strftime("%H:%M")

    dawn_time = sunrise_time - (twilight_duration / 2.4)
    dawn_time_str = dawn_time.strftime("%H:%M")
    
    print(f"\nTotal night duration is: {night_duration}")
    print(f"Twilight duration is: {twilight_duration}")

    print(f"\nThe time of astronomical night ends tomrrow is {dawn_time_str}")
    print(f"The local sunrise time tomorrow is {sunrise_time_str}")

get_astronomical_night()