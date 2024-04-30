import suntime
from datetime import datetime, timedelta
import pytz
from geopy.geocoders import Nominatim

def calculate_sunrise_sunset(latitude, longitude, date):
    # Get the next astronomical sunrise and sunset for the given date
    sunrise_utc, sunset_utc = suntime.sun(latitude, longitude, date)
    
    # Convert UTC time to local timezone
    geolocator = Nominatim(user_agent="your_user_agent_string")
    location = geolocator.reverse("{},{}".format(latitude, longitude))
    address = location.raw['address']

    if 'timezone' in address:
        timezone = address['timezone']
        tz = pytz.timezone(timezone)
        local_sunrise = sunrise_utc.astimezone(tz)
        local_sunset = sunset_utc.astimezone(tz)
        timezone_label = "Local time"
    else:
        # If timezone is not available, fallback to UTC
        local_sunrise = sunrise_utc
        local_sunset = sunset_utc
        timezone_label = "UTC"

    print("Timezone: ", timezone_label)
    return local_sunrise, local_sunset


if __name__ == "__main__":
    latitude = 37.7749    # Replace with your latitude
    longitude = -122.4194  # Replace with your longitude
    date = '2025-04-25'    # Replace with your desired date (YYYY-MM-DD format)

    sunrise, sunset = calculate_sunrise_sunset(latitude, longitude, date)

    print("Astronomical Sunrise: ", sunrise)
    print("Astronomical Sunset: ", sunset)