import datetime
from dateutil import tz
from suntime import Sun, SunTimeException

warsaw_lat = 51.21
warsaw_lon = 21.01

latitude = 33.697472      # Latitude    
longitude = -117.72314     # Longitude

sun = Sun(latitude, longitude)

# Get today's sunrise and sunset in UTC
today_sr = sun.get_sunrise_time()
today_ss = sun.get_sunset_time()
print('Today at HOME the sun rose at {} and went down at {} UTC'.
      format(today_sr.strftime('%H:%M'), today_ss.strftime('%H:%M')))

# On a special date in your machine's local time zone
abd = datetime.date(2024, 4, 26)
abd_sr = sun.get_sunrise_time(abd, tz.gettz('America/Los_Angeles'))
abd_ss = sun.get_sunset_time(abd, tz.gettz('America/Los_Angeles'))
print('On {} the sun in Los Angeles rose at {} and went down at {}.'.
      format(abd, abd_sr.strftime('%H:%M'), abd_ss.strftime('%H:%M')))

# Error handling (no sunset or sunrise on given location)
latitude = 87.55
longitude = 0.1
sun = Sun(latitude, longitude)
try:
    abd_sr = sun.get_sunrise_time(abd)
    abd_ss = sun.get_sunset_time(abd)
    print('On {} at somewhere in the north the sun raised at {} and get down at {}.'.
          format(abd, abd_sr.strftime('%H:%M'), abd_ss.strftime('%H:%M')))
except SunTimeException as e:
    print("Error: {0}.".format(e))