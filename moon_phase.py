import math

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
    # Calculate the Julian Day Number
    jdn = int(365.25 * (year + 4716)) + int(30.6001 * (month + 1)) + day + (hour / 24) - 1524.5
    
    # Calculate the moon's age in days
    moon_age = (jdn - 2451549.5) / 29.530588853
    
    # Determine the moon phase
    if moon_age < 1:
        phase = "New Moon"
    elif moon_age < 7.38:
        phase = "Waxing Crescent"
    elif moon_age < 7.38 + 7.38:
        phase = "First Quarter"
    elif moon_age < 14.77:
        phase = "Waxing Gibbous"
    elif moon_age < 14.77 + 7.38:
        phase = "Full Moon"
    elif moon_age < 22.15:
        phase = "Waning Gibbous"
    elif moon_age < 22.15 + 7.38:
        phase = "Last Quarter"
    else:
        phase = "Waning Crescent"
    
    return phase

# Example usage
year = 2024
month = 4
day = 16
hour = 10.5
latitude = 40.73  # New York City
longitude = -74.00  # New York City

moon_phase = calculate_moon_phase(year, month, day, hour, latitude, longitude)
print(f"On {month}/{day}/{year} at {hour:02.2f} hours, the moon will be in the {moon_phase} phase.")
