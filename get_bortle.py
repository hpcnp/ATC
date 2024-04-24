import requests
from geopy.geocoders import Nominatim

def get_bortle_class(latitude, longitude):
    """
    Look up the Bortle class of a location based on its latitude and longitude.
    
    Parameters:
    latitude (float): Latitude of the location.
    longitude (float): Longitude of the location.
    
    Returns:
    bortle_class (int): Bortle class (1-9) of the location.
    """
    # Geocode the location to get the address
    geolocator = Nominatim(user_agent="my_app")
    location = geolocator.reverse(f"{latitude}, {longitude}")
    address = location.address
    
    # Use the Dark Sky Atlas API to look up the Bortle class
    api_key = "YOUR_API_KEY_HERE"
    url = f"https://darksky.atlas.de/api/v1/bortle?address={address}&key={api_key}"
    response = requests.get(url)
    data = response.json()
    
    if "bortle" in data:
        bortle_class = data["bortle"]
        return bortle_class
    else:
        return None

# Example usage
latitude = 37.7749
longitude = -122.4194
bortle_class = get_bortle_class(latitude, longitude)

if bortle_class is not None:
    print(f"The Bortle class of the location is: {bortle_class}")
else:
    print("Unable to determine the Bortle class for the given location.")
