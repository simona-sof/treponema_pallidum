import pandas as pd
from geopy.geocoders import Nominatim

def get_lat_long(location_name):
    geolocator = Nominatim(user_agent="geoapiExercises")
    location = geolocator.geocode(location_name)
    if location:
        return location.latitude, location.longitude
    else:
        return None, None
    
def long_lat(df):
    latitudes = []
    longitudes = []
    for index, row in df.iterrows():
        location_name = row.iloc[0]
        if location_name == 'taï national park (africa)':
            location_name = 'taï national park'
        lat, long = get_lat_long(location_name)
        latitudes.append(lat)
        longitudes.append(long)
    return longitudes, latitudes

def main():
    # PATH to csv file
    df = pd.read_csv("pre-processing/meta/meta.csv")
    # Sort by accessions
    df = df.sort_values(by=['accession'])
    # adjust date to format XXXX-XX-XX
    df['date'] = df['date'].str[:10]
    # PATH to save new meta
    df.to_csv('pre-processing/meta/sorted_meta.csv', index=False)

    # Add location coordinates
    locations = pd.DataFrame(df['location'].unique())
    long, lat = long_lat(locations)
    locations['latitude'] = lat
    locations['longitude'] = long
    # PATH to save lat_longs
    locations.to_csv("config/lat_longs.csv", index=False, header=False)

if __name__ == "__main__":
    main()
