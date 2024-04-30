import pandas as pd
import math

#print("SCOPES")
#scopes_dataframe = pd.read_csv('/Users/christopherporter/Desktop/ASTRO/ATC/SCOPES.dat')
#print(scopes_dataframe)
#print("\n\n")
#
#print("Edge HD 925 FL:",scopes_dataframe.loc[0,'focal_length_mm'])


# =========
#
# READ INPUT FILES ===================================================================================================================
#
print("SCOPES")
scopes_dataframe = pd.read_csv('/Users/christopherporter/Desktop/ASTRO/ATC/SCOPES.dat')
print(scopes_dataframe)
print("\n\n")

# Scope Fields
#  1- Label, 
#  2- Focal Length (mm), 
#  3- aperture (mm), 
#  4- Reducer (factor), 
#  5- Reducer On Scope (YES/NO)

print("CAMERAS")
cameras_dataframe = pd.read_csv('/Users/christopherporter/Desktop/ASTRO/ATC/CAMERAS.dat')
print(cameras_dataframe)
print("\n\n")

# Camera Fields
#  1- Label, 
#  2- Vertical Pixels, 
#  3- Horizontal Pixels, 
#  4- Pixel Size (Microns), 
#  5- Type (Mono/OSC), 
#  6- QE (%), 
#  7- Cooled (YES,NO)
#  8- Imaging (YES,NO)

print(" ..Reading objects DB")
objects_dataframe = pd.read_csv('/Users/christopherporter/Desktop/ASTRO/ATC/astro_db.dat')
print(" ..Done Reading objects DB\n\n")

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

#
# Calculate Filed of View for each scope and camera combination
# Loop through cameras
for camera, row in cameras_dataframe.iterrows():
    if str(cameras_dataframe.loc[camera,'imaging']) == "NO":
        print("  DEBUG: Skipping, not an imaging camera")
        continue
    # Loop through scopes
    for scope, row in scopes_dataframe.iterrows():
        # FOV (in arc-min) = 3436 * D / L
        # D is the sensor dimension in mm
        # L is the FL in mm
        vert_sensor_size = 0.001 * cameras_dataframe.loc[camera,'pix_size'] * cameras_dataframe.loc[camera,'vert_pix']
        horz_sensor_size = 0.001 * cameras_dataframe.loc[camera,'pix_size'] * cameras_dataframe.loc[camera,'horz_pix']
        field_of_view =  (3436.0 * math.sqrt(vert_sensor_size**2 + horz_sensor_size**2)) / scopes_dataframe.loc[scope,'focal_length_mm']
        print(f"   DEBUG: Camera: {cameras_dataframe.loc[camera,'label']},  Scope: {scopes_dataframe.loc[scope,'label']}, field of view (arcmin): {field_of_view:.3f}")

        # Arcsec per pixel =(PS/FL)∗206.265
        # PS is the pixel size (microns)
        # FL is the focal length (mm)
        arcsec_per_pixel = (206.265 * cameras_dataframe.loc[camera,['pix_size']]) / scopes_dataframe.loc[scope,'focal_length_mm']
        print(f"   DEBUG: Camera: {cameras_dataframe.loc[camera,'label']},  Scope: {scopes_dataframe.loc[scope,'label']}, arcsec per pixel: {arcsec_per_pixel}")

        #tmp_data = {cameras_dataframe.loc[camera,0]: [scopes_dataframe.loc[scope,0], field_of_view, arcsec_per_pixel]}