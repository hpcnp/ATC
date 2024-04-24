# README for ATC
Astronomy Target Chooser (ATC)

A tool to help amateur astronomers optimize their photography time by using an astronomical database
and the astronomer's location, setup, and interests to recommend targets, gain, exposure, filters, etc.

INPUTS
------

- Date & Time
- Latitude & Longitude
- DB of astronomical targets
- FOV (scope and camera combination)
- Priority list of targets
- Override list of targets
- Priority list of target types
- Minimum visible time
- Minimum elevation 
- Minimum object size (or fraction of the FOV)
- Integration time by apparent magnitude and object type
- Borttle class for location 
- List of cameras, sensor sizes, pixel size, QE, and mm/mc
- Plan output type (asiair, Nina, kstars)

OUTPUTS
-------

- Nightly target list
- Recommended camera, scope combo
- Predicted temperature and dew point
- Lens heater power
- Exposure count
  - Gain
  - Filter
  - Temperature target for camera
- Log

CONSTRAINTS
-----------

- Angular distance from moon
- Moon phase
- Astro night start
- Astro night end
- Light pollution no-go zone
- Target type priority 
- Storage availability (GB)


ASSUMPTIONS
-----------

- The mount is polar aligned
- Cloud and weather detection system (which can override everything and shut down)
- No need to take darks, flats, bias frames
- Mount starts from home position
- Automation has the ability to start / stop imaging session, including turn things on and off

