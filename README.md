## Instructions
run `get_local_noons` (this should only need to be done once per year. It's very slow.

run get_offset_points.py to convert from polar coordinates (solar elongation and angle) to offset from the Sun in ecliptic latitude and longitude.

run calc_ecliptic_offsets.py to get position in ecliptic and equatorial coordinates for each pointing for each day OUT:targets.csv

If you need to include different types of targets, these can be combined together using topcat or similar.

Next run schedule_interpolated.py (target?.csv is an input). This produces obs_ha*.csv
Note this sets the priority to Boss, 2, 1, 3, 5, 4, 6. (to which we will add 7,8,9)
However, the actual numbering in the output file doesn't change. 

Next run generate_observations2.py (obs_ha*.csv is an input). This produces ips_*.sh with single_observation.py commands.

Notes for 190405
================


We scheduled 60, 90, 120 (and same on other side) up to 1st april inclusive (targets_boss.csv) and 90, 120, 150 from 2nd April (onwards targets_boss2.csv)

Therefore we only need to update targets_boss2.csv. Current observations should be given the same priority, then we can add 60 and 300 at lower priority, and 180 at lowest priority.

Then we need a targets3.csv. There will be no boss field to worry about.

Observations2
-------------

NB Error in current schedule!!! (only applies to targets_boss2.csv) 

last three pointings are:
[-16.1, -25.7],
[-30.0,   0.0],
[-26.6,  14.5],
should be
[-16.1, -25.7],
[-30.0,   0.0],
[-26.6, -14.5],
In other words, instead of scheduling SW, W and SSW, we have SSW, W and NW!

Won't be an issue once we schedule the new ones, except that they will be incorrectly labeled.

Actually, the labeling is all wrong anyway:

TARGETS = ('1', '2', '3', '4', '5', '6', 'boss')
POINTINGS = ['SE', 'E', 'SSE', 'SW', 'W', 'SSW', 'BOSS']

Actually this should have been
POINTINGS = ['SE', 'E', 'SSE', 'SSW', 'W', 'NW', 'BOSS']

TARGETS = ('1', '2', '3', '4', '5', '6', 'boss', '7', '8', '9')
So now it will be
POINTINGS = ['SE', 'E', 'SSE', 'SSW', 'W', 'NW', 'BOSS', 'SW', 'NE', S]

Finally, run generation_observations2.py and simply grep out those observations that are already scheduled.
- might be hard to grep out the Sun ones

- just modify generate_observations2.py not to print them out.

* check if we clash with solar observations!
  - looks like we just need to exclude observations within 1 hour of local solar noon.

  - also some low-luminosity QSOs in week starting June 16th.

Observations3
-------------

TARGETS = ('1', '2', '3', '4', '5', '6', '7', '8', '9')
So now it will be
POINTINGS = ['SE', 'E', 'SSE', 'NE', 'SW', 'W', 'SSW', 'NW', S]
Priority = ['1', '5', '2', '6', '3', '7', '4', '8', '9']
