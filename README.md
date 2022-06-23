## Instructions
This assumes that you have a beam file in the correct format (ha vs dec not az vs el. as used for lookup_beam)

All of the metadata required to schedule the entire observing run (as well as the locations of generic data files) are all recorded in a `yaml`-format configuration file. See the example in the codebase.

Each command reads from a single file and produces a single output file. These files are all specified in the configuration file.

Run `get_local_noons [yaml_file]` - this determines the UTC time (and Local Sidereal time) of local noon for each day. All scheduling is done relative to this LST.

Run `calc_ecliptic_offsets.py [yaml_file]` to get position in ecliptic and equatorial coordinates for each pointing for each day

Run `schedule_interpolated.py [yaml_file]` This is where the actual optimised scheduling happens.

Run `generate_observations.py [yaml_file]` (obs_ha*.csv is an input). This produces ips_*.sh with the commands that the MWA operations team can use to schedule the observations.
