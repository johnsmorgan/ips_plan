
## Instructions
This assumes that you have a beam file in the correct format (ha vs dec not az vs el. as used for lookup_beam)

Run `get_local_noons [yaml_file]` (this should only need to be done once per year. It's very slow.

Run `calc_ecliptic_offsets.py [yaml_file]` to get position in ecliptic and equatorial coordinates for each pointing for each day OUT:targets.csv

Run `schedule_interpolated.py [yaml_file]` (target?.csv is an input). This produces obs_ha*.csv

Run `generate_observations.py [yaml_file]` (obs_ha*.csv is an input). This produces ips_*.sh with single_observation.py commands.

## Generating the beams

Currently this is done in a rather ad-hoc way: first, all of the beams (except those with rotational symmetry) are generated using the scripts written for the IPS5 paper. Then these are re-gridded into HA vs Decl. We should redo this so we can use lookup_beam.py directly.
