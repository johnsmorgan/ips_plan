from yaml import safe_load
from datetime import datetime
from astropy.time import Time
from astropy.coordinates import EarthLocation, Longitude
from astroplan import Observer
import astropy.units as u

INTERVAL=1 # days
conf = safe_load(open("MWA_IPS_2020A_full.yaml"))

location = EarthLocation.from_geodetic(lat=conf['lat']*u.deg,
                                       lon=conf['lon']*u.deg,
                                       height=conf['alt']*u.m)
obs = Observer(location)

t = Time(datetime.combine(conf['startDate'], datetime.min.time()))  + 12*u.hour - conf['timezone']*u.hour
t_end = Time(datetime.combine(conf['stopDate'], datetime.min.time()))  + 12*u.hour - conf['timezone']*u.hour

with open(conf['files']['noon'], 'w') as f:
    print("local_noon_str,local_noon_lst", file=f)
    while t <= t_end:
        local_noon_str = obs.noon(t, which='nearest').isot
        local_noon_lst = Longitude(obs.local_sidereal_time(local_noon_str), wrap_angle=180*u.deg).deg
        print(local_noon_str, end=' ')
        print(local_noon_lst)
        print()
        print("%s,%.3f" % (local_noon_str[:19],local_noon_lst), file=f)
        t += INTERVAL*u.day