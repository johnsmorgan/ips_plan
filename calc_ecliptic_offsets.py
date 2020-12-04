"""
Calculate coordinates of each observation
"""
import argparse
from datetime import datetime
from numpy import radians, sin, cos, arcsin, arctan2, ones
from yaml import safe_load

from astropy import io
from astropy import units as u
from astropy.coordinates import SkyCoord, Longitude, get_sun, GeocentricTrueEcliptic, GCRS
from astropy.table import Column
from astropy.time import Time

def destination(theta, d, phi1=0., lambda1=0.):
    """
    from http://www.movable-type.co.uk/scripts/latlong.html

    Given a starting lat and long (phi1, lambda1), bearing theta and angular distance d

    Returns new lat and long
    """
    phi2 = arcsin(sin(phi1)*cos(d) + cos(phi1)*sin(d)*cos(theta))
    lambda2 = lambda1 + arctan2(sin(theta)*sin(d)*cos(phi1), cos(d)-sin(phi1)*sin(phi2))
    return phi2, lambda2

parser = argparse.ArgumentParser()
parser.add_argument('infile', nargs=1, help='Input yaml file')
args = parser.parse_args()
conf = safe_load(open(args.infile))

noons = io.ascii.read(conf['files']['noons'])

times = Time(noons['local_noon_str'])
start_time = Time(datetime.combine(conf['startDate'], datetime.min.time()))
stop_time = Time(datetime.combine(conf['stopDate'], datetime.min.time()))
good_times = (times > start_time) & (times < stop_time)
times = times[good_times]

out_table = noons[good_times]

sun_3d = get_sun(times)
sun_eq = SkyCoord(sun_3d.ra, sun_3d.dec)

sun_ra = Column(data=sun_eq.ra, name='ra_sun')
sun_dec = Column(data=sun_eq.dec, name='dec_sun')
out_table.add_columns([sun_ra, sun_dec])

sun_ecliptic = sun_eq.geocentrictrueecliptic
#print sun_eq.geocentrictrueecliptic
for target in conf['priority']:
    if conf['fields'][target]['system'] == "Heliocentric":
        l1, phi = conf['fields'][target]['coordinates']
        target_offset = destination(radians(phi), radians(l1))
        print(target_offset)
        target_ecliptic = GeocentricTrueEcliptic(sun_ecliptic.lon + u.rad*target_offset[1], sun_ecliptic.lat + u.rad*target_offset[0])
        #print(target_ecliptic)
        target_eq = target_ecliptic.transform_to(GCRS)
        #print target_eq
        lon = Column(data=target_ecliptic.lon.deg, name='lon_%s' % (target))
        lat = Column(data=target_ecliptic.lat, name='lat_%s' % (target))
        ra = Column(data=target_eq.ra.deg, name='ra_%s' % (target))
        ha = Column(data=Longitude(target_eq.ra - sun_eq.ra, wrap_angle=180*u.deg).deg, name='ha_%s' % (target))
        dec = Column(data=target_eq.dec, name='dec_%s' % (target))
        out_table.add_columns([lon, lat, ra, dec, ha])
    elif conf['fields'][target]['system'] == "Ecliptic":
        ra_coord, dec_coord = conf['fields'][target]['coordinates']
        ra = Column(data=ra_coord*ones(times.shape)*u.deg, name='ra_%s' % (target))
        ha = Column(data=Longitude(ra_coord*ones(times.shape)*u.deg - sun_eq.ra, wrap_angle=180*u.deg).deg, name='ha_%s' % (target))
        dec = Column(data=dec_coord*ones(times.shape)*u.deg, name='dec_%s' % (target))
        out_table.add_columns([ra, dec, ha])
out_table.write(conf['files']['targets'], format='csv', overwrite=True)
