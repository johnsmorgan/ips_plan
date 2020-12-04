import csv, json, argparse
import numpy as np
from yaml import safe_load
from astropy.coordinates import Angle, AltAz, EarthLocation
from astropy.io import ascii
from astropy import units as u
from astropy.time import Time

S = slice(None, None, None)

parser = argparse.ArgumentParser()
parser.add_argument('infile', help='Input yaml file')
args = parser.parse_args()
conf = safe_load(open(args.infile))

SUN_OBS_STR = "schedule_observation.py --starttime={pre_time_comma} --stoptime=++16s --freq='{coarse_channels}' --obsname={obs_name_prefix}Sun --source=Sun --mode=HW_LFILES --inttime={inttime} --freqres={freqres} --creator=jmorgan --project={project}"
OBSERVATION_STR = "schedule_observation.py --starttime={time_comma} --stoptime=++{duration}s --freq='{coarse_channels}'  --obsname={obs_name_prefix}{field} --shifttime={shifttime} --mode=HW_LFILES --inttime={inttime} --freqres={freqres} --creator={creator} --project={project} --azimuth={az} --elevation={el} --nomode"

NO_WRITE = []

azel = {int(k): v for k, v in json.load(open("azel.json")).items()}

obs_ha = ascii.read(conf['files']['observations'])
noons = Time(obs_ha['local_noon_str'][S])

observations = []

target_time = 0
sun_target_time = 0
schedule_time = 0

for t in conf['priority']:
    times = noons + Angle(np.nan_to_num(obs_ha['ha_%s' % t][S])*u.deg).cycle*u.sday - u.second*conf['obs']['duration']/2.
    # round to nearest 8s (MWA observations must start and stop on these boundaries)
    dt = (86400*times.jd2 % 8)*u.second
    times = np.where(dt<4*u.s, times-dt, times+(8*u.s-dt))
    for j in range(len(noons)):
        if np.isnan(obs_ha['ha_%s' % t][S][j]):
            continue
        print(t)
        out_dict = []
        out_dict = conf['obs']
        out_dict['sweetspot'] = obs_ha['beam_%s' % t][S].data[j]
        out_dict['time'] = times[j].isot[:19]
        out_dict['time_comma'] = times[j].isot[:19].replace('T', ',')
        out_dict['pre_time_comma'] = (times[j]-32*u.second).isot[:19].replace('T', ',')
        out_dict['az'], out_dict['el'] = azel[out_dict['sweetspot']]
        out_dict['obs_name_prefix'] = conf['obsName']
        out_dict['field'] = t
        out_dict['project'] = conf['project']
        out_dict['ha'] = obs_ha['ha_%s' % t][S].data[j]
        out_dict['sun_attenuation'] = obs_ha['sun_attenuation_%s' % t][S].data[j]
        out_dict['target_sensitivity'] = obs_ha['ha_%s' % t][S].data[j]
        observations.append(out_dict.copy())

observations = sorted(observations, key=lambda o: o['time'])

#with open('observations2.json', 'w') as jsonfile:
    #json.dump(observations, jsonfile)

t1 = None
with open(conf['files']['schedule'], 'w') as outfile:
    for o in observations:
        if o['field'] in NO_WRITE:
            t1 = Time(o['time'])
            continue
        t2 = Time(o['time'])
        if t1:
            #print('# ' + str((t2-t1).sec))
            if t2-t1 < u.second*(conf['obs']['duration']+32):
                if t2-t1 < u.second*(conf['obs']['duration']):
                    print("# ERROR Clashing observations %fs apart" % ((t2-t1).sec), file=outfile)
                print("# Observations %.0fs apart -- skipping Sun observation" % ((t2-t1).sec), file=outfile)
            else:
                print(SUN_OBS_STR.format(**o), file=outfile)
                sun_target_time += 16
                schedule_time += 32

        else:
            print(SUN_OBS_STR.format(**o), file=outfile)
            #pass
        #print(o)
        print(OBSERVATION_STR.format(**o), file=outfile)
        target_time += conf['obs']['duration']
        sun_target_time += conf['obs']['duration']
        schedule_time += conf['obs']['duration']
        t1 = t2

    print("# target time: %d sun/target time: %d schedule time: %d" % (target_time, sun_target_time, schedule_time), file=outfile)
with open('observations_long.csv', 'w') as csvfile:
    fieldnames = ('time', 'field', 'sweetspot', 'ha', 'az', 'el', 'sun_attenuation', 'target_sensitivity')
    writer = csv.DictWriter(csvfile, fieldnames=fieldnames, extrasaction='ignore')
    writer.writeheader()
    for out_dict in observations:
        writer.writerow(out_dict)
