import sys
from yaml import safe_load
from h5py import File
import csv
import numpy as np
#from matplotlib import pyplot as plt
from scipy.interpolate import interp1d
#PLOT=False
INTERVAL_DEGREES = 3.3424596978682920e-02 # interval in degrees for fine HA search
N=75 #observation length in units of INTERVAL_DEGREES
HA_TWEAK = 2 # degrees to allow the hour angle to drift in fine optimisation for solar nulling
conf = safe_load(open("MWA_IPS_2020B.yaml"))

COLS = conf['priority']
#LABELS = {k: v for k, v in zip(('1', '2', '3', '4', '5', '6', 'boss', '7', '8', '9'), ('SE', 'E', 'SSE', 'SSW', 'W', 'NW', 'BOSS', 'SW', 'NE', 'S'))}


targets = csv.DictReader(line for line in open(conf['files']['targets']))
df = File("/home/jmorgan/Projects/allsky/sweetspots/beams.hdf5", 'r')
fine_scale = np.arange(np.min(df['beams'].dims[2][0]), np.max(df['beams'].dims[2][0]), INTERVAL_DEGREES)

obs_ha = []

arg_closest = lambda x, y: np.argmin(np.abs((x-y)))


def neighbours(arr, val):
    """
    return two closest values in arr to val
    assumes arr is sorted (lowest value first)
    """
    closest = arg_closest(arr, val)
    if arr[closest] > val:
        closest1 = closest-1
        closest2 = closest
    else:
        closest1 = closest
        closest2 = closest+1
    return closest1, closest2


def lin_interp(y1, y2, dx):
    """
    return linear interpolationn of y1 and y2

    y1 and y2 are y(x1) and y(x2)
    
    dx controls the relative weighting of y1 and y2 in the interpolation
    """
    return y1*(1-dx)+y2*dx

for target in targets:
    print(target['local_noon_str'][:10])
    sun_dec = float(target['dec_sun'])
    out_dict = {}
    for key in ('local_noon_str', 'local_noon_lst'):
        out_dict[key] = target[key]
    sun_dec_idx = neighbours(df['beams'].dims[1][0][...], sun_dec)
    sun_beam = lin_interp(df['beams'][:, sun_dec_idx[0], :],
                          df['beams'][:, sun_dec_idx[1], :],
                          sun_dec-df['beams'].dims[1][0][sun_dec_idx[0]])
    #print(sun_beam)
    ha_grid = None
    day_has = []
    for c in COLS:
        if target['ha_%s' % c] == '':
            out_dict['ha_%s' % c] = np.nan
            out_dict['beam_%s' % c] = -1
            out_dict['sun_attenuation_%s' % c] = np.nan
            out_dict['target_sensitivity_%s' % c] = np.nan
            continue
        daily_ha = []
        target_dec = float(target['dec_%s' % (c)])
        target_dec_idx = neighbours(df['beams'].dims[1][0][...], target_dec)

        #only roll to integer degree
        target_ha_idx = int(round(float(target['ha_%s' % c])))
        # produce sun_beam and target beam
        # both are 2D arrays (pointing, ha)

        target_beam = lin_interp(df['beams'][:, target_dec_idx[0], :],
                                 df['beams'][:, target_dec_idx[1], :],
                                 target_dec - df['beams'].dims[1][0][target_dec_idx[0]])
                                
        target_beam = np.roll(target_beam, target_ha_idx, axis=1)/np.expand_dims(df['broadness'][...], 1)
        sun_filter = sun_beam < 10**conf['solarAttenuationCutoff']

        # applying sun_filter to target_beam will return a ravelled array.
        # we need to be able to identify the original location of our peak sensitivity within target_beam
        if ha_grid is None:
            ha_grid, beam_grid = np.meshgrid(np.arange(target_beam.shape[1]), np.arange(target_beam.shape[0]))
            # also need a filter to stop subsequent observations from clashing
            target_filter = np.ones(target_beam.shape, np.bool)

        try:
            flat_idx = np.nanargmax(target_beam[sun_filter&target_filter])
        except ValueError:
            print("Warning, no observation meets criteria for %s target %s" % (target['local_noon_str'], c))
            out_dict['ha_%s' % c] = np.nan
            out_dict['beam_%s' % c] = -1
            out_dict['sun_attenuation_%s' % c] = np.nan
            out_dict['target_sensitivity_%s' % c] = np.nan
            #continue
            break

        ha_idx = ha_grid[sun_filter&target_filter][flat_idx]
        beam_idx = beam_grid[sun_filter&target_filter][flat_idx]
        #FIXME interpolate to 8s resolution and choose minimum within 4-minute window
        ha =  df['beams'].dims[2][0][ha_idx]
        beam = df['beams'].dims[0][0][beam_idx]
        # refine ha 
        # 
        fine_ha_interp = interp1d(df['beams'].dims[2][0][ha_idx-2:ha_idx+3],
                                  sun_beam[beam_idx, ha_idx-2:ha_idx+3],
                                  kind='quadratic')
        fine_has = fine_scale[np.abs(fine_scale-ha)<2.0]
        #print("ha=%f min=%f max=%f" % (ha, fine_has[0], fine_has[-1]), end=' ')
        fine_sun_beam = fine_ha_interp(fine_has)
        for ha_ in day_has:
            fine_sun_beam = np.where(np.abs(fine_has-ha_)<N*INTERVAL_DEGREES, np.inf, fine_sun_beam)
        #if np.any(fine_sun_beam==np.inf):
            #print()
            #print(fine_sun_beam)
            #print(day_has)
        ha = fine_has[np.argmin(fine_sun_beam)]
        #print(ha)
        flag_range = neighbours(df['beams'].dims[2][0], ha)
        target_filter[:, flag_range[0]-2:flag_range[1]+2] = False
        day_has.append(ha)

        out_dict['beam_%s' % c] = beam
        out_dict['ha_%s' % c] = ha
        out_dict['sun_attenuation_%s' % c] = np.log10(np.min(fine_sun_beam))
        out_dict['target_sensitivity_%s' % c] = target_beam[beam_idx, ha_idx]
#        if c in PLOT:
#            beam = df['beams'][beam_idx]
#            plt.imshow(np.log10(beam[::-1]), extent=(-179.5, 179.5, -90, 90), zorder=1, vmin=-4, vmax=0)
#            plt.plot(ha, sun_dec, 'o', color='yellow', zorder=2)
#            plt.plot(-float(target['ha_%s' % c]) + ha, float(target['dec_%s' % c]), '+', color='black', zorder=2)
#            plt.title("pointing=%s az=%d alt=%d sweetspot=%d\nlog10(sun_attenuation)=%+.1f target_sensitivity=%.3f" % (LABELS[c], df['az'][beam_idx], df['el'][beam_idx], df['sweetspot_number'][beam_idx], np.log10(sun_beam[beam_idx, ha_idx]),target_beam[beam_idx, ha_idx]))
#            plt.xlabel("HA / degrees")
#            plt.ylabel("Decl. / degrees")
#            plt.savefig('%s_%s_hdf5.png' % (target['local_noon_str'][:10], c))
#            plt.close()
    obs_ha.append(out_dict)

with open(conf['files']['observations'], 'w') as csvfile:
    fieldnames = sorted(obs_ha[0].keys(), key = lambda k: k[::-1])
    fieldnames.insert(0, fieldnames.pop(-1))
    fieldnames.insert(0, fieldnames.pop(-1))
    writer = csv.DictWriter(csvfile, fieldnames=fieldnames)

    writer.writeheader()
    for out_dict in obs_ha:
        writer.writerow(out_dict)
