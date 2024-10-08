import pandas as pd
from config import Configuration
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('TkAgg')
import numpy as np
from astropy.io import fits
from astropy.stats import sigma_clipped_stats
# get the star list
star_list = pd.read_csv(Configuration.MASTER_DIRECTORY + Configuration.FIELD + '_star_list.txt',
                        delimiter=' ', header=0)
master = fits.getdata(Configuration.MASTER_DIRECTORY + Configuration.FIELD + "_master.fits")
mn_sky, md_sky, st_sky = sigma_clipped_stats(master, sigma=5)

lc = pd.read_csv(Configuration.LIGHTCURVE_DIRECTORY + Configuration.DATE + "/" + Configuration.FIELD + "/2084.lc", sep=' ')

plt.figure(figsize=(8,6))
plt.errorbar(lc.jd - 2460579, lc.cln + 2.5*np.log10(300), yerr=lc.er, c='k', fmt='none')
plt.errorbar([0.18,0.18], [16.3, 16.3], yerr=[np.mean(lc.er), np.mean(lc.er)], c='k', fmt='none')
plt.scatter(lc.jd - 2460579, lc.cln + 2.5*np.log10(300), c='k')
plt.title('SX For')
plt.xlabel('JD - 2460579 [days]')
plt.ylabel('T$_V$')
plt.gca().invert_yaxis()
plt.savefig(Configuration.ANALYSIS_DIRECTORY + 'toros_commissioning_sx_for.png', bbox_inches='tight', pad_inches=0)
nstars = len(star_list)
rms_cln = np.zeros(nstars)
rms_mag = np.zeros(nstars)
for idx in range(0, nstars):
    if idx >= 1000:
        star_id = str(idx)
    elif (idx < 1000) & (idx >= 100):
        star_id = '0' + str(idx)
    elif (idx < 100) & (idx >= 10):
        star_id = '00' + str(idx)
    else:
        star_id = '000' + str(idx)

    lc = pd.read_csv(Configuration.LIGHTCURVE_DIRECTORY + Configuration.DATE + "/" + Configuration.FIELD + "/" + star_id + ".lc", sep=' ')

    rms_cln[idx] = lc[lc.cln > 0].cln.std()
    rms_mag[idx] = lc[lc.cln > 0].mag.std()

plt.figure(figsize=(8,6))

plt.subplot(2, 1, 1)
plt.scatter(star_list.phot_g_mean_mag, rms_cln, marker='.', c='k', label='Data')
plt.xlabel('G')
plt.ylabel('rms')
plt.yscale('log')

plt.subplot(2, 1, 2)
plt.scatter(star_list.phot_g_mean_mag, rms_cln / star_list.master_mag_er, marker='.', c='k')
plt.plot(star_list.phot_g_mean_mag, np.ones(len(star_list)), c='r', linewidth=2)
plt.xlabel('G')
plt.ylabel('rms / model')
plt.yscale('log')

plt.savefig(Configuration.ANALYSIS_DIRECTORY + 'toros_commissioning_rms.png', bbox_inches='tight', pad_inches=0)

plt.close()

plt.figure(figsize=(8,6))

xx = star_list.phot_bp_mean_mag.to_numpy() - star_list.phot_rp_mean_mag.to_numpy()
yy = star_list.phot_g_mean_mag.to_numpy() - (star_list.master_mag.to_numpy() + 2.5*np.log10(300))

plt.subplot(2, 1, 1)
plt.scatter(star_list.phot_g_mean_mag, star_list.phot_g_mean_mag - (star_list.master_mag + 2.5*np.log10(300)),
            marker='.', c='k')
zpt = np.median(yy[(~np.isnan(yy)) & (star_list.phot_g_mean_mag < 18)])
plt.plot(np.arange(6, 24), np.zeros(len(np.arange(6, 24))) + zpt, c='r', label='Z$_{pnt}$='+str(np.around(zpt, decimals=2)))
plt.xlabel('G')
plt.ylabel('Zeropoint [G - T$_V$]')
plt.xlim([8, 22])
plt.ylim([-10, -2])
plt.legend()

plt.subplot(2, 1, 2)

p = np.polyfit(xx[~np.isnan(xx) & ~np.isnan(yy)], yy[~np.isnan(xx) & ~np.isnan(yy)], 1)
v = np.poly1d(p)
plt.scatter(star_list.phot_bp_mean_mag - star_list.phot_rp_mean_mag, star_list.phot_g_mean_mag - (star_list.master_mag + 2.5*np.log10(300)),
            marker='.', c='k')
plt.plot(np.arange(-1, 5), v(np.arange(-1, 5)), c='r', label='Z$_{pnt}$='+str(np.around(p[0], decimals=2))+'x + '+str(np.around(p[1], decimals=2)))
plt.xlabel('G$_B$ - G$_R$')
plt.ylabel('Zeropoint [G - T$_V$]')
plt.ylim([-10, -2])
plt.xlim([0, 3.])
plt.legend()
plt.savefig(Configuration.ANALYSIS_DIRECTORY + 'toros_commissioning_zeropoint.png', bbox_inches='tight', pad_inches=0)

plt.close()

actual_sky = 25 - 2.5 * np.log10(md_sky / 300 / Configuration.PIXEL_SIZE) + zpt

print('The sky is ' + str(np.around(actual_sky, decimals=2)) + ' mag/sq arcsec.')