import pandas as pd
import numpy as np
from IPython.core.pylabtools import figsize
from config import Configuration
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('TkAgg')
import numpy as np
from astropy.io import fits
from astropy.stats import sigma_clipped_stats

lc = pd.read_csv(Configuration.LIGHTCURVE_DIRECTORY + "FIELD_36.007/5909.lc", sep=' ')
ref1 = pd.read_csv(Configuration.LIGHTCURVE_DIRECTORY + "FIELD_36.007/0523.lc", sep=' ')
ref2 = pd.read_csv(Configuration.LIGHTCURVE_DIRECTORY + "FIELD_36.007/0448.lc", sep=' ')
ztf = pd.read_csv(Configuration.ANALYSIS_DIRECTORY+ "snr_ztf.txt", sep=' ')

lc['day'] = lc.apply(lambda x: int(x.jd), axis=1)

dys = lc.day.unique()
clns = np.zeros(3)
errs = np.zeros(3)
chks = np.zeros(3)
for idx, dy in enumerate(dys):
    clns[idx] = lc[(lc.day == dy) & (lc.zpt.abs() < 0.2)].cln.mean()
    errs[idx] = lc[(lc.day == dy) & (lc.zpt.abs() < 0.2)].cln.std()
    chks[idx] = ref1[(lc.day == dy) & (lc.zpt.abs() < 0.2)].cln.mean()

plt.errorbar(dys-2400000., clns-(chks - np.mean(chks)) + 4, yerr=errs, fmt='none', c='k')
plt.scatter(dys-2400000., clns-(chks - np.mean(chks)) + 4, c='k')

plt.errorbar(ztf.mjd, ztf.mag, yerr=ztf.er, fmt='none', c='r')
plt.scatter(ztf.mjd, ztf.mag, c='r')
plt.gca().invert_yaxis()
plt.show()


# get the star list
star_list = pd.read_csv(Configuration.MASTER_DIRECTORY + "FIELD_2b.022_star_list.txt", sep=' ')

lc = pd.read_csv(Configuration.LIGHTCURVE_DIRECTORY + "2024-09-24/FIELD_2b.022/2084.lc", sep=' ')

plt.figure(figsize=(6,6))
plt.errorbar((lc.jd - lc.jd.min()) * 24 * 60, lc.cln + 1.2, yerr=lc.er, c='k', fmt='none')
# plt.errorbar([0.18,0.18], [16.3, 16.3], yerr=[np.mean(lc.er), np.mean(lc.er)], c='k', fmt='none')
plt.scatter((lc.jd - lc.jd.min()) * 24 * 60, lc.cln + 1.2, c='k')
plt.title('SX For')
plt.xlabel('Time Since First Exposure [min]')
plt.ylabel('T$_G$')
plt.gca().invert_yaxis()
# plt.show()
plt.savefig(Configuration.ANALYSIS_DIRECTORY + 'toros_commissioning_sx_for.png', bbox_inches='tight', pad_inches=0.1)

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

    lc = pd.read_csv(Configuration.LIGHTCURVE_DIRECTORY + "2024-09-24/FIELD_2b.022/" + star_id + ".lc", sep=' ')

    rms_cln[idx] = lc[lc.cln > 0].cln.std()
    rms_mag[idx] = lc[lc.cln > 0].mag.std()

plt.figure(figsize=(10,5))

plt.scatter(star_list[star_list.phot_g_mean_mag < 19].master_mag + 1.2, rms_cln[star_list.phot_g_mean_mag < 19], marker='.', c='k', label='Data')
# plt.scatter(star_list[star_list.phot_g_mean_mag < 19].master_mag + 1.2, star_list.master_mag_er[star_list.phot_g_mean_mag < 19]/np.sqrt(300), marker='.', c='r', label='Model')
plt.xlabel('T$_{G}$', fontsize=15)
plt.xticks(fontsize=13)
plt.ylabel('rms', fontsize=15)
plt.yticks(fontsize=13)
plt.yscale('log')

plt.show()
plt.savefig(Configuration.ANALYSIS_DIRECTORY + 'fig2.png', bbox_inches='tight', pad_inches=0.5)

plt.close()

master, mhead = fits.getdata(Configuration.MASTER_DIRECTORY + "FIELD_2b.022_master.fits", header=True)
mn_sky, md_sky, st_sky = sigma_clipped_stats(master, sigma=5)

diff, dhead = fits.getdata(Configuration.DIFFERENCED_DIRECTORY + '2024-09-24/FIELD_2b.022/FIELD_2b.022_12000x10600_170_bkcfspad.fits', header=True)
sci, shead = fits.getdata(Configuration.CLEAN_DIRECTORY + '2024-09-24/FIELD_2b.022/FIELD_2b.022_12000x10600_170_bkcfsp.fits', header=True)

noi = np.sqrt((sci-np.mean(sci)) ** 2 + (master-np.mean(master)) ** 2)
ndiff = diff / noi
print(sigma_clipped_stats(ndiff[5280:5280 + 5280, 3960:3960 + 1320]))
plt.figure(figsize=[6,6])
plt.hist(np.ravel(ndiff), bins=50, range=[-5,5], linewidth=3, color='k', histtype='step')
plt.xlim([-3,3])
plt.ylabel('Count', fontsize=15)
plt.yticks(fontsize=13)
plt.xticks(fontsize=13)
plt.xlabel('Normalized Pixel Value', fontsize=15)
plt.savefig("/home/oelkerrj/Desktop/hist.png", bbox='tight', padinches=0.1)
plt.show()

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

